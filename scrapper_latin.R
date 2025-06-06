library(rvest)
library(httr)
library(jsonlite)
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(proxy)
library(stringi)
library(text2vec)
library(entropy) 
library(udpipe)
setwd("C:\\Users\\wagbr\\OneDrive\\Documentos\\Estilometria")
#-------------------------------------------------
safe_read_html <- function(url) {
  RETRY("GET", url,
        user_agent("Mozilla/5.0 (compatible; R)"),
        times = 4, pause_base = 1, pause_cap = 4) |>
    content(as = "text", encoding = "UTF-8") |>
    read_html()
}

#-------------------------------------------------
get_pages <- function(author_id, work_id) {
  html0  <- safe_read_html(
    sprintf("https://latin.packhum.org/loc/%s/%s/0",
            author_id, work_id))
  
  js_txt  <- html_text2(html_node(html0, "div#dataObj > script"))
  json_raw <- str_match(js_txt,
                        regex("var\\s+locInfo\\s*=\\s*(\\{.*?\\});", dotall = TRUE))[ ,2]
  if (is.na(json_raw))
    stop("locInfo não encontrado — verifique IDs de autor/obra.")
  
  loc <- fromJSON(json_raw, simplifyVector = TRUE)
  
  tibble(idx = seq_along(loc$pages) - 1,
         citation = loc$pages)
}

#-------------------------------------------------
scrape_page <- function(author_id, work_id, idx) {
  pg <- safe_read_html(
    sprintf("https://latin.packhum.org/loc/%s/%s/%s",
            author_id, work_id, idx))
  
  rows  <- html_nodes(pg, "#loctext tr")
  class <- html_attr(rows, "class")     # NA se não houver atributo
  
  tibble(
    kind     = ifelse(str_detect(class, "title"), "title", "text"),
    text_raw = html_text(rows %>% html_node("td:first-child")),
    ref_raw  = html_text(rows %>% html_node("td:nth-child(2)"))
  ) |>                                        # nada filtrado ainda!
    filter(str_squish(text_raw) != "")
}

#-------------------------------------------------
scrape_work <- function(author_id, work_id, polite = TRUE) {
  
  pages <- get_pages(author_id, work_id)
  
  raw_df <- map_dfr(seq_len(nrow(pages)), function(i) {
    if (polite) Sys.sleep(runif(1, .4, .9))
    cat("Página", pages$idx[i], "raspada.\n")
    df <- scrape_page(author_id, work_id, pages$idx[i])
    cbind(pages[i, "citation"], df)
  })
  
  raw_df |>
    separate(citation, into = c("book","chapter","section"),
             sep = "\\.", fill = "right", remove = FALSE) |>
    mutate(
      book    = str_remove(book, "\\(.*\\)"),
      line_no = suppressWarnings(as.integer(str_trim(ref_raw)))
    ) |>
    select(book, chapter, section, line_no, kind, text_raw)
}

#-------------------------------------------------
# EXEMPLO: baixar Annales (autor 1351, obra 5)
annales_df <- scrape_work(1351, 5)
head(annales_df)

#---------------------------------------------------------
# 0. Carregar o modelo udpipe
#---------------------------------------------------------
udpipe_download_model("latin-proiel")
modelo <- udpipe_load_model("latin-proiel-ud-2.5-191206.udpipe")

#---------------------------------------------------------
# 0. Função para remover hifenização de quebra de linha
#---------------------------------------------------------
fix_line_hyphens <- function(txt) {
  str_replace_all(txt, "-\\s+([[:lower:]])", "\\1")
}

#---------------------------------------------------------
# 1. Concatena cada CAPÍTULO (book.chapter)
#---------------------------------------------------------
chapters_df <- annales_df %>%
  filter(is.na(kind)) %>%                              # só linhas de texto
  mutate(ch_id = paste(book, chapter, sep = ".")) %>%     # ex. "15.44"
  group_by(ch_id) %>%
  summarise(text = fix_line_hyphens(paste(text_raw, collapse = " ")),
            .groups = "drop") %>%
  filter(str_count(text, "\\S+") > 100) 
head(chapters_df)



extract_features <- function(ann) {
  
  # ── IDENTIFICADOR DO DOCUMENTO ─────────────────────────
  this_doc <- ann$doc_id[1]      # assume 1 doc_id único por data.frame
  
  # parâmetros configuráveis
  char_ng_max <- 5
  min_denom   <- 5
  
  ## ── syntactic 2-grams (TREE2) ─────────────────────────
  tree_bi <- ann %>%
    filter(!is.na(head_token_id) & upos != "PUNCT") %>%
    mutate(head_upos = ann$upos[match(head_token_id, token_id)],
           feat      = paste0("TREE2_", upos, "_", head_upos)) %>%
    count(feat, name = "n")
  
  ## ── lexical lemma+upos trigrams (T3) ──────────────────
  lemupos <- ann %>%
    filter(upos != "PUNCT") %>%
    transmute(lu = paste0(tolower(lemma), "_", upos))
  
  tok_tri <- tibble(
    feat = zoo::rollapply(lemupos$lu, 3,
                          \(x) paste0("T3_", paste(x, collapse = "_")),
                          by = 1, align = "left", fill = NA)
  ) %>%
    filter(!is.na(feat)) %>%
    count(feat, name = "n")
  
  ## ── token & sentence stats ────────────────────────────
  ntokens  <- nrow(ann %>% filter(upos != "PUNCT"))
  nsents   <- length(unique(ann$sentence_id))
  mean_len <- ntokens / nsents
  
  ## ── lexical richness (TTR, entropy) ───────────────────
  lem_counter <- table(lemupos$lu)
  ttr   <- length(lem_counter) / ntokens
  ent   <- entropy::entropy(lem_counter, unit = "log2")
  
  ## ── char n-grams 1-5 ──────────────────────────────────
  txt_clean <- paste(tolower(ann$token[ann$upos != "PUNCT"]), collapse = " ")
  char_ng_vec <- integer()
  for (n in 1:char_ng_max) {
    grams <- stri_sub(txt_clean,
                      1:(stri_length(txt_clean) - n + 1),
                      length = n)
    char_ng_vec <- c(char_ng_vec,
                     table(paste0("C", n, "_", grams)))
  }
  
  ## ── vetor final de features ───────────────────────────
  feat_vec <- tibble(
    doc_id  = this_doc,                                      # <<<< novo
    feature = c(tree_bi$feat, tok_tri$feat, names(char_ng_vec),
                "TTR_LEMMA", "ENT_LEMMA", "MEAN_SENT_LEN", "N_TOK"),
    value   = c( tree_bi$n / ntokens,
                 tok_tri$n / ntokens,
                 char_ng_vec / ntokens,
                 ttr, ent, mean_len, ntokens)
  ) %>%
    mutate(value = ifelse(grepl("^TREE2|^T3", feature) &
                            ntokens < min_denom, NA, value))
  
  return(feat_vec)
}


# anotação em loop, cada cap. recebe seu próprio doc_id ===================
ann_list <- map2(
  chapters_df$text,
  chapters_df$ch_id,
  ~ udpipe_annotate(modelo, x = .x, doc_id = .y) %>%
    as.data.frame()
)
head(as.data.frame(ann_list[1]))

# extrair features para cada elemento da lista ============================
feats <- map_dfr(ann_list, extract_features)

dtm <- feats %>%
  pivot_wider(names_from = feature,
              values_from = value,
              values_fill = 0) %>%
  column_to_rownames("doc_id") %>%
  as.matrix()
head(dtm[,1:10])
dim(dtm)

# opcional: retire docs/cols totalmente zerados
dtm <- dtm[rowSums(dtm) > 0, colSums(dtm) > 0]
dim(dtm)

linha_testimonium_taciteus <- which(rownames(dtm) == "15.44")

library(ggplot2)
library(ggrepel) 

# ───────────────────────── PCA ───────────────────────────
pca   <- prcomp(dtm, center = TRUE, scale. = TRUE)
coords <- as.data.frame(pca$x[, 1:2])        # PC1 e PC2
coords$chapter   <- rownames(coords)
coords$highlight <- coords$chapter == "15.44"

# ──────────────────────── plot ───────────────────────────
ggplot(coords, aes(PC1, PC2)) +
  # pontos “normais”
  geom_point(data = subset(coords, !highlight),
             color = "grey75", size = 2.5, alpha = .4) +
  # ponto destacado
  geom_point(data = subset(coords, highlight),
             color = "red", size = 5, alpha = .9) +
  # rótulo do ponto destacado
  geom_text_repel(data = subset(coords, highlight),
                  aes(label = chapter),
                  color = "red", size = 5,
                  box.padding = .4) +
  labs(title    = "Annales 15 – PCA dos capítulos (feature set Laken)",
       subtitle = "Capítulo 15.44 realçado em vermelho",
       x = "PC1",
       y = "PC2") +
  theme_minimal(base_size = 14)

# ─────────────────────────────────────────────────────────────
# Agora trabalhamos com a distância de Mahalanobis
# ─────────────────────────────────────────────────────────────

# 2. “Branqueamento”: divide cada coluna pelo desvio-padrão (sdev)
Z <- sweep(pca$x, 2, pca$sdev, FUN = "/")           # docs × (n−1)

# 3. Distância Mahalanobis ao centroide (que vira 0 no espaço Z)
mah <- sqrt(rowSums(Z^2))                           # vetor por capítulo

dist_df <- data.frame(chapter = rownames(dtm),
                      mah     = mah,
                      highlight = rownames(dtm) == "15.44")
dist_df <- dist_df[order(dist_df$mah),]
dist_df$position <- 1:nrow(dist_df)
head(dist_df[order(-dist_df$mah),], 10)

# 4. Visualização – barra + PCA bi-dimensional ------------------
## 4a. Barras
rotular <- c("15.44",
             dist_df %>%
               arrange(desc(mah)) %>%
               slice_head(n = 5) %>%
               pull(chapter))

# ── Construa o gráfico como barras finas e sem eixo x ─────────
ggplot(dist_df, aes(reorder(chapter, mah), mah,
                    fill = highlight)) +
  geom_col(width = .6, show.legend = FALSE,
           alpha = ifelse(dist_df$highlight, .9, .6)) +
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "grey70")) +
  coord_flip() +
  # rótulos apenas nos capítulos escolhidos
  geom_text(aes(label = ifelse(chapter %in% rotular,
                               chapter, "")),
            hjust = -0.1, size = 3.8, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "Annales 15 – Distância de Mahalanobis ao centroide",
       subtitle = "Capítulo 15.44 (vermelho) vs. 5 mais distantes",
       x = NULL, y = "Distância (D)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

save.image("testimonium_taciteum.RData")
