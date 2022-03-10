## code to prepare `cell_signatures` dataset goes here

library(tidyverse)
library(tmod)
panglo <- read.delim("PanglaoDB_markers_27_Mar_2020.tsv.gz", sep="\t")

sel <- panglo %>% group_by(cell.type) %>% summarise(n=n()) %>% filter(n >= 10) %>% pull(cell.type)

panglo.df <- panglo %>% filter(cell.type %in% sel) %>%
  filter(grepl("Hs", species)) %>%
  select(gene=official.gene.symbol, Title=cell.type) %>%
  mutate(Title=gsub("cells$", "cell", Title)) %>%
  mutate(source="PanglaoDB") %>% 
  mutate(ID=sprintf("PG%05d", as.numeric(factor(Title))))



## source: http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt
cmark <- read.delim("Human_cell_markers.txt.gz", sep="\t") %>%
  filter(cellType == "Normal cell")

cmark.sel <- cmark %>% filter(cellType == "Normal cell") %>% 
  group_by(cellName) %>%
  summarise(n=n()) %>%
  filter(n >= 10) %>%
  pull(cellName)

cmark.df <- cmark %>% filter(cellName %in% cmark.sel) %>%
  filter(speciesType == "Human") %>%
  select(gene=geneSymbol, Title=cellName) %>%
  mutate(source="CellMarker") %>%
  separate_rows(gene, sep = ", *") %>%
  mutate(gene=gsub("[^[:alnum:]]", "", gene)) %>%
  filter(!duplicated(paste(gene, Title))) %>% 
  mutate(ID=sprintf("CM%05d", as.numeric(factor(Title))))

## source: CIBERSORT, LM22 signature file
csort <- read_tsv("LM22.txt")


## we take top 25 markers for each cell
csort.ctypes <- set_names(colnames(csort)[-1])

csort.markers <- map_dfr(c(10, 25, 50), ~ {
  ntop <- .
  map(csort.ctypes, ~ {
    csort[ order(-csort[[.]]), ][[ "Gene symbol" ]][1:ntop]
  }) %>% imap_dfr(~ {
    data.frame(gene=.x, Title=paste0(.y, " top", ntop), source="CIBERSORT")
  })
}) %>% mutate(ID=sprintf("CS%05d", as.numeric(factor(Title))))

cell_signatures_df <- rbind(cmark.df, panglo.df, csort.markers) 
cell_signatures <- makeTmodFromDataFrame(cell_signatures_df, feature_col = "gene",
                                         module_col = "ID", title_col = "Title",
                                         extra_module_cols="source")
                                         


usethis::use_data(cell_signatures, overwrite = TRUE)
