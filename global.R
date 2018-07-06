## 3/25/2017

## global object loading
## load libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(dbplyr)
library(rlang)
library(ggthemes)
library(tidyr)
library(reshape2)
library(forcats)
library(ggrepel)
library(RSQLite)
library(DT)
source("plotFunctions.R")
source("helpers.R")
# the only global object I really need to load is the expression information and libraries

ig_db <- src_sqlite("./data/ig_db.filt.sqlite3", create = F)
tscan.full <- tbl(ig_db, "tscan.filt")
exp.lm <- tbl(ig_db, "exp.lm")
exp.l2m <- tbl(ig_db, "exp.l2m")
utr <- tbl(ig_db, "utr")

# read in vector that will be used to preserve sample order
sample_order <- readRDS("./data/sample_order.rds")
