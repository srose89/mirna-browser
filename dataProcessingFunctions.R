## 3/25/2017

## Data processing functions

## Melt Expression data frame
meltExp <- function(exp){
  e.m <- melt(as.matrix(exp))
  colnames(e.m) <- c('miRNA', 'Sample', 'expression')
  e.m$cell_type <- gsub("_[1-9]$", "", e.m$Sample)
  return(e.m)
}

# capitalize first letter of string
## to match gene names consistently
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}