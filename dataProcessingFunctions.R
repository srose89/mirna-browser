## 3/25/2017

## Data processing functions

## Melt Expression data frame
meltExp <- function(exp){
  e.m <- reshape2::melt(as.matrix(exp))
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

# for expression filtering
filterExpressedAQ <- function(exp, cutoff = 0, mod = NULL){
  mod <- if(is.null(mod)) as.factor(gsub("_[1-9]$", "", colnames(exp))) else as.character(mod)
  lt30 <- exp > cutoff
  stable <- table(mod)
  lt30cell <- t(apply(lt30, 1, function(x){ tapply(x, mod, sum) }))
  keep <- apply(lt30cell, 1, function(x){any(x - stable == 0)})
  keep.names <- names(keep[keep == TRUE])
  exp.f <- exp[keep.names,]
  return(exp.f)
}
