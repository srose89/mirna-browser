## helpers.R for immgen mirna browser

## 3/25/2017


# capitalize first letter of string
## to match gene names consistently
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}