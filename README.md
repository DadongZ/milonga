# milonga
Multiple imputation for binary surveillance data with dropouts

## Author
Dadong Zhang, Josyf C Mychaleckyj

## Description
Multiple imputation for binary surveillance data with dropouts based on permutations.

## Installation
```r
install.packages('devtools')
devtools::install_github("DadongZ/milonga")

#example
library(milonga)
?milonga::milonga
data(polio)
milonga(polio, cols=2:7)
```

## Questions?
If you have questions or encounter problems when using milonga, please report directly [here](https://github.com/DadongZ/milonga/issues) using the Issues tab by GitHub.
