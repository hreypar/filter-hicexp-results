# filter-hicexp-results

Filter hicexp object that has been both  normalized and compared (difference detection) using `multiHiCCompare`.

## Module description ##

## Input ##

A file called `*.qlf.cycnorm.hicexp.Rds`

## Output ##

- One `*.significantpairs.Rds` file for each input:
- **Plain text files in `.csv` format:**
- **Descriptive Plots:**

### Requirements ###

- [`mk`](http://doc.cat-v.org/bell_labs/mk/mk.pdf "A successor for `make`.")

- [`multiHiCCompare`](https://github.com/dozmorovlab/multiHiCcompare) 

- [`ggplot2`](https://ggplot2.tidyverse.org/)

### Configuration file ###

The file config.mk contains variables that can be set according to the user; these are passed to the `mkfile` and feed the parameters of the functions called in the recipes.

### Author ###

[Helena RG](hreyes@inmegen.edu.mx)

