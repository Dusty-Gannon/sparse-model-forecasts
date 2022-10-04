
# SParse mOdeling of Non-Gaussian timE Series #

## Repo organization ##

This repo is organized like an R package in order to facilitate code documentation. To get started, use `devtools::load_all()` to load
all the existing functions into the environment.
All user-defined R functions are written in an R script and placed inside the R directory with RoxyGen2 function documentation. 
This allows others to clone the repo and use the usual R syntax
```?some_function()```
to render the documentation for the function and learn how to use it.

Other directories and files (not common to all R packages) can be added as needed.

To add functions to the 'package':

1. Write a new function in an R script (preferably a new script, but multiple functions can go inside a single script as well).

2. RStudio makes it easy to document a function. Placing your cursor inside the function, use the Code dropdown menu, then Insert Roxygen Skeleton. Otherwise, follow the formatting for [Roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html). 

3. Fill out the fields to document the function.

4. Inside the R console, use `devtools::document()` to add the function to the man pages.

5. Use `devtools::load_all()` to load all the functions in the 'package'