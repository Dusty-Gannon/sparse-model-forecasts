# Sparse modeling for timeseries and forecasting

## Repo organization

This repo is organized like an R package in order to facilitate code documentation. To get started, use `devtools::load_all()` to load all the existing functions into the environment. All user-defined R functions are written in an R script and placed inside the R directory with RoxyGen2 function documentation. This allows others to clone the repo and use the usual R syntax

``` r
?some_function()
```

to render the documentation for the function and learn how to use it.

Other directories and files (not common to all R packages) can be added as needed.

### General flow

The general workflow in the repo is:

1.  Modular functions for generating data and fitting models are written in `.R`, `.stan`, or `.jags` files and stored in the `/R/`, `/Stan/`, or `/JAGS/` directories (respectively).

2.  For more complicated simulations, the modular functions are used to generate / simulate data that are then stored in an untracked data directory on Beartooth. The data directory is untracked due to the size of some data files. These wrapper scripts are stored in `/Simulations/`.

3.  Code used to fit sparse models to simulated data and visualize results are stored in `/Analyses/`.

### Adding functions to the 'package'

1.  Write a new function in an R script (preferably a new script, but multiple functions can go inside a single script as well).

2.  RStudio makes it easy to document a function. Placing your cursor inside the function, use the Code dropdown menu, then Insert Roxygen Skeleton. Otherwise, follow the formatting for [Roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html).

3.  Fill out the fields to document the function.

4.  Inside the R console, use `devtools::document()` to add the function to the man pages.

5.  Use `devtools::load_all()` to load all the functions in the 'package'

## Revision 1 notes

### Meeting 5 August 2026

#### Prediction vs forecasting

- Care around challenge 1

  - Hypothesis generation across the methods, not inference on parameters

- Add sims with full model to all challenges

## Motivation

- Motivate the choice of HS priors right away

  - flexibility of the approach to sparse modeling objectives

  - Less computationally demanding than all-subsets regression or similar

  - Still requires some adaptation for uses case in time series modeling

## Next steps (by next meeting 25 Aug 2026)

- [ ] Amy -- Add full model fits to sims

- [ ] Kait and Josh -- Tree ring data methods

- [ ] Lauren and Alice -- Conceptual figure

- [ ] Dusty and Topher -- Framing and motivation

## Progress 25 Aug 2026

- Simulations still in progress

  - Add full model (just lm or arima) to all simulations

- Concept Figure

  - If we include prior maths, include arrows/colors to highlight connections

  - Keeping separate panels give an "easy" and "challenge" avenue, which could be nice

- General organization

  - Maybe title includes "with and without covariates"

  - How do we tie all this together from the beginning? reiterate throughout

  - 

## 

- No transparency change
- Same order for colors
- Shift/no shift
- Drop the 0.1 and 0.5 from cov shift fig
- Rotate shift/no shift text
- Middle plot only
- 

## 

- 

- [ ] Topher will take intro

- [ ] harmonize the use of 'prediction' vs 'forecast' throughout

- [ ] estimates vs predictions

- [ ] Double-check Burnham and Anderson ref and the calculations

- [x] Double-check the expectation in the bias-variance tradeoff

### Bookkeeping

- [ ] Fig ref `fig:fourier-rmse-multilength`

  - [ ] Sec04 line 87

  - 
  
## Notes, coworking session Nov 6

### Section 3

- [ ] Issue with collinearity and relation to BV-tradeoff in 3.0 (currently in 3.1 lines 238 - 244).

- [ ] How collinearity and variable selection interacts with covariate shifts

- [ ] Make phenomenon of selected variable effect inflation clearer

- [ ] Introduce covariate shift methods at line 268, after correlation description

- Move lines 277 - 281 elsewhere (maybe discussion?)

- Cut some of intro text in 3.2, keep line 318 and 322 - 328, but move to 3.0

- Cut almost all of lines 329 - 348, but keep 339 - 342 and move to line 268

### Section 4

- Introduce using endogenous ideas -> AR errors -> methods development

- Move ARMA details to appendix

- 

## Notes coworking Nov 10

- [ ] Section 2.1 needs a concluding/transition sentence or two

- [ ] Section 3 methods needs info on PredInt calculations

- [ ] Coverage up to section 3.0


## Notes 17 Nov

- [ ] Tree figure colors

- [ ] Tree fig Thicker vert lines

- [ ] 

- [ ] annotations for pre- and post-shift

- [ ] Change labels

### Discussion ideas

- Remind about the new stuff

- Plain language recommendations/summary of costs and benefits

- Forecasting is critical for ecology

    - Global change is inherently a future-looking idea
    
    - Covariate shifts are inevitable

- Exogenous vs. endogenous. What can we do and what can't we?

- Tree results in the context of the literature

    - Bring back the ecology
    
    - Really big picture
    
- Bias-variance trade-off circle back

    - For trees as well
    
    - Higher weight on the variance?

- 












