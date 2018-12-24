# hfufs
## Hybrid Fu's Fs calculator

Fu's Fs is a population genetics statistic that may be useful for detecting population expansion due to appearance of a recent mutation conferring higher fitness.
Calculation of Fu's Fs can involve Stirling numbers of the first kind, which grow very large very fast, leading to potential overflow problems. It also involves a logit transformation which raises possible floating point underflow issues.

hfufs accounts for floating point overflow and underflow issues in the calculation of Fu's Fs, enabling calculation on very large alignments, even those that have extreme positive or negative values of Fu's Fs.

Fu's Fs requires three parameters:
1. n - number of individuals sampled
2. k - number of different alleles/haplotypes observed
3. theta - average pairwise differences between different individuals
These parameters can be calculated from an alignment or other representation of sequence alleles using many different software programs, such as PopGenome (https://CRAN.R-project.org/package=PopGenome).

## Installation (in R)
Make sure your library paths are all set for R, and that you have write access to them.
```
library(devtools)
devtools::install_github("swainechen/hfufs")
```
## Basic usage
```
n <- 100
k <- 30
theta <- 12.345
hfufs(n, k, theta)
# -0.7374915
```

## Under Active Construction
This works now if you can calculate n, k, and theta. There are quite a few other packages to do this from an alignment. I'm in the process of making this package and set of functions easier to use and integrate with one of these packages (PopGenome), so please check back often!
In addition, the Stirling number estimator is done and should work well for non-population genetics applications also!
