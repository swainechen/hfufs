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

## Basic Usage (in R)
```
source("hfufs.R")
n <- 100
k <- 30
theta <- 12.345
hfufs(n, k, theta)
# -0.7374915
```

## Under Active Construction
I'm in the process of converting this to an R package for convenience. If you are reviewing the manuscript for this and need access to the original code, you can find it here (you can wget this):
https://raw.githubusercontent.com/swainechen/hfufs/c247bbc3fc91ad240e5b4fb91fd9256357ca5516/src/hfufs.R
Then the above "Basic Usage" will work.
