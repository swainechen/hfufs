# hfufs
## Hybrid Fu's Fs calculator

Fu's Fs is a population genetics statistic that may be useful for detecting population expansion due to appearance of a recent mutation conferring higher fitness.
Calculation of Fu's Fs can involve Stirling numbers of the first kind, which grow very large very fast, leading to potential overflow problems. It also involves a logit transformation which raises possible floating point underflow issues.

hfufs accounts for floating point overflow and underflow issues in the calculation of Fu's Fs, enabling calculation on very large alignments, even those that have extreme positive or negative values of Fu's Fs.

An additional improvement is to do a single term estimate - hfufs uses a Stirling number estimator for each term of the sum for calculating Fu's Fs. Nico Temme developed a single asymptotic estimator that can be used to directly calculate the appropriate sum, increasing both accuracy and speed. This is implemented in the afufs function and is the recommended way to calculate Fu's Fs.

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
afufs(n, k, theta)
# -0.7368616
```

## More practical usage
The hfufs/afufs functions require you to calculate n, k, and theta from your data yourself. There are quite a few other packages to do this from an alignment. One of these is PopGenome (https://cran.r-project.org/web/packages/PopGenome/index.html). If you have that package installed, you can simply do:
```
library(PopGenome)
library(devtools)
devtools::install_github("swainechen/hfufs")
library(hfufs)
fasta_file <- "/path/to/aligned.fasta"
pg.object <- hf.readData(fasta_file)
pg.dataframe <- hf.alignment.stats(pg.object)
```
This interface is a bit easier for reading in single fasta files. The pg.object variable will still hold all the PopGenome information, and parts will be extracted for convenience into pg.dataframe.

In addition, the Stirling number estimator and the actual hfufs/afufs functions are independent and should work well for non-population genetics applications also!
