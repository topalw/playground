Practical 2: Hardy Weinberg Equilibrium and multiple testing
================
Jerome Goudet and Bruce Weir
2022-07-06































``` r
### TIPS ###
# We first load all the libraries we will need for this practical. 
# Make sure you have those installed 
# If not then please visit https://www2.unil.ch/popgen/teaching/SISG22/
library(gaston)
library(hierfstat)
library(JGTeach)
library(HardyWeinberg)
```

# HW ![\\chi^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cchi%5E2 "\chi^2") tests

1.  With the
    [pan.txt](https://www2.unil.ch/popgen/teaching/SISGData/pan.txt),
    create a bed object using `hierfstat::ms2bed`, and calculate for
    each locus its inbreeding coefficient; use it to test whether the
    loci are in Hardy Weinberg Equilibrium; plot
    ![-log10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-log10 "-log10")
    of these p-values against their expectations under the null
    hypothesis.

<details>
<summary>
Answer
</summary>

``` r
# F=1-Ho/He
pan <- ms2bed("pan.txt")
inb.coeff <- 1 - pan@snps$hz/(2*pan@p*(1-pan@p)) # He = 2p(1-p)
# nb inds
ni <- dim(pan)[1] # or nrow(pan)
x2 <- ni*inb.coeff^2 
p.val.x2 <- pchisq(x2,df=1,lower=FALSE) # check ?pchisq
nl <- dim(pan)[2] # or ncol(pan)
# theoretical dist p.val under null
theo.pval <- (1:nl) / nl
plot(-log10(theo.pval),-log10(sort(p.val.x2)),
     col="red",cex=0.5,xlab="Theo p val dist",
     ylab="emp p-val dist x2", pch=16)
abline(c(0,1))
```

![](pract2_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

</details>

Is it what you would have expected?

2.  redo the same but for loci with minor allele count of at least 10
    (maf
    ![\\geq 0.01](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgeq%200.01 "\geq 0.01");
    this filtering “rule” is often used in genomic analysis). Identify,
    using e.g. a different color, the loci that are not in HWE according
    to this test on the plot of observed heterozygosity against allele
    frequencies. Do you see a pattern?

<details>
<summary>
Answer
</summary>

``` r
x <- (0:1000) / 1000 # same as seq(0,1,by=0.001)
# indexing on maf >= 1%
maf01 <- which(pan@snps$maf >= 0.01) # indexing of snps with maf >= 0.01
nl <- length(maf01)
theo.pval <- 1:nl / nl
plot(-log10(theo.pval),-log10(sort(p.val.x2[maf01])), 
     col="red",cex=0.5,xlab="Theo p val dist",
     ylab="emp p-val dist x2", pch=16)
abline(c(0,1))
```

![](pract2_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# highlight outliers in het ~ p plot 
plot(pan[,maf01]@p, pan[,maf01]@snps$hz, col="black", pch=16, cex=0.6, ylab = 'Het', xlab='p')
lines(x,2*x*(1-x),col="orange") 
outliers <- which(-log10(p.val.x2[maf01])>4) 
# two step indexing (first maf01 and then those points that are p-val outliers)
points(pan[,maf01][,outliers]@p, pan[,maf01][,outliers]@snps$hz, col="red", pch=16, cex=0.6)
```

![](pract2_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

</details>

-   A rule often stated for the validity of a
    ![\\chi^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cchi%5E2 "\chi^2")-test
    is

> The minimum expected number per cell is 1 (5) \[the proportion of
> cells with expected counts lower than 5 should not exceed 20%\]

Is this rule working here?

<details>
<summary>
Answer
</summary>

``` r
par(mfrow=c(1,2))

# what frequency leads to np^2==1 e.g. p=(1/n)^0.5
# nb inds
ni <- dim(pan)[1]
xi <- 1
mafn1 <- which(pan@snps$maf >= (xi/ni)^.5) 
nl <- length(mafn1)
theo.pval <- 1:nl/nl
plot(-log10(theo.pval),-log10(sort(p.val.x2[mafn1])),
     col="red",cex=0.5,xlab="Theo p val dist", pch = 16,
     ylab="emp p-val dist x2",main=expression(np^2>=1));abline(c(0,1))

#what frequency leads to np^2==5
xi <- 5
mafn5 <- which(pan@snps$maf >= (xi/ni)^.5)
nl <- length(mafn5)
theo.pval <- 1:nl/nl
plot(-log10(theo.pval),-log10(sort(p.val.x2[mafn5])),
     col="red",cex=0.5,xlab="Theo p val dist", pch=16,
     ylab="emp p-val dist x2",main=expression(np^2>=5));abline(c(0,1))
```

![](pract2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
par(mfrow=c(1,1))
```

</details>

# HW exact tests

3.  Now use the function `HardyWeinberg::HWExactStats` to obtain the
    exact p-values for these loci, using first the argument `midp` set
    to FALSE and then set to TRUE.

<details>
<summary>
Answer
</summary>

``` r
# what are we cbinding here? 
hw.ex <- HWExactStats(cbind(pan@snps$N0, pan@snps$N1, pan@snps$N2), midp=FALSE)
hw.mp <- HWExactStats(cbind(pan@snps$N0, pan@snps$N1, pan@snps$N2), midp=TRUE)
nl <- length(hw.mp)

par(mfrow=c(1,2))

plot(-log10(1:nl/nl), -log10(sort(hw.ex)), cex=0.6, pch=16, col='red',
     xlab='theoretical p-val dist', ylab='emp p-val dist exact') ; abline(c(0,1))
plot(-log10(1:nl/nl), -log10(sort(hw.mp)), cex=0.6, pch=16, col='red',
      xlab='theoretical p-val dist', ylab='emp p-val dist exact-mp') ; abline(c(0,1))
```

![](pract2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
par(mfrow=c(1,1))
### TIPS 
# ';' allows you to specify two commands in one line  
```

</details>

4.  using the East Asian samples from the 1000 genome project, plot the
    SNPs heterozygosity against the frequency of the alternate allele
    for chr22:0-20M. Add to the plot a line showing the expected
    heterozygosity under Hardy Weinberg Equilibrium; Discuss the
    resulting figures, compared to what we saw with the simulated data.

<details>
<summary>
Answer
</summary>

``` r
# read vcf (care for path!)
ch22 <- read.VCF("chr22_Mb0_20.recode.vcf.gz") 
```

    ## ped stats and snps stats have been set. 
    ## 'p' has been set. 
    ## 'mu' and 'sigma' have been set.

``` r
# description file
samp.desc.file <- "https://www2.unil.ch/popgen/teaching/SISG18/integrated_call_samples_v3.20130502.ALL.panel"
samp.desc <- read.table(samp.desc.file,header=TRUE)
# subset EAS
EAS <- which(samp.desc$super_pop=="EAS")
plot(ch22[EAS,]@p, ch22[EAS,]@snps$hz, col="black", pch=16, cex=0.6, xlab='p',ylab='het',
     main='Het~p for EAS')
lines(x, 2*x*(1-x), col="orange")
```

![](pract2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

</details>

5.  Test for Hardy Weinberg using the exact mid-pvalue test for all
    these loci.

<details>
<summary>
Answer
</summary>

``` r
hw.mp.EAS <- HWExactStats(cbind(ch22[EAS,]@snps$N0, ch22[EAS,]@snps$N1, ch22[EAS,]@snps$N2), midp=TRUE)
```

</details>

6.  plot the p-values of the previous tests against their expectation
    under the null (Rather than the p-values, -log10 of the p-values is
    more informative). Identify on the plot of heterozygosity against
    allele frequencies the loci not conforming to HWE, and discuss the
    results.

<details>
<summary>
Answer
</summary>

``` r
nl <- length(hw.mp.EAS) # length of of p-values
plot(-log10(1:nl/nl), -log10(sort(hw.mp.EAS)), cex=0.6, pch=16);abline(c(0,1))
```

![](pract2_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# index outliers 
outliers <- which(-log10(hw.mp.EAS) > 6)
plot(ch22[EAS,]@p, ch22[EAS,]@snps$hz, col="black", pch=16, cex=0.6, xlab='p', ylab='Het')
lines(x, 2*x*(1-x), col="orange")
points(ch22[EAS,outliers]@p, ch22[EAS,outliers]@snps$hz, col="red",pch=16,cex=0.6)
```

![](pract2_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->
</details>

# Power to detect HW departure \[optional\]

7.  \[optional\] How likely are we to detect departure from HW if
    ![f=0.125](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%3D0.125 "f=0.125")
    with a sample of
    ![100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100 "100")
    individuals?

<details>
<summary>
Answer
</summary>

``` r
ni <- 100
f <- 0.125
pchisq(qchisq(0.95,df=1), df=1, ncp=ni*f^2, lower=FALSE)
```

    ## [1] 0.239527

``` r
#density of chisq with ncp nf2

x <- seq(0.2,20,0.1)
plot(x, dchisq(x,df=1), type="h", col='#FF000080',lwd =2,
     xlab=expression(chi^2), ylab="probability density") # chisq prob dens
lines(x, dchisq(x,df=1,ncp=ni*f^2), type="h", col='#0000FF80', lwd=2) # with ncp
legend('topright', lty=rep(1,2), col=c('#FF000080','#0000FF80'), lwd=c(2,2),legend=c('no ncp','ncp = nf^2'))
abline(v=qchisq(0.95,df=1)) # 95th centile of chisq dist 
```

![](pract2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
</details>

-   How many individuals would be needed to detect departures from HW
    when
    ![f=0.125](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%3D0.125 "f=0.125")
    with a power of
    ![0.8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.8 "0.8")?

<details>
<summary>
Answer
</summary>

``` r
ns <- 1:10*100 
round( pchisq(qchisq(0.95,df=1), df=1, ncp=ns*f^2, lower=FALSE), digits=3)
```

    ##  [1] 0.240 0.424 0.581 0.705 0.798 0.865 0.911 0.942 0.963 0.977

``` r
ni <- 500
f <- 0.125
pchisq( qchisq(0.95,df=1), df=1, ncp=ni*f^2, lower=FALSE)
```

    ## [1] 0.7981762

``` r
plot(x, dchisq(x,df=1), type="h", col="#FF000080", lwd=2,
     xlab=expression(chi^2), ylab="probability density")
#density of chisq with ncp nf2
lines(x, dchisq(x,df=1,ncp=ni*f^2), type="h", col="#0000FF80",lwd=2) 
legend('topright', lty=rep(1,2), col=c('#FF000080','#0000FF80'),lwd=c(2,2), legend=c('no ncp','ncp = nf^2'))
abline(v=qchisq(0.95,df=1)) # 95th centile of chisq dist 
```

![](pract2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
</details>

END OF PRACTICAL 2.