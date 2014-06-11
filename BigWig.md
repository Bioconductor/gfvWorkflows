# BigWig example


```r
library(GenomicFiles)
fls <- list.files("data", "*bigWig", full=TRUE)
exp <- c("wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig", 
         "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig", 
         "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig", 
         "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep2.bigWig")
if (!identical(exp, basename(fls)))
    stop("expected files:\n  ", paste(exp, collapse="\n  "))
```

### Parallelization using BiocParallel


```r
register(SerialParam())
```

### Object construction is lightweight


```r
seqinfo <- seqinfo(BigWigFile(fls[1]))
rng <- tileGenome(seqlengths(BigWigFile(fls[1])), tilewidth=1e5,
                  cut.last.tile.in.chrom=TRUE)
bwv <- BigWigFileViews(fls, fileRange=rng)
print(object.size(bwv),units="Mb")
```

```
## 0.3 Mb
```

### How to reduce by range

First build a 'map' function:


```r
MAP <- function(FILE, RANGE, ...) {
  import(FILE, selection=RANGE, as="NumericList")[[1]]
}
z <- MAP(fls[1], rng[1])
class(z)
```

```
## [1] "numeric"
```

```r
length(z)
```

```
## [1] 100000
```

```r
summary(z)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     0.0     0.0     0.0     8.6     0.0   410.0
```

Now build a reduce function:


```r
REDUCE <- function(MAPPED, ...) {
  m <- simplify2array(MAPPED)
  rs <- rowSums(m)
  Rle(rs)
}
REDUCE(list(1:4,5:8,9:12))
```

```
## numeric-Rle of length 4 with 4 runs
##   Lengths:  1  1  1  1
##   Values : 15 18 21 24
```


```r
res <- reduceByRange(bwv[1:8,], MAP, REDUCE)
res[1:2]
```

```
## [[1]]
## numeric-Rle of length 100000 with 9004 runs
##   Lengths: 10041     1     4     1     5 ...    12    17     7     8     8
##   Values :     0     1     2     3     7 ...     0     6     8    12    13
## 
## [[2]]
## numeric-Rle of length 100000 with 6426 runs
##   Lengths:    26     8     2     7     4 ...    49     1    75     1 22583
##   Values :    13    17    20    14    18 ...     0     1     2     1     0
```

### The built-in summary method


```r
x <- summary(bwv[1:8,])
assay(x)
```

```
##        [,1]   [,2]     [,3]     [,4]
## [1,] 36.978 29.794    4.218    2.433
## [2,]  9.579  7.454    2.670    1.637
## [3,]  7.971  6.353    3.345    1.515
## [4,]  2.054  2.406    5.258    4.434
## [5,]  3.104  2.217    5.248    3.539
## [6,] 62.175 58.487 4220.426 3200.795
## [7,] 12.348 10.571    6.414    3.008
## [8,]  9.941  7.042   11.596    8.243
```

### The built-in coverage method


```r
x <- coverage(bwv[1:8,])
assay(x)[1,1]
```

```
## [[1]]
## RleList of length 1
## $chr1
## numeric-Rle of length 100000 with 6018 runs
##   Lengths: 10537     1     4     2     1 ...     6    12    24     8     8
##   Values :     0     2     6    10    14 ...     2     0     6    10    11
```
