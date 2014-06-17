# BigWigFileViews

Access and analyze data from a number of **BigWig** files using the **GenomicFiles** package.

A **BigWig** file is a compressed binary indexed file which contains read
coverage information along the genome. The coverage information is
stored at different resolutions to optimize the speed at which
coverage can be visualized. The paper describing the format can be
found here:

http://bioinformatics.oxfordjournals.org/content/26/17/2204.long

We start by loading the **GenomicFiles** library, and referencing four
files stored in a directory `data`. These are BigWig files storing the
coverage of four RNA-Seq experiments from the ENCODE project. The
files can be downloaded from:

ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/


```r
library(GenomicFiles)
fls <- list.files("data", "*bigWig", full=TRUE)
expect <- c("wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig", 
         "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig", 
         "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig", 
         "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep2.bigWig")
if (!identical(expect, basename(fls)))
    stop("expected files:\n  ", paste(expect, collapse="\n  "))
```

### Parallelization using BiocParallel

The **GenomicFiles** package internally uses **BiocParallel** to allow
parallelization across either *files* or *ranges*. Here, we will just
register a serial back-end. For other settings, see `?register`.


```r
register(SerialParam())
```

### Object construction

We query the first **BigWigFile** for the information about the different
sequences (chromosomes). The functions for querying and importing data
from BigWig files are in the **rtracklayer** Bioconductor package.


```r
( seqinfo <- seqinfo(BigWigFile(fls[1])) )
```

```
## Seqinfo of length 25
## seqnames seqlengths isCircular genome
## chr1      249250621       <NA>   <NA>
## chr10     135534747       <NA>   <NA>
## chr11     135006516       <NA>   <NA>
## chr12     133851895       <NA>   <NA>
## chr13     115169878       <NA>   <NA>
## ...             ...        ...    ...
## chr8      146364022       <NA>   <NA>
## chr9      141213431       <NA>   <NA>
## chrM          16571       <NA>   <NA>
## chrX      155270560       <NA>   <NA>
## chrY       59373566       <NA>   <NA>
```

We can then use the `tileGenome` function from the **GenomicRanges**
package to tile the genome with 100 kb ranges.


```r
( rng <- tileGenome(seqlengths(BigWigFile(fls[1])), tilewidth=1e5,
                  cut.last.tile.in.chrom=TRUE) )
```

```
## GRanges with 30971 ranges and 0 metadata columns:
##           seqnames               ranges strand
##              <Rle>            <IRanges>  <Rle>
##       [1]     chr1     [     1, 100000]      *
##       [2]     chr1     [100001, 200000]      *
##       [3]     chr1     [200001, 300000]      *
##       [4]     chr1     [300001, 400000]      *
##       [5]     chr1     [400001, 500000]      *
##       ...      ...                  ...    ...
##   [30967]     chrY [58900001, 59000000]      *
##   [30968]     chrY [59000001, 59100000]      *
##   [30969]     chrY [59100001, 59200000]      *
##   [30970]     chrY [59200001, 59300000]      *
##   [30971]     chrY [59300001, 59373566]      *
##   ---
##   seqlengths:
##         chr1     chr10     chr11 ...      chrM      chrX      chrY
##    249250621 135534747 135006516 ...     16571 155270560  59373566
```

Now we are ready to construct a **BigWigFileViews**. This object will
have *columns* which are the files specified by `fileList` and
*rows* which are the ranges specified by `fileRange`.


```r
( bwv <- BigWigFileViews(fileList=fls, fileRange=rng) ) 
```

```
## BigWigFileViews dim: 30971 ranges x 4 samples 
## names: wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep2.bigWig 
## detail: use fileList(), fileSample(), fileRange(), ...
```

Note that the object is lightweight; we simply store the ranges and a
pointer to the files stored on disk.


```r
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
