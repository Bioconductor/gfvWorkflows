# BigWig example


```r
library(RUnit)
set.seed(123L)
library(GenomicFiles)
fls <- list.files("data", "*bigWig", full = TRUE)
exp <- c("wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig", 
    "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig", "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig", 
    "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep2.bigWig")
if (!identical(exp, basename(fls))) stop("expected files:\n  ", paste(exp, collapse = "\n  "))

seqinfo <- seqinfo(BigWigFile(fls[1]))
rng <- tileGenome(seqlengths(BigWigFile(fls[1])), tilewidth = 1e+05, cut.last.tile.in.chrom = TRUE)
bwv <- BigWigFileViews(fls, fileRange = rng)

register(SerialParam())
x <- summary(bwv[1:8, ])
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

