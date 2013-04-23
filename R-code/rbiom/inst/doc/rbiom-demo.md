<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{BIOM support in R main vignette}
-->

<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# rbiom package for BIOM file format support in R

This is an [R Markdown document](http://www.rstudio.com/ide/docs/r_markdown). Markdown is a simple formatting syntax for authoring web pages. Further details on [R markdown here](http://www.rstudio.com/ide/docs/r_markdown).

The BIOM file format (canonically pronounced "biome") is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for [the Earth Microbiome Project](http://www.earthmicrobiome.org/) and is a [Genomics Standards Consortium](http://gensc.org/) candidate project. Please see [the biom-format home page](http://biom-format.org/) for more details.

This demo is designed to provide an overview of the rbiom package to get you started using it quickly. The rbiom package itself is intended to be a utility package that will be depended-upon by other packages in the future. It provides I/O functionality, and functions to make it easier to with data from biom-format files. It does not (and probably should not) provide statistical analysis functions. However, it does provide tools to access data from BIOM format files in ways that are extremely common in R (such as `"data.frame"`, `"matrix"`, and `"Matrix"` classes).

**Package versions** at the time (Mon Apr 22 19:53:01 2013) of this build:

```r
library("rbiom")
packageVersion("rbiom")
```

```
## [1] '0.3.5'
```


---

# Read BIOM format
Here is an example importing BIOM formats of different types into R using the `read_biom` function. The resulting data objects in R are given names beginning with `x`.


```r
min_dense_file = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
min_sparse_file = system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
min_dense_file = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
rich_dense_char = system.file("extdata", "rich_dense_char.biom", package = "rbiom")
rich_sparse_char = system.file("extdata", "rich_sparse_char.biom", package = "rbiom")
x1 = read_biom(min_dense_file)
x2 = read_biom(min_sparse_file)
x3 = read_biom(rich_dense_file)
x4 = read_biom(rich_sparse_file)
x5 = read_biom(rich_dense_char)
x6 = read_biom(rich_sparse_char)
x1
```

```
## biom object. 
## type: OTU table 
## matrix_type: dense 
## 5 rows and 6 columns
```


It would be hard to interpret and wasteful of RAM to stream all the data from a BIOM format file to the standard out when printed with `print` or `show` methods. Instead, a brief summary of the contents BIOM data is shown. 


---

# Access BIOM data
To get access to the data in a familiar form appropriate for many standard R functions, we will need to use accessor functions, also included in the rbiom package.

### Core observation data
The core "observation" data is stored in either sparse or dense matrices in the BIOM format file, and sparse matrix support is carried through on the R side via [the Matrix package](http://cran.r-project.org/web/packages/Matrix/index.html). The variables `x1` and `x2`, assigned above from BIOM files, have similar (but not identical) data. They are small enough to observe directly as tables in standard output in an R session:


```r
biom_table(x1)
```

```
## 5 x 6 Matrix of class "dgeMatrix"
##          Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
## GG_OTU_1       0       0       1       0       0       0
## GG_OTU_2       5       1       0       2       3       1
## GG_OTU_3       0       0       1       4       2       0
## GG_OTU_4       2       1       1       0       0       1
## GG_OTU_5       0       1       1       0       0       0
```

```r
biom_table(x2)
```

```
## 5 x 6 sparse Matrix of class "dgCMatrix"
##          Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
## GG_OTU_1       .       .       1       .       .       .
## GG_OTU_2       5       1       .       2       3       1
## GG_OTU_3       .       .       1       4       .       2
## GG_OTU_4       2       1       1       .       .       1
## GG_OTU_5       .       1       1       .       .       .
```


As you can see above, `x1` and `x2` are represented in R by slightly different matrix classes, assigned automatically based on the data. However, most operations in R will understand this automatically and you should not have to worry about the precise matrix class. However, if the R function you are attempting to use is having a problem with these fancier classes, you can easily coerce the data to the simple, standard `"matrix"` class using the `as` function:


```r
as(biom_table(x2), "matrix")
```

```
##          Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
## GG_OTU_1       0       0       1       0       0       0
## GG_OTU_2       5       1       0       2       3       1
## GG_OTU_3       0       0       1       4       0       2
## GG_OTU_4       2       1       1       0       0       1
## GG_OTU_5       0       1       1       0       0       0
```


### Observation Metadata
Observation metadata is metadata associated with the individual units being counted/recorded in a sample, as opposed to measurements of properties of the samples themselves. For microbiome census data, for instance, observation metadata is often a taxonomic classification and anything else that might be known about a particular OTU/species. For other types of data, it might be metadata known about a particular genome, gene family, pathway, etc. In this case, the observations are counts of OTUs and the metadata is taxonomic classification, if present. The absence of observation metadata is also supported, as we see for the minimal cases of `x1` and `x2`, where we see  instead.


```r
observ_meta(x1)
```

```
## Error: could not find function "observ_meta"
```

```r
observ_meta(x2)
```

```
## Error: could not find function "observ_meta"
```

```r
observ_meta(x3)
```

```
## Error: could not find function "observ_meta"
```

```r
observ_meta(x4)[1:2, 1:3]
```

```
## Error: could not find function "observ_meta"
```

```r
class(observ_meta(x4))
```

```
## Error: could not find function "observ_meta"
```


### Sample Metadata
Similarly, we can access metadata -- if available -- that describe properties of the samples. We access this sample metadata using the `sample_meta` function.


```r
sample_meta(x1)
```

```
## Error: could not find function "sample_meta"
```

```r
sample_meta(x2)
```

```
## Error: could not find function "sample_meta"
```

```r
sample_meta(x3)
```

```
## Error: could not find function "sample_meta"
```

```r
sample_meta(x4)[1:2, 1:3]
```

```
## Error: could not find function "sample_meta"
```

```r
class(sample_meta(x4))
```

```
## Error: could not find function "sample_meta"
```



### Plots
The data really is accessible to other R functions.

```r
plot(biom_table(x4))
```

![plot of chunk plot](figure/plot1.png) 

```r
plot(as(biom_table(x4), "vector"), as(biom_table(x1), "vector"), type = "o")
```

![plot of chunk plot](figure/plot2.png) 

```r
boxplot(as(biom_table(x4), "vector"))
```

![plot of chunk plot](figure/plot3.png) 

```r
heatmap(as(biom_table(x4), "matrix"))
```

![plot of chunk plot](figure/plot4.png) 



---

# Write BIOM format
The biom objects in R can be written to a file/connection using the `write_biom` function. If you modified the biom object, this may still work as well, but no guarantees about this as we are still working on internal checks. The following example writes `x4` to a temporary file, then reads it back using `read_biom` and stores it as variable `y`. The exact comparison of these two objects using the `identical` function shows that they are exactly the same in R.

```r
outfile = tempfile()
write_biom(x4, outfile)
y = read_biom(outfile)
identical(x4, y)
```

```
## [1] TRUE
```


Furthermore, it is possible to invoke standard operating system commands through the R `system` function -- in this case to invoke the `diff` command available on Unix-like systems or the `FC` command on Windows -- in order to compare the original and temporary files directly. Note that this is shown here for convenience, but not automatically run with the rest of the script because of the OS-dependence. During development, though, this same command is tested privately and no differences are reported between the files.


```r
# On Unix OSes
system(paste0("diff ", rich_sparse_file, outfile))
# On windows
system(paste0("FC ", rich_sparse_file, outfile))
```



