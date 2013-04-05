
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# rbiom package for BIOM file format support in R

This is an [R Markdown document](http://www.rstudio.com/ide/docs/r_markdown). Markdown is a simple formatting syntax for authoring web pages. Further details on [R markdown here](http://www.rstudio.com/ide/docs/r_markdown).

The BIOM file format (canonically pronounced "biome") is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for [the Earth Microbiome Project](http://www.earthmicrobiome.org/) and is a [Genomics Standards Consortium](http://gensc.org/) candidate project. Please see [the biom-format home page](http://biom-format.org/) for more details.

This demo is designed to provide an overview of the rbiom package to get you started using it quickly. The rbiom package itself is intended to be a utility package that will be depended-upon by other packages in the future. It provides I/O functionality, and functions to make it easier to with data from biom-format files. It does not (and probably should not) provide statistical analysis functions. However, it does provide tools to access data from BIOM format files in ways that are extremely common in R (such as `"data.frame"`, `"matrix"`, and `"Matrix"` classes).

**Package versions** at the time (Thu Apr  4 18:48:44 2013) of this build:

```r
library("rbiom")
packageVersion("rbiom")
```

```
## [1] '0.2.3'
```


---

# Read BIOM format
Here is an example importing a BIOM format file into R using the `read_biom` function.


```r
min_dense_file = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
min_sparse_file = system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
min_dense_file
```

```
## [1] "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/rbiom/extdata/min_dense_otu_table.biom"
```

```r
x1 = read_biom(min_dense_file)
x2 = read_biom(min_sparse_file)
x3 = read_biom(rich_dense_file)
x4 = read_biom(rich_sparse_file)
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
## NULL
```

```r
observ_meta(x2)
```

```
## NULL
```

```r
observ_meta(x3)
```

```
##                   V1                V2                     V3
## GG_OTU_1 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
## GG_OTU_2 k__Bacteria  p__Cyanobacteria    c__Nostocophycideae
## GG_OTU_3  k__Archaea  p__Euryarchaeota     c__Methanomicrobia
## GG_OTU_4 k__Bacteria     p__Firmicutes          c__Clostridia
## GG_OTU_5 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
##                            V4                    V5                V6
## GG_OTU_1 o__Enterobacteriales f__Enterobacteriaceae    g__Escherichia
## GG_OTU_2        o__Nostocales        f__Nostocaceae g__Dolichospermum
## GG_OTU_3 o__Methanosarcinales f__Methanosarcinaceae g__Methanosarcina
## GG_OTU_4   o__Halanaerobiales   f__Halanaerobiaceae  g__Halanaerobium
## GG_OTU_5 o__Enterobacteriales f__Enterobacteriaceae    g__Escherichia
##                                       V7
## GG_OTU_1                             s__
## GG_OTU_2                             s__
## GG_OTU_3                             s__
## GG_OTU_4 s__Halanaerobiumsaccharolyticum
## GG_OTU_5                             s__
```

```r
observ_meta(x4)[1:2, 1:3]
```

```
##                   V1                V2                     V3
## GG_OTU_1 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
## GG_OTU_2 k__Bacteria  p__Cyanobacteria    c__Nostocophycideae
```

```r
class(observ_meta(x4))
```

```
## [1] "data.frame"
```


### Sample Metadata
Similarly, we can access metadata -- if available -- that describe properties of the samples. We access this sample metadata using the `sample_meta` function.


```r
sample_meta(x1)
```

```
## NULL
```

```r
sample_meta(x2)
```

```
## NULL
```

```r
sample_meta(x3)
```

```
##         BarcodeSequence  LinkerPrimerSequence BODY_SITE Description
## Sample1    CGCTTATCGAGA CATGCTGCCTCCCGTAGGAGT       gut   human gut
## Sample2    CATACCAGTAGC CATGCTGCCTCCCGTAGGAGT       gut   human gut
## Sample3    CTCTCTACCTGT CATGCTGCCTCCCGTAGGAGT       gut   human gut
## Sample4    CTCTCGGCCTGT CATGCTGCCTCCCGTAGGAGT      skin  human skin
## Sample5    CTCTCTACCAAT CATGCTGCCTCCCGTAGGAGT      skin  human skin
## Sample6    CTAACTACCAAT CATGCTGCCTCCCGTAGGAGT      skin  human skin
```

```r
sample_meta(x4)[1:2, 1:3]
```

```
##         BarcodeSequence  LinkerPrimerSequence BODY_SITE
## Sample1    CGCTTATCGAGA CATGCTGCCTCCCGTAGGAGT       gut
## Sample2    CATACCAGTAGC CATGCTGCCTCCCGTAGGAGT       gut
```

```r
class(sample_meta(x4))
```

```
## [1] "data.frame"
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



---

# Write BIOM format
The biom objects in R can be written to a file/connection using the `write_biom` function. If you modified the biom object, this may still work as well, but no guarantees about this as we are still working on internal checks. The following example writes `x4` to a temporary file, then reads it back using `read_biom` and stores it as variable `y`. The exact comparison of these two objects using the `identical` function shows that they are exactly the same.

```r
outfile = tempfile()
write_biom(x4, outfile)
y = read_biom(outfile)
identical(x4, y)
```

```
## [1] TRUE
```



