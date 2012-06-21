# Use testthat to test file import and resulting class (and values)
min_dense_file   <- "inst/extdata/min_dense_otu_table.biom"
min_sparse_file  <- "inst/extdata/min_sparse_otu_table.biom"
rich_dense_file  <- "inst/extdata/rich_dense_otu_table.biom"
rich_sparse_file <- "inst/extdata/rich_sparse_otu_table.biom"

# Test read_biom
x1 <- read_biom(min_dense_file)
x2 <- read_biom(min_sparse_file)
x3 <- read_biom(rich_dense_file)
x4 <- read_biom(rich_sparse_file)

# Test the classes are all "BIOM"
expect_that(x1, is_a("BIOM"))
expect_that(x2, is_a("BIOM"))
expect_that(x3, is_a("BIOM"))
expect_that(x4, is_a("BIOM"))

# Test that abundance values are what you expect
expect_that(abundance(x2), is_identical_to(abundance(x4)))
expect_that(abundance(x1), is_identical_to(abundance(x3)))

# Test that empty stuff is NULL
expect_that(taxonomy(x1, FALSE), is_a("NULL"))
expect_that(sampleData(x1, FALSE), is_a("NULL"))
expect_that(tree(x1, FALSE), is_a("NULL"))
