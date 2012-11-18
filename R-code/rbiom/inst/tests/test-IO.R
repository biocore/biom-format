# Use testthat to test file import and resulting class (and values)
min_dense_file   <- system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
min_sparse_file  <- system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
rich_dense_file  <- system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
rich_sparse_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")

# Test read_biom
x1 <- read_biom(min_dense_file)
x2 <- read_biom(min_sparse_file)
x3 <- read_biom(rich_dense_file)
x4 <- read_biom(rich_sparse_file)



# # # # TESTS!

test_that("Classes are all biom", {
	expect_that(x1, is_a("biom"))
	expect_that(x2, is_a("biom"))
	expect_that(x3, is_a("biom"))
	expect_that(x4, is_a("biom"))
})

test_that("min/rich files have same abundances", {
	expect_that(abundance(x2), is_identical_to(abundance(x4)))
	expect_that(abundance(x1), is_identical_to(abundance(x3)))
})

test_that("abundances can be manipulated mathematically", {
	expect_that(2*abundance(x2), is_identical_to(4*abundance(x4)/2))
	expect_that(2*abundance(x1)-abundance(x1), is_identical_to(abundance(x3)))
})

test_that("empty stuff is NULL", {
	expect_that(taxonomy(x1, FALSE), is_a("NULL"))
	expect_that(sampleData(x1, FALSE), is_a("NULL"))
	expect_that(tree(x1, FALSE), is_a("NULL"))
})

test_that("Expected classes of non-empty components", {
	expect_that(taxonomy(x3, FALSE), is_a("data.frame"))
	expect_that(sampleData(x3, FALSE), is_a("data.frame"))	
	expect_that(abundance(x3, FALSE), is_a("Matrix"))
	expect_that(header(x3, FALSE), is_a("list"))
})

test_that("imported biom files are S4", {
	expect_that(isS4(x1), is_true())	
	expect_that(isS4(x2), is_true())
	expect_that(isS4(x3), is_true())
	expect_that(isS4(x4), is_true())
})

test_that("show method output tests",{
	expect_that(x1, prints_text("format: Biological Observation Matrix 1.0.0-dev"))
	expect_that(x4, prints_text("format: Biological Observation Matrix 1.0.0-dev"))
})


