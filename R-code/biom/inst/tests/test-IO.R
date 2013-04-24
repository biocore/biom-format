################################################################################
# Use testthat to test file import and resulting class (and values)
################################################################################
library("biom"); library("testthat")
packageVersion("biom")
# # # # TESTS!
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "biom")
min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "biom")
rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "biom")
rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "biom")
rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "biom")
rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "biom")
# Test read biom
x1 = read_biom(min_dense_file)
x2 = read_biom(min_sparse_file)
x3 = read_biom(rich_dense_file)
x4 = read_biom(rich_sparse_file)
x5 = read_biom(rich_dense_char)
x6 = read_biom(rich_sparse_char)



# # # # TESTS!

test_that("Classes are all biom", {
	expect_is(x1, "biom")
	expect_is(x2, "biom")
	expect_is(x3, "biom")
  expect_is(x4, "biom")
	expect_is(x5, "biom")
	expect_is(x6, "biom")
})

test_that("min/rich files have same biom_data", {
	expect_that(biom_data(x2), is_identical_to(biom_data(x4)))
	expect_that(biom_data(x1), is_identical_to(biom_data(x3)))
})

test_that("biom_datas can be manipulated mathematically", {
	expect_that(2*biom_data(x2), is_identical_to(4*biom_data(x4)/2))
	expect_that(2*biom_data(x1)-biom_data(x1), is_identical_to(biom_data(x3)))
})

test_that("empty stuff is NULL", {
	expect_is(sample_metadata(x1), "NULL")
	expect_is(sample_metadata(x1, 2:4), "NULL")
	expect_is(observation_metadata(x1), "NULL")
})

test_that("Expected classes of non-empty components", {
	expect_is(observation_metadata(x3), "data.frame")
	expect_is(observation_metadata(x3, 2:4), "data.frame")
	expect_is(observation_metadata(x3, 3), "data.frame")
	expect_is(observation_metadata(x1), "NULL")
	expect_that(sample_metadata(x3), is_a("data.frame"))	
	expect_that(biom_data(x3), is_a("Matrix"))
	expect_that(header(x3), is_a("list"))
})

test_that("imported biom files are S4", {
	expect_that(isS4(x1), is_true())	
	expect_that(isS4(x2), is_true())
	expect_that(isS4(x3), is_true())
	expect_that(isS4(x4), is_true())
})

test_that("show method output tests",{
	expect_output(x1, "biom object. type:")
	expect_output(x4, "biom object. type:")
})


