################################################################################
# Use testthat to test data accessors
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

# Test ncol, nrow, colnames, rownames
test_that("Test that ncol, nrow, colnames, rownames, all work as expected", {
  expect_equivalent(ncol(x1), 6L)
  expect_equivalent(ncol(x2), 6L)
  expect_equivalent(ncol(x3), 6L)
  expect_equivalent(ncol(x4), 6L)
  expect_equivalent(ncol(x5), 6L)
  expect_equivalent(ncol(x6), 6L)
  
  expect_equivalent(nrow(x1), 5L)
  expect_equivalent(nrow(x2), 5L)
  expect_equivalent(nrow(x3), 5L)
  expect_equivalent(nrow(x4), 5L)
  expect_equivalent(nrow(x5), 5L)
  expect_equivalent(nrow(x6), 5L)
  
  expect_equivalent(colnames(x1), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"))
  expect_equivalent(colnames(x2), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"))
  expect_equivalent(colnames(x3), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"))
  expect_equivalent(colnames(x4), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"))
  expect_equivalent(colnames(x5), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"))
  expect_equivalent(colnames(x6), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"))
  
  expect_equivalent(rownames(x1), c("GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"))
  expect_equivalent(rownames(x2), c("GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"))
  expect_equivalent(rownames(x3), c("GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"))
  expect_equivalent(rownames(x4), c("GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"))
  expect_equivalent(rownames(x5), c("GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"))
  expect_equivalent(rownames(x6), c("GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"))
  
})

# Read tables
T1 = biom_data(x1)
T2 = biom_data(x2)
T3 = biom_data(x3)
T4 = biom_data(x4)
T5 = biom_data(x5)
T6 = biom_data(x6)

# # # # TESTS!

test_that("Test that the results of biom_data are matrix classes", {
  expect_is(T1, "Matrix")
  expect_is(T2, "Matrix")
  expect_is(T3, "Matrix")
  expect_is(T4, "Matrix")
  expect_is(T5, "matrix")
  expect_is(T6, "matrix")
})

test_that("Some arbitrary test values match expected", {
  expect_equal(T1[5, 1], 0L)
  expect_equal(T1[3, 4], 4L)
  expect_equal(T2[3, 4], 4L)
  expect_equal(T3[3, 4], 4L)
  expect_equal(T4[3, 4], 4L)
  expect_equal(T5[3, 4], "4")
  expect_equal(T6[3, 4], "4")
  expect_equal(T2[5, 1], 0L)
  expect_equal(T3[4, 6], 1L)
  expect_equal(T3[2, 3], 0L)  
  expect_equal(T4[4, 5], 0L)
  expect_equal(T4[1, 3], 1L)  
  expect_equal(T5[1, 5], "clouds")
  expect_equal(T5[5, 1], "0")
  expect_equal(T5[3, 2], "lightning")  
  expect_equal(T6[3, 4], "4")
  expect_equal(T6[2, 5], "bottle")
  expect_equal(T6[4, 5], NA_character_)
})


test_that("Test pre-access biom_data subsetting", {
  label = "multiple rows and multiple columns"
  expect_equal(T1[1:3, 3:6], biom_data(x1, 1:3, 3:6), label=label)
  expect_equal(T2[1:3, 3:6], biom_data(x2, 1:3, 3:6), label=label)
  expect_equal(T3[1:3, 3:6], biom_data(x3, 1:3, 3:6), label=label)
  expect_equal(T4[1:3, 3:6], biom_data(x4, 1:3, 3:6), label=label)
  expect_equal(T5[1:3, 3:6], biom_data(x5, 1:3, 3:6), label=label)
  expect_equal(T6[1:3, 3:6], biom_data(x6, 1:3, 3:6), label=label)

  label = "single row and column (single value)"
  expect_equivalent(T1[3, 4], biom_data(x1, 3, 4), label=label)
  expect_equivalent(T2[3, 4], biom_data(x2, 3, 4), label=label)
  expect_equivalent(T3[3, 4], biom_data(x3, 3, 4), label=label)
  expect_equivalent(T4[3, 4], biom_data(x4, 3, 4), label=label)
  expect_equivalent(T5[3, 4], biom_data(x5, 3, 4), label=label)
  expect_equivalent(T6[3, 4], biom_data(x6, 3, 4), label=label)

  label = "single rows and multiple cols"
  expect_equal(T1[1, 3:6], biom_data(x1, 1, 3:6), label=label)
  expect_equal(T2[1, 3:6], biom_data(x2, 1, 3:6), label=label)
  expect_equal(T3[1, 3:6], biom_data(x3, 1, 3:6), label=label)
  expect_equal(T4[1, 3:6], biom_data(x4, 1, 3:6), label=label)
  expect_equal(T5[1, 3:6], biom_data(x5, 1, 3:6), label=label)
  expect_equal(T6[1, 3:6], biom_data(x6, 1, 3:6), label=label)

  label = "single column and multiple rows"
  expect_equal(T1[2:5, 3], biom_data(x1, 2:5, 3), label=label)
  expect_equal(T2[2:5, 3], biom_data(x2, 2:5, 3), label=label)
  expect_equal(T3[2:5, 3], biom_data(x3, 2:5, 3), label=label)
  expect_equal(T4[2:5, 3], biom_data(x4, 2:5, 3), label=label)
  expect_equal(T5[2:5, 3], biom_data(x5, 2:5, 3), label=label)
  expect_equal(T6[2:5, 3], biom_data(x6, 2:5, 3), label=label)
})


test_that("Test observation metadata extraction", {
	expect_equal(as(observation_metadata(x3, 1)[1:3], "character"),
			c("k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria"))
	expect_equal(as(observation_metadata(x3, 2)[3:5], "character"),
							 c("c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae"))
	expect_equal(as(observation_metadata(x4, 5)[4], "character"), "o__Enterobacteriales")
})


test_that("Test sample metadata extraction", {
	expect_equal(sample_metadata(x3, 1)$Description, "human gut")
	expect_equal(sample_metadata(x3, 2)$BODY_SITE, "gut")
	expect_equal(sample_metadata(x3, 4)$BODY_SITE, "skin")	
	expect_equal(sample_metadata(x4, 4)$BODY_SITE, "skin")	
})


test_that("Test metadata bad-index warnings", {
	expect_warning(out<-sample_metadata(x3, 2:8))
	label = "sample_metadata() output after corrected index has incorrect dimensions"
	expect_equal(dim(out), c(5, 4), label=label)
	expect_warning(out<-observation_metadata(x3, 2:8))
	label = "observation_metadata() output after corrected index has incorrect dimensions"
	expect_equal(dim(out), c(4, 7), label=label)
	expect_warning(out<-observation_metadata(x3, 8:10))
	expect_error(out<-observation_metadata(x3, c("non-id", "also-not-there")))
	expect_error(out<-sample_metadata(x3, c("non-id", "also-not-there")))
})
