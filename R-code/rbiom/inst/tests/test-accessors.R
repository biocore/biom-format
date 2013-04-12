################################################################################
# Use testthat to test data accessors
################################################################################
library("rbiom"); library("testthat")
packageVersion("rbiom")
# # # # TESTS!
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "rbiom")
rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "rbiom")
# Test read biom
x1 = read_biom(min_dense_file)
x2 = read_biom(min_sparse_file)
x3 = read_biom(rich_dense_file)
x4 = read_biom(rich_sparse_file)
x5 = read_biom(rich_dense_char)
x6 = read_biom(rich_sparse_char)

# Read tables
T1 = biom_table(x1)
T2 = biom_table(x2)
T3 = biom_table(x3)
T4 = biom_table(x4)
T5 = biom_table(x5)
T6 = biom_table(x6)

# # # # TESTS!

test_that("Test that the results of biom_table are matrix classes", {
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


test_that("Test pre-access biom_table subsetting", {
  # multiple values
  expect_equal(T1[1:3, 3:6], biom_table(x1, 1:3, 3:6), label="multiple rows and cols")
  expect_equal(T2[1:3, 3:6], biom_table(x2, 1:3, 3:6), label="multiple rows and cols")
  expect_equal(T3[1:3, 3:6], biom_table(x3, 1:3, 3:6), label="multiple rows and cols")
  expect_equal(T4[1:3, 3:6], biom_table(x4, 1:3, 3:6), label="multiple rows and cols")
  expect_equal(T5[1:3, 3:6], biom_table(x5, 1:3, 3:6), label="multiple rows and cols")
  expect_equal(T6[1:3, 3:6], biom_table(x6, 1:3, 3:6), label="multiple rows and cols")
  # single value
  label = "single row and column"
  expect_equivalent(T1[3, 4], biom_table(x1, 3, 4), label=label)
  expect_equivalent(T2[3, 4], biom_table(x2, 3, 4), label=label)
  expect_equivalent(T3[3, 4], biom_table(x3, 3, 4), label=label)
  expect_equivalent(T4[3, 4], biom_table(x4, 3, 4), label=label)
  expect_equivalent(T5[3, 4], biom_table(x5, 3, 4), label=label)
  expect_equivalent(T6[3, 4], biom_table(x6, 3, 4), label=label)
  # 1 row, multiple cols
  label = "single rows and multiple cols"
  expect_equal(T1[1, 3:6], biom_table(x1, 1, 3:6), label=label)
  expect_equal(T2[1, 3:6], biom_table(x2, 1, 3:6), label=label)
  expect_equal(T3[1, 3:6], biom_table(x3, 1, 3:6), label=label)
  expect_equal(T4[1, 3:6], biom_table(x4, 1, 3:6), label=label)
  expect_equal(T5[1, 3:6], biom_table(x5, 1, 3:6), label=label)
  expect_equal(T6[1, 3:6], biom_table(x6, 1, 3:6), label=label)
  # 1 column, multiple rows
  label = "single column and multiple rows"
  expect_equal(T1[2:5, 3], biom_table(x1, 2:5, 3), label=label)
  expect_equal(T2[2:5, 3], biom_table(x2, 2:5, 3), label=label)
  expect_equal(T3[2:5, 3], biom_table(x3, 2:5, 3), label=label)
  expect_equal(T4[2:5, 3], biom_table(x4, 2:5, 3), label=label)
  expect_equal(T5[2:5, 3], biom_table(x5, 2:5, 3), label=label)
  expect_equal(T6[2:5, 3], biom_table(x6, 2:5, 3), label=label)
})

