test_that("partition_seq_idx contiguous preserves order and duplicates", {
  ids <- c(5, 2, 5, 1, 9, 9, 4)
  out <- partition_seq_idx(ids, n_parts = 3, strategy = "contiguous")
  expect_length(out, 3)
  expect_equal(out[[1]], c(5L, 2L, 5L))
  expect_equal(out[[2]], c(1L, 9L))
  expect_equal(out[[3]], c(9L, 4L))
})

test_that("partition_seq_idx round_robin preserves per-partition order", {
  ids <- c(10, 2, 7, 2, 8, 5, 1)
  out <- partition_seq_idx(ids, n_parts = 3, strategy = "round_robin")
  expect_length(out, 3)
  expect_equal(out[[1]], c(10L, 2L, 1L))
  expect_equal(out[[2]], c(2L, 8L))
  expect_equal(out[[3]], c(7L, 5L))
})

test_that("partition_seq_idx handles empty input and drop_empty", {
  keep_empty <- partition_seq_idx(integer(), n_parts = 4, drop_empty = FALSE)
  expect_length(keep_empty, 4)
  expect_true(all(vapply(keep_empty, length, integer(1L)) == 0L))

  drop_empty <- partition_seq_idx(integer(), n_parts = 4, drop_empty = TRUE)
  expect_length(drop_empty, 0)
})

test_that("partition_seq_idx validates inputs", {
  expect_error(partition_seq_idx(c(1, 2.5), n_parts = 2), "positive whole")
  expect_error(partition_seq_idx(1:4, n_parts = 0), "`n_parts`")
  expect_error(
    partition_seq_idx(1:4, n_parts = 2, drop_empty = NA),
    "`drop_empty`"
  )
})
