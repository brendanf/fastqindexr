/*
 * Catch unit tests for fastqindexr (same pattern as optimotu: tests in src TUs,
 * harness in test-runner.cpp, run_testthat_tests registered in RcppExports.cpp).
 */
#include <testthat.h>

context("Catch harness smoke") {
  test_that("Catch is wired through testthat::run_cpp_tests()") {
    expect_true(1 + 1 == 2);
  }
}
