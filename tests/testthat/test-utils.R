test_that("range_scroller works", {
  expect_equal(range_scroller(1), "18S_v1")
  expect_equal(range_scroller(3:6)[2], "18S_v4")
})
