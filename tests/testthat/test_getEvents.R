library(amplican)
library(testthat)

test_that("mismatch labels stay with their events across insertion columns", {
  actual <- getEvents(
    Biostrings::DNAStringSet(c("ATTAT", "AAGT")),
    Biostrings::DNAStringSet(c("AC-GT", "ACGT")),
    scores = c(11, 22), ID = "fixture")
  expected <- GenomicRanges::GRanges(
    seqnames = "fixture", ranges = IRanges::IRanges(c(3, 2, 2, 3), width = 1),
    strand = "+", originally = c("", "C", "C", "G"),
    replacement = c("T", "A", "T", "A"),
    type = c("insertion", "mismatch", "mismatch", "mismatch"),
    read_id = c("1", "2", "1", "1"), score = c(11, 22, 11, 11))
  names(expected) <- expected$read_id

  expect_identical(actual, expected)
})

test_that("an insertion cannot reorder mismatch labels within one read", {
  actual <- getEvents(
    Biostrings::DNAStringSet("ATTATAT"),
    Biostrings::DNAStringSet("AC-GTAA"), scores = 17, ID = "fixture")
  expected <- GenomicRanges::GRanges(
    seqnames = "fixture", ranges = IRanges::IRanges(c(3, 6, 2, 3), width = 1),
    strand = "+", originally = c("", "A", "C", "G"),
    replacement = c("T", "T", "T", "A"),
    type = c("insertion", "mismatch", "mismatch", "mismatch"),
    read_id = "1", score = 17)
  names(expected) <- expected$read_id

  expect_identical(actual, expected)
})

test_that("event content is independent of read order and chunk boundaries", {
  pattern <- Biostrings::DNAStringSet(c("ATTAT", "AAGT", "ACGT"))
  subject <- Biostrings::DNAStringSet(c("AC-GT", "ACGT", "ACGT"))
  scores <- c(11, 22, 33)
  extract <- function(ids) {
    events <- getEvents(pattern[ids], subject[ids], scores[ids], ID = "fixture")
    if (length(events)) {
      events$read_id <- as.character(ids[as.integer(events$read_id)])
      names(events) <- events$read_id
    }
    events
  }
  canonical <- function(events) {
    events[order(as.integer(events$read_id), GenomicRanges::start(events),
                 GenomicRanges::end(events), events$type)]
  }
  whole <- canonical(extract(1:3))

  expect_identical(canonical(extract(3:1)), whole)
  expect_identical(canonical(c(extract(1L), extract(2L), extract(3L))), whole)
  expect_identical(canonical(c(extract(1:2), extract(3L))), whole)
  expect_identical(extract(3L), GenomicRanges::GRanges())
})

test_that("mismatch extraction preserves per-read offsets and strand metadata", {
  actual <- getEvents(
    Biostrings::DNAStringSet(c("ATTAT", "AAGT")),
    Biostrings::DNAStringSet(c("AC-GT", "ACGT")),
    scores = c(11, 22), ID = "fixture", ampl_shift = c(4L, 8L),
    ampl_start = c(2L, 5L), strand_info = "-")
  expected <- GenomicRanges::GRanges(
    seqnames = "fixture", ranges = IRanges::IRanges(c(7, 13, 6, 7), width = 1),
    strand = "-", originally = c("", "C", "C", "G"),
    replacement = c("T", "A", "T", "A"),
    type = c("insertion", "mismatch", "mismatch", "mismatch"),
    read_id = c("1", "2", "1", "1"), score = c(11, 22, 11, 11))
  names(expected) <- expected$read_id

  expect_identical(actual, expected)
})
