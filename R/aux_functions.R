# Clean version of match.arg() for a single choice
.match_arg <- function(arg, choices) {
  arg_name <- deparse1(substitute(arg))

  if (is.null(arg)) return(choices[1L])

  chk::chk_string(arg, arg_name)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)

  if (all(i == 0L)) {
    if (length(choices) == 2L) {
      chk::err("`", arg_name, "` should be one of ", paste(dQuote(choices, FALSE), collapse = " or "))
    }
    else {
      chk::err("`", arg_name, "` should be one of ", paste(dQuote(choices[-length(choices)], FALSE), collapse = ", "),
               ", or ", dQuote(choices[length(choices)], FALSE))
    }
  }

  choices[i]
}

# Similar to na.omit() without appending attributes
.drop_na <- function(x) {
  x[!is.na(x)]
}

# Similar to pbapply::pblapply() for handling cluster argument, but
# with no progress bar.
.lapply_selector <- function(..., cl = NULL) {
  if (is.null(cl)) {
    lapply(...)
  }
  else if (inherits(cl, "cluster")) {
    if (isTRUE(getOption("pboptions")$use_lb))
      parallel::parLapplyLB(cl, ...)
    else parallel::parLapply(cl, ...)
  }
  else {
    parallel::mclapply(..., mc.cores = as.integer(cl))
  }
}

# Similar to rbind() or dplyr::bind_rows(), but specifically to append a df
# to another with fewer columns
.rbind_larger <- function(x, y) {

  if (is.data.frame(x) || is.data.frame(y)) {
    if (is.null(x)) return(y)
    if (is.null(y)) return(x)

    x <- as.data.frame(x)
    y <- as.data.frame(y)

    nam <- {
      if (ncol(x) > ncol(y)) names(x)
      else names(y)
    }

    .expand <- function(d, nam) {
      nam_not_in_d <- setdiff(nam, names(d))
      to_add <- as.data.frame(matrix(nrow = nrow(d), ncol = length(nam_not_in_d),
                                     dimnames = list(rownames(d), nam_not_in_d)))
      cbind(d, to_add)[nam]
    }
  }
  else if (is.matrix(x) && is.matrix(y)) {
    if (is.null(x)) return(y)
    if (is.null(y)) return(x)

    nam <- {
      if (ncol(x) > ncol(y)) names(x)
      else names(y)
    }

    .expand <- function(d, nam) {
      nam_not_in_d <- setdiff(nam, colnames(d))
      to_add <- matrix(nrow = nrow(d), ncol = length(nam_not_in_d),
                       dimnames = list(rownames(d), nam_not_in_d))
      cbind(d, to_add)[, nam, drop = FALSE]
    }
  }
  else {
    stop("`x` and `y` must be data frames or matrices")
  }

  if (ncol(x) > ncol(y)) {
    rbind(x, .expand(y, nam))
  }
  else if (ncol(x) < ncol(y)) {
    rbind(.expand(x, nam), y)
  }
  else {
    rbind(x, y)
  }
}

#When given a vector of strings, creates a string of the form "a and b"
#or "a, b, and c"
.word_list <- function(word.list = NULL, and.or = "and") {

  word.list <- word.list[!word.list %in% c(NA_character_, "")]
  L <- length(word.list)

  if (L == 0) return("")

  if (L == 1) return(word.list)


  and.or <- .match_arg(and.or, c("and", "or"))
  if (L == 2) {
    out <- paste(word.list, collapse = paste0(" ", and.or, " "))
  }
  else {
    out <- paste(paste(word.list[seq_len(L - 1)], collapse = ", "),
                 word.list[L], sep = paste0(", ", and.or, " "))

  }

  out
}

# Checks if supplied argument is a color; vectorized
.is_color <- function(x) {
  vapply(x, function(z) {
    tryCatch(is.matrix(grDevices::col2rgb(z)),
             error = function(e) FALSE)
  }, logical(1L))
}

# Extracts vertebra positions for remaining observations from a
# regions_data dataset
.get_pos <- function(x, subset = TRUE) {
  if (inherits(x, "regions_pco")) {
    x <- attr(x, "data")
  }

  pos_ind <- attr(x, "pos_ind")

  pos <- unlist(lapply(x, `[[`, pos_ind))

  if (!subset) return(pos)

  pos[pos %in% attr(x, "eligible_vertebrae")]
}

# Extracts measurements from
# regions_data dataset
.get_data_without_pos <- function(x, subset = TRUE) {
  if (inherits(x, "regions_pco")) {
    x <- attr(x, "data")
  }

  pos_ind <- attr(x, "pos_ind")

  out <- do.call("rbind", x)[-pos_ind]

  if (!subset) return(out)

  out[.get_pos(x, FALSE) %in% .get_eligible_vertebrae(x),]
}

.get_eligible_vertebrae <- function(x, subset = TRUE) {
  if (inherits(x, "regions_pco")) {
    x <- attr(x, "data")
  }

  if (subset) return(attr(x, "eligible_vertebrae"))

  pos <- .get_pos(x, subset = FALSE)

  sort(unique(pos))
}

# Fast lm(), uses .lm.fit() but accommodates weights
.fast_lm <- function(x, y, w = NULL) {
  if (is.null(w)) return(.lm.fit(x, y))

  w <- sqrt(w)
  out <- .lm.fit(x * w, y * w)
  out$residuals <- out$residuals / w

  out
}

# Get last element in a vector
.last <- function(x) {
  x[length(x)]
}

# Test whether each row of matrix contains any elements in vec
.any_mat_in <- function(m, vec) {
  comp <- {
    if (length(vec) > 1) function(x, y) {x %in% y}
    else function(x, y) {x == y}
  }

  .rowSums(matrix(comp(m, vec), ncol = ncol(m)),
           m = nrow(m), n = ncol(m)) > 0
}

# Test whether each row of matrix contains all elements in vec
.all_mat_in <- function(m, vec) {
  p <- ncol(m)

  out <- .rowSums(matrix(m == vec[1], ncol = p),
                  m = nrow(m), n = p) > 0

  for (i in seq_along(vec)[-1]) {
    out[out][.rowSums(matrix(m[out,] == vec[i], ncol = p),
                      m = sum(out), n = p) == 0] <- FALSE
  }

  out
}

# Fast, minimal version of pmax() that takes in vector and scalar
.pmax2 <- function(x, n = 0) {
  x[x < n] <- n
  x
}

# Version of chk::wrn() that prints warnings immediately instead of waiting till
# the end
.wrn_immediate <- function(...) {
  op <- options()
  on.exit(options(op))
  w <- getOption("warn")
  if (!chk::vld_whole_number(w) || w < 2) {
    options(warn = 1)
  }
  chk::wrn(...)
}

# Checks if two vector of numbers are equal within tolerance
.equiv <- function(x, y, tol = sqrt(.Machine$double.eps)) {
  chk::chk_numeric(x)
  chk::chk_numeric(y)
  chk::chk_number(tol)
  chk::chk_gte(tol, 0)

  abs(x - y) < tol
}