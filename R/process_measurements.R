#' Process vertebra measurements
#'
#' `process_measurements()` initializes the analysis workflow by processing a dataset of vertebra measurements into an object usable by \pkg{MorphoRegions}. Such processing includes identifying the vertebra indices and the measurements and filling in missing values.
#'
#' @param data a data.frame containing a column of vertebra indices and measurements for each vertebra, or a list thereof for multiple specimens.
#' @param pos the name or index of the variable in `data` containing the vertebra indices. Default is to use the first column.
#' @param measurements the names or indices of the variables in `data` containing the relevant vertebra measurements. If unspecified, will use all variables other than that specified in `pos`.
#' @param fillNA `logical`; whether to fill in missing values using a simple linear imputation. Default is `TRUE`. See Details.
#'
#' @returns A `regions_data` object, which is a list of data.frames (one for each specimen) with attributes containing metadata.
#'
#' @details
#' Any rows with missing values for all measurements will be removed. When missing values in non-removed rows are present and `fillNA` is set to `TRUE`, `process_measurements()` fills them in if the sequence of missing values is no greater than 2 in length. For numeric variables, it uses a linear interpolation, and for categorical variables, it fills in the missing values with the surrounding non-missing values if they are identical and leaves them missing otherwise. Otherwise, missing values are left as they are.
#'
#' When a list of data frames is supplied to `data`, only the variables named in `measurements` that are common across datasets will be stored as measurement variables.
#'
#' @seealso [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#'
#' @example man/examples/example-process_measurements.R

#' @export
process_measurements <- function(data, pos = 1L, measurements, fillNA = TRUE) {

  if (is.matrix(data) || is.data.frame(data)) {
    data <- list(as.data.frame(data))
  }
  else if (is.list(data) && all(vapply(data, function(i) is.data.frame(i) || is.matrix(i), logical(1L)))) {
    for (i in seq_along(data)) {
      data[[i]] <- as.data.frame(data[[i]])
    }
  }
  else {
    chk::err("`data` must be a matrix, dataframe, or list thereof")
  }

  chk::chk_scalar(pos)
  chk::chk_flag(fillNA)

  pos_names <- character(length(data))
  pos_var_list <- vector("list", length(data))

  for (i in seq_along(data)) {
    if (chk::vld_whole_number(pos) && chk::vld_gte(pos, 1) && chk::vld_lte(pos, ncol(data[[i]]))) {
      pos_names[i] <- names(data[[i]])[as.integer(pos)]
    }
    else if (chk::vld_string(pos) && chk::vld_subset(pos, names(data[[i]]))) {
      pos_names[i] <- pos
    }
    else {
      chk::err("`pos` must be a single value indicating the column in `data` containing the vertebra positions")
    }

    pos_var_list[[i]] <- data[[i]][[pos_names[i]]]
  }

  if (length(unique(pos_names)) > 1) {
    chk::err("the variable identified by `pos` must have the same name in all supplied datasets")
  }

  pos <- pos_names[1]

  for (i in seq_along(pos_var_list)) {
    if (!chk::vld_whole_numeric(pos_var_list[[i]])) {
      chk::err("`pos` must refer to a variable of whole numbers identifying vertebra positions")
    }
  }

  measurements_names_list <- vector("list", length(data))

  for (i in seq_along(data)) {
    if (missing(measurements)) {
      measurements_names_list[[i]] <- setdiff(names(data[[i]]), pos)
    }
    else {
      if (length(measurements) == 0) {
        measurements_names_list[[i]] <- character()
      }
      else if (chk::vld_whole_numeric(measurements) && chk::vld_subset(measurements, seq_len(ncol(data[[i]])))) {
        measurements_names_list[[i]] <- names(data[[i]])[as.integer(measurements)]
      }
      else if (chk::vld_character(measurements) && chk::vld_subset(measurements, names(data[[i]]))) {
        measurements_names_list[[i]] <- measurements
      }
      else {
        chk::err("if supplied, `measurements` must indicate the columns in `data` containing the measurement values")
      }

      if (pos %in% measurements_names_list[[i]]) {
        chk::err("`pos` and `measurements` cannot overlap")
      }
    }
  }

  # Get only measurements that are common across datasets
  measurements <- Reduce(union, measurements_names_list)

  for (i in seq_along(data)) {

    #Subset to specified variables
    data[[i]] <- data[[i]][names(data[[i]]) %in% c(pos, measurements)]

    #Reorder columns to be consistent
    if (i > 1)
      data[[i]] <- data[[i]][names(data[[1]])]

    pos_ind <- match(pos, names(data[[i]]))

    #Rows that are all NA
    all_NA_rows <- which(apply(data[[i]], 1, function(x) all(is.na(x[-pos_ind]))))

    if (length(all_NA_rows) > 0)
      data[[i]] <- data[[i]][-all_NA_rows,]

    if (fillNA) {
      #Fill in missing values
      data[[i]] <- .missingval(data[[i]], pos_ind)

      if (anyNA(data[[i]])) {
        chk::wrn(sprintf("missing values remain in %s because there were sequences of missing values greater than %s in length", if (length(data) == 1) "the dataset" else paste("dataset", i), 2))
      }
    }
  }

  attr(data, "pos_ind") <- match(pos, names(data[[1]]))
  attr(data, "eligible_vertebrae") <- sort(unique(unlist(lapply(data, `[[`, pos))))

  class(data) <- "regions_data"

  data
}

# Takes in a dataframe and fills in missing values
# `max_NA_seq_len` controls max number of consecutive missing values that can
# be imputed
.missingval <- function(data, pos_ind, max_NA_seq_len = 2) {

  if (!anyNA(data[-pos_ind])) return(data)

  chk::chk_count(max_NA_seq_len)

  # Add missing vertebrae
  missing_vertebrae <- setdiff(seq(min(data[[pos_ind]]),
                                   max(data[[pos_ind]])),
                               data[[pos_ind]])

  if (length(missing_vertebrae) > 0) {
    missing_rows <- as.data.frame(matrix(NA, nrow = length(missing_vertebrae),
                                         ncol = ncol(data)))
    names(missing_rows) <- names(data)
    missing_rows[[pos_ind]] <- missing_vertebrae

    data <- rbind(data, missing_rows)
  }

  # Order by vertebra
  o <- order(data[[pos_ind]])

  pos <- data[[pos_ind]]

  # Find strings of 1-2 NAs
  variables_with_missing <- which(vapply(data[!(pos %in% missing_vertebrae), -pos_ind],
                                         anyNA, logical(1L)))

  for (i in variables_with_missing) {
    # Extract variable ordered by vertebra
    dat <- data[o, -pos_ind][[i]]

    miss.par <- which(is.na(dat)) #find which ones are missing
    seqs <- split(miss.par, cumsum(c(1, diff(miss.par) != 1))) #split them into sequences

    l.seq <- which(lengths(seqs) <= max_NA_seq_len) #which strings are two or less
    if (length(l.seq) == 0) next #if no short strings skip
    seqs <- seqs[l.seq]

    for (fill in seqs) { #Fill each string
      before <- min(fill) - 1
      after <- max(fill) + 1
      if (before < 1) before <- after #if at the beginning, use the end points
      if (after > length(dat)) after <- before #if at the end, use beginning points

      if (!is.numeric(dat) && dat[before] != dat[after]) next

      if (is.numeric(dat)) {
        val <- seq(dat[before], dat[after], length.out = length(fill) + 2) #calculate missing as mean of adjacent
        dat[fill] <- val[-c(1, length(val))] #fill in the missing
      }
      else {
        val <- rep(dat[before], length(fill))
        dat[fill] <- val #fill in the missing
      }
    }

    data[o, -pos_ind][[i]] <- dat
  }

  # Remove extra fully missing rows
  subset(data, !pos %in% missing_vertebrae)
}
