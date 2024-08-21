#' Plot a vertebra map
#'
#' `plotvertmap()` plots a map of the supplied vertebrae, optionally adding colors, marks, and text to identify existing and estimated features of the vertebrae.
#'
#' @param x a `regions_data`, `regions_pco`, or `regions_sim` object; the output of a call to [process_measurements()], [svdPCO()], or [simregions()], respectively.
#' @param type string; the labeling of the x-axis of the plot. Either `"count"` to identify the vertebra index (or absolute position when `centraL` is supplied) or `"percent"` to identify the percent vertebra count (or percent total length when `centraL` is supplied). Abbreviations allowed. Default is `"count"`.
#' @param bps an optional vector containing the region breakpoints. One of `bps`, `modelsupport`, or `bpvar` should be specified to display regions breakpoints. See Details.
#' @param modelsupport an optional `regions_modelsupport` object; the output of a call to [modelsupport()]. One of `bps`, `modelsupport`, or `bpvar` should be specified to display regions breakpoints. See Details.
#' @param criterion string; the criterion to use to select the best model for which breakpoints are to be displayed when `modelsupport` is specified. Ignored otherwise. Allowable options include `"aic"` to use the AICc and `"bic"` to use the BIC. Abbreviations allowed. Default is `"aic"`.
#' @param model `numeric`; from which model among the best as determined by `criterion` should breakpoints be selected when `modelsupport` is supplied. Ignored otherwise. 1 is the best model, 2 the second best, etc. Default is 1.
#' @param bpvar an optional `regions_BPvar` object; the output of a call to [calcBPvar()]. One of `bps`, `modelsupport`, or `bpvar` should be specified to display regions breakpoints. See Details.
#' @param bp.sd an optional vector of the standard deviations of the breakpoints (e.g., as calculated by [calcBPvar()]). When `bpvar` is supplied, the weighted standard deviations are used.
#' @param sd.col when `bp.sd` is specified, the color of the mark on plot indicating the standard deviations. Default is black.
#' @param dropNA `logical`; when some vertebrae are missing, e.g., due to subsampling or starting the analysis at a vertebra beyond the first, whether to remove the missing vertebrae from the plot (`TRUE`) or retain them and label them as missing (i.e., lacking a region) (`FALSE`). Default is `FALSE` to retain them.
#' @param text `logical`; whether to print the vertebra index on each vertebra. Default is `FALSE`.
#' @param name an optional string containing a label used on the left side of the plot.
#' @param centraL an optional numeric vector containing centrum length for each vertebra, which is used to change the size of the plotted vertebrae, or a string containing the name of the variable in the original dataset containing centrum length. Should be of length equal to the number of included vertebrae (i.e., the length of the original dataset). Any vertebrae with centrum length of 0 will be omitted.
#' @param reg.lim a vector of breakpoints indicating other region limits, e.g., anatomic regions.
#' @param lim.col when `reg.lim` is specified, the color of the lines separating the regions. Default is black.
#' @param block.cols when breakpoints are specified (i.e., using `bps`, `modelsupport`, or `bpvar`) and `block.lim` is not specified, a vector of color names or hex codes, one for each region. If not specified, [RColorBrewer::brewer.pal()] with `name = "paired"` will be used to generate colors. When `block.lim` is specified, a named list of vectors of color names or hex codes. See Details.
#' @param block.lim a vector of breakpoints indicating the limits of traditional regions, which will be colored using `block.cols`. See Details.
#'
#' @returns A `ggplot` object that can be manipulated using `ggplot2` syntax.
#'
#' @details
#'
#' `plotvertmap()` uses [ggplot2::geom_rect()] to create the plot. The plots are best viewed with a short height and a long width.
#'
#' ### Specifying breakpoints:
#'
#' There are three ways to specify regions in `plotvertmap()`. First is to supply the vector of breakpoints directly to `bps`. Second is to supply a `regions_modelsupport` object to `modelsupport`. When supplied, the `criterion` and `model` arguments can be used to select which of the sets of breakpoints in the object is to be used. `model` selects which breakpoint model is to be used (1 for first best, 2 for second best, etc.), and `criterion` determines which criterion (AICc or BIC) is used to rank the models. Third is to supply ` regions_BPvar` object to `bpvar`. The weighted average breakpoints will be used after rounding (e.g., a weighted average breakpoint of 3.3 will place vertebrae 1, 2, and 3 in a region, and a weighted average breakpoint of 3.9 will place vertebrae 1, 2, 3, and 4 in a region).
#'
#' ### Using `block.cols`:
#'
#' When `block.lim` is specified, `block.cols` must be specified as a list of vectors of colors, with an entry for each "block". Blocks are predefined regions separate from those specified using the above arguments, e.g., traditional regions. For each region, the most common block is found and assigned to that region. A color of that block as supplied in `block.cols` is used to color that region. So, each block needs as many colors as there are regions assigned to it. For example, if regions 1 and 2 are both assigned to block 1 (i.e., because block 1 is the most common block in those regions), the entry in `block.cols` for that block must have (at least) 2 colors. If an incorrect number of colors per block is supplied, an error will be thrown identifying which blocks are lacking colors. See Examples.
#'
#' @example man/examples/example-plotvertmap.R

#' @export
plotvertmap <- function(x, type = "count",
                        bps = NULL,
                        modelsupport = NULL, criterion = "aic", model = 1,
                        bpvar = NULL, bp.sd = NULL, sd.col = "black",
                        dropNA = FALSE, text = FALSE, name = NULL,
                        centraL = NULL,
                        reg.lim = NULL, lim.col = "black",
                        block.cols = NULL, block.lim = NULL) {

  if (inherits(x, "regions_pco") || inherits(x, "regions_data")) {
    Xvar <- .get_pos(x, subset = FALSE)
    eligible_vertebrae <- .get_eligible_vertebrae(x)
    unique_vert <- .get_eligible_vertebrae(x, subset = FALSE)
  }
  else if (inherits(x, "regions_sim")) {
    Xvar <- x$Xvar
    eligible_vertebrae <- sort(unique(Xvar))
    unique_vert <- eligible_vertebrae
  }
  else {
    chk::err("`x` must inherit from class 'regions_pco', 'regions_data', or 'regions_sim")
  }

  chk::chk_flag(dropNA)
  chk::chk_flag(text)

  chk::chk_string(type)
  type <- tolower(type)
  type <- .match_arg(type, c("count", "percent"))

  if (!is.null(name)) {
    chk::chk_string(name)
  }

  border.col <- "white"

  if (!is.null(bps)) {
    if (chk::vld_atomic(bps) && all(is.na(bps))) {
      bps <- NA_real_
    }
    else {
      chk::chk_numeric(bps)
      chk::chk_range(bps, range(unique_vert))
    }

    BPs <- sort(.drop_na(bps))

    names(BPs) <- paste0("breakpoint", seq_along(BPs))
  }
  else if (!is.null(modelsupport)) {
    chk::chk_is(modelsupport, "regions_modelsupport")

    chk::chk_string(criterion)
    criterion <- tolower(criterion)
    criterion <- .match_arg(criterion, c("aic", "bic"))
    model_support_crit <- modelsupport[[switch(criterion, aic = "Model_support", bic = "Model_support_BIC")]]

    if (is.null(model)) {
      model <- 1L
    }
    else {
      chk::chk_whole_number(model)
      chk::chk_range(model, c(1, nrow(model_support_crit)))
    }

    keep <- which(startsWith(names(model_support_crit), "breakpoint"))
    BPs <- .drop_na(unlist(model_support_crit[model, keep]))

    names(BPs) <- paste0("breakpoint", seq_along(BPs))
  }
  else if (!is.null(bpvar)) {
    chk::chk_is(bpvar, "regions_BPvar")
    BPs <- round(drop(bpvar$WeightedBP["wMean",]))

    names(BPs) <- paste0("breakpoint", seq_along(BPs))

    if (is.null(bp.sd)) {
      bp.sd <- drop(bpvar$WeightedBP["wSD",])
    }
  }
  else if (inherits(x, "regions_sim")) {
    BPs <- x$BPs
  }
  else {
    BPs <- numeric()
  }

  nreg <- length(BPs) + 1

  col.by.block <- !is.null(block.lim)
  if (col.by.block) {
    chk::chk_numeric(block.lim)
    chk::chk_not_any_na(block.lim)

    if (anyDuplicated(block.lim) != 0) {
      chk::err("`block.lim` cannot contain duplicated values")
    }

    if (is.null(block.cols)) {
      chk::err("`block.cols` must be specified when `block.lim` is specified")
    }

    chk::chk_list(block.cols)

    if (length(block.lim) != length(block.cols) - 1) {
      chk::err("`block.cols` must have length equal to one greater than that of `block.lim` (i.e., there should be one more region specified than there are region limits)")
    }

    if (!all(vapply(block.cols, function(x) is.character(x) && all(.is_color(x)), logical(1L)))) {
      chk::err("all entries in `block.cols` must be character vectors containing color names or hex codes")
    }

    if (anyDuplicated(lapply(unlist(block.cols), grDevices::col2rgb)) != 0) {
      chk::wrn("duplicated colors were found in `block.cols`")
    }
  }

  max.pos <- .last(unique_vert)
  vert.all <- seq_len(max.pos)

  vertmap <- data.frame(vname = vert.all)

  if (!is.null(centraL)) {
    vertmap$centraL <- 0

    if (!inherits(x, "regions_sim"))
      data_combined <- .get_data_without_pos(x, subset = FALSE)

    if (is.numeric(centraL)) {
      if (!inherits(x, "regions_sim") && length(centraL) == nrow(data_combined)) {
        for (i in unique_vert) {
          vertmap$centraL[vert.all == i] <- mean(centraL[Xvar == i], na.rm = TRUE)
        }
      }
      else if (length(centraL) == max.pos) {
        vertmap$centraL <- centraL
      }
      else if (length(centraL) == length(unique_vert)) {
        vertmap$centraL[match(vert.all, unique_vert)] <- centraL
      }
      else if (length(centraL) == length(eligible_vertebrae)) {
        vertmap$centraL[match(vert.all, eligible_vertebrae)] <- centraL
      }
      else {
        chk::err(sprintf("`centraL` must have length equal to %s",
                         .word_list(c(
                           sprintf("the total number of vertebrae (in this case, %s)",
                                   max.pos),
                           if (max.pos != length(unique_vert))
                             sprintf("the number of unique vertebrae present in the original dataset (in thise case, %s)",
                                     length(unique_vert)),
                           if (length(unique_vert) != length(eligible_vertebrae))
                             sprintf("the number of available vertebrae (in this case, %s)",
                                     length(eligible_vertebrae)),
                           if (!inherits(x, "regions_sim") && length(unique_vert) != nrow(data_combined))
                             sprintf("the number of observations across all specimens (in this case, %s)",
                                     nrow(data_combined))
                         ), "or")))
      }

      chk::chk_not_any_na(vertmap$centraL, "`centraL`")
      chk::chk_gte(vertmap$centraL, 0, "`centraL`")
    }
    else if (chk::vld_string(centraL)) {
      if (inherits(x, "regions_sim")) {
        chk::err("when `x` is a 'regions_sim' object, `centraL` cannot be specified as a string")
      }

      if (!centraL %in% names(data_combined)) {
        chk::err("the value supplied to `centraL` is not the name of a variable in the dataset")
      }

      # Compute mean centrum length across specimens for each vertebra
      for (i in unique_vert) {
        vertmap$centraL[vertmap$vname == i] <- mean(data_combined[[centraL]][Xvar == i], na.rm = TRUE)
      }

      chk::chk_not_any_na(vertmap$centraL, "the variable named in `centraL`")
      chk::chk_gte(vertmap$centraL, 0, "the variable named in `centraL`")
    }
    else {
      chk::err("`centraL` must be a vector of centrum lengths or a string containing them name of the centrum length variable")
    }
  }

  regions <- paste0("Region", seq_len(nreg))

  # Place 1st region
  vertmap$reg <- regions[1]
  # Add other regions
  for (i in seq_along(BPs)) {
    vertmap$reg[vert.all > BPs[i]] <- regions[i + 1]
  }

  # Assign vertebrae not in data value of "Missing"
  vertmap$reg[!vert.all %in% eligible_vertebrae] <- "Missing"
  if (dropNA) {
    # Drop vertebrae not in data
    vertmap <- vertmap[vertmap$reg != "Missing",, drop = FALSE]
  }

  # Indices of remaining vertebrae
  vertmap$ind <- seq_len(nrow(vertmap))

  if (is.null(centraL)) {
    # Set beginning position of each vert in ind:
    vertmap$ind.beg <- vertmap$ind - .5
    vertmap$ind.end <- vertmap$ind + .5

    if (type == "percent") {
      # Calculate % vertebral count:
      vertmap$pct.beg <- (vertmap$ind - 1) / nrow(vertmap)
      vertmap$pct.end <- vertmap$ind / nrow(vertmap)
    }
  }
  else {
    vertmap$L.end <- vertmap$L.beg <- NA_real_
    vertmap$L.pct.end <- vertmap$L.pct.beg <- NA_real_
    L.pos <- cumsum(vertmap$centraL)
    max.L.pos <- L.pos[length(L.pos)]

    vertmap$L.end <- L.pos
    vertmap$L.beg <- c(0, L.pos[-length(L.pos)])

    if (type == "percent") {
      vertmap$L.pct.beg <- vertmap$L.beg / max.L.pos
      vertmap$L.pct.end <- vertmap$L.end / max.L.pos
    }
  }

  # Add indiv identity:
  vertmap$indiv <- name

  ## Define position of blocks if requested:
  if (col.by.block) {
    missing_vert <- vertmap$reg == "Missing"
    if (is.null(names(block.cols))) {
      names(block.cols) <- paste0("block", seq_along(block.cols))
    }
    block.names <- names(block.cols)
    vertmap$trad.reg <- NA_character_

    # Add traditional regions to table:
    for (i in seq_along(block.names)) {
      if (i == 1L) {
        vertmap$trad.reg[!missing_vert &
                           vertmap$vname <= block.lim[i]] <- block.names[i]
      }
      else if (i == length(block.names)) {
        vertmap$trad.reg[!missing_vert &
                           vertmap$vname > block.lim[i-1]] <- block.names[i]
      }
      else {
        vertmap$trad.reg[!missing_vert &
                           vertmap$vname <= block.lim[i] &
                           vertmap$vname > block.lim[i-1]] <- block.names[i]
      }
    }

    # Match traditional regions boundaries and new regions by finding most
    # common traditional region in each estimated region. E.g., if traditional
    # region A is most common traditional region in estimated region 1, all
    # vertebrae in estimated region 1 get assigned traditional region A
    regs <- setdiff(levels(as.factor(vertmap$reg[])), "Missing")
    vertmap$tradreg.corr <- NA_character_
    for (i in seq_along(regs)) {
      t <- table(vertmap$trad.reg[vertmap$reg == regs[i]])
      vertmap$tradreg.corr[vertmap$reg == regs[i]] <- names(which.max(t))
    }

    # Rename estimated regions based on their appartenance to traditional regions
    reg.corr <- unique(vertmap$tradreg.corr)
    reg.ct <- rep.int(NA_integer_, sum(!missing_vert))
    reg.nm <- rep.int(NA_character_, sum(!missing_vert))
    vertmap$Region <- NA_character_

    for (i in seq_along(reg.corr)) {
      reg_sub <- unique(vertmap$reg[!missing_vert & vertmap$tradreg.corr == reg.corr[i]])
      for (j in seq_along(reg_sub)) {
        reg.ct[vertmap$reg[!missing_vert] == reg_sub[j]] <- j
        reg.nm[vertmap$reg[!missing_vert] == reg_sub[j]] <- reg.corr[i]
      }
    }

    vertmap$Region[!missing_vert] <- paste(reg.nm, reg.ct, sep = "_")
    vertmap$Region[missing_vert] <- "Missing"
  }
  else {
    vertmap$Region <- vertmap$reg
  }

  ## Traditional regions limit position ##
  if (!is.null(reg.lim)) {
    if (is.null(lim.col)) {
      lim.col <- "black"
    }

    # Adjust limit position
    if (is.null(centraL)) {
      if (type == "count") {
        reg.lim <- vapply(reg.lim, function(m) {
          vertmap$ind.end[.last(which(vertmap$vname <= m))]
        }, numeric(1L))
      }
      else if (type == "percent") {
        reg.lim <- vapply(reg.lim, function(m) {
          vertmap$pct.end[.last(which(vertmap$vname <= m))]
        }, numeric(1L))
      }
    }
    else {
      if (type == "count") {
        reg.lim <- vapply(reg.lim, function(m) {
          vertmap$L.end[.last(which(vertmap$vname <= m))]
        }, numeric(1L))
      }
      else if (type == "percent") {
        reg.lim <- vapply(reg.lim, function(m) {
          vertmap$L.pct.end[.last(which(vertmap$vname <= m))]
        }, numeric(1L))
      }
    }
  }

  # Correct BP sd position according to % count and % length:
  if (!is.null(bp.sd)) {
    if (length(BPs) == 0) {
      chk::err("`bp.sd` cannot be supplied when no breakpoints are present")
    }
    chk::chk_numeric(bp.sd)
    chk::chk_not_any_na(bp.sd)
    if (length(bp.sd) != length(BPs)) {
      chk::err(sprintf("`bp.sd` must have length equal to the number of breakpoints (in this case, %s)",
                       length(BPs)))
    }
    chk::chk_gte(bp.sd, 0)

    bp.inds <- vapply(BPs, function(b) {
      vertmap$ind[.last(which(vertmap$vname <= b))]
    }, numeric(1L))

    if (is.null(centraL)) {
      bp.sd <- data.frame(
        bp = vertmap$ind.end[bp.inds],
        beg = vertmap$ind.end[bp.inds] - bp.sd,
        end = vertmap$ind.end[bp.inds] + bp.sd)

      if (type == "percent") {
        bp.sd <- data.frame(
          bp = (bp.sd$bp - .5) / nrow(vertmap),
          beg = (bp.sd$beg - .5) / nrow(vertmap),
          end = (bp.sd$end - .5) / nrow(vertmap))
      }
    }
    else {
      bp.sd <- data.frame(
        bp = vertmap$L.end[bp.inds],
        beg = vapply(seq_along(BPs), function(b) {
          ind <- bp.inds[b]
          p <- vertmap$L.end[ind]
          sd <- bp.sd[b]

          if (sd > 1) {
            sd.floor <- floor(sd)
            p <- p - sum(vertmap$centraL[(ind - sd.floor + 1):(ind)])
            sd <- sd - sd.floor
            ind <- ind - sd.floor
          }

          max(0, p - sd * (vertmap$L.end[ind] - vertmap$L.beg[ind]))
        }, numeric(1L)),
        end = vapply(seq_along(BPs), function(b) {
          ind <- bp.inds[b]
          p <- vertmap$L.end[ind]
          sd <- bp.sd[b]

          if (sd > 1) {
            sd.floor <- floor(sd)
            p <- p + sum(vertmap$centraL[(ind + 1):(ind + sd.floor)])
            sd <- sd - sd.floor
            ind <- ind + sd.floor
          }

          min(max.L.pos, p + sd * (vertmap$L.end[ind + 1] - vertmap$L.beg[ind + 1]))
        }, numeric(1L))
      )

      if (type == "percent") {
        bp.sd <- data.frame(
          bp = bp.sd$bp / max.L.pos,
          beg = bp.sd$beg / max.L.pos,
          end = bp.sd$end / max.L.pos)
      }
    }

    # Assign vertical position and deal with overlaps
    y <- rep.int(0, length(BPs))
    y[1] <- 1

    for (i in seq_along(BPs)[-1]) {
      num_ovl <- vapply(seq_len(max(y)), function(l) {
        sum(bp.sd$beg[y == l] < bp.sd$end[i] & bp.sd$end[y == l] >= bp.sd$beg[i])
      }, numeric(1L))

      y[i] <- {
        if (all(num_ovl > 0)) max(y) + 1
        else min(which(num_ovl == 0))
      }
    }

    # Rescale to keep first position at 0 and alternate above and below
    y <- floor(y / 2) * (-1) ^ y

    sd.jit <- .15

    # Rescale to keep centered at 1.5 and differ by sd.jit
    bp.sd$y <- sd.jit * (y - mean(y)) + 1.5
  }

  ## Add vertebra number; dodge SDs if requested
  if (text) {
    vertmap$text.y <- {
      if (is.null(bp.sd)) 1.5
      else .5 + .5 * min(bp.sd$y)
    }
  }

  ## Set color for each region
  if (col.by.block) {
    colors <- unlist(block.cols)
    names(colors) <- unlist(lapply(names(block.cols), function(i) paste(i, seq_along(block.cols[[i]]), sep = "_")))

    needed.colors <- vapply(unique(reg.nm), function(n) max(reg.ct[reg.nm == n]), integer(1L))
    if (!all(lengths(block.cols) >= needed.colors[names(block.cols)])) {
      deficient.regions <- names(block.cols)[lengths(block.cols) < needed.colors[names(block.cols)]]
      d <- data.frame(lengths(block.cols),
                      needed.colors[names(block.cols)],
                      row.names = names(block.cols))
      names(d) <- c("# supplied", "# needed")
      chk::err("not enough colors were supplied for each region specified in `block.cols`:\n\n",
               paste(utils::capture.output(print(d)), collapse = "\n"),
               "\n\nPlease supply additional colors for ", .word_list(deficient.regions))
    }
  }
  else {
    if (!is.null(block.cols)) {
      chk::chk_character(block.cols)

      if (!all(.is_color(block.cols))) {
        chk::err("all values in `block.cols` must be color names or hex codes")
      }

      if (length(block.cols) < nreg) {
        chk::err(sprintf("`block.cols` must have length greater than or equal to the number of regions (in this case, %s)", nreg))
      }

      if (anyDuplicated(lapply(block.cols, grDevices::col2rgb)) != 0) {
        chk::wrn("duplicated colors were found in `block.cols`")
      }
      pal <- block.cols
    }
    else if (nreg < 3) {
      pal <- RColorBrewer::brewer.pal(3, "Paired")
    }
    else if (nreg %% 2 != 0){
      pal <- RColorBrewer::brewer.pal((nreg + 1), "Paired")
      pal <- pal[-(length(pal)/2)]
    }
    else {
      pal <- RColorBrewer::brewer.pal(nreg, "Paired")
    }

    colors <- setNames(pal[seq_len(nreg)], regions)
  }

  colors <- c(Missing = "grey", colors)		# Add color for missing vertebrae

  if (is.null(centraL)) {
    # Initialize plot
    p <- ggplot(data = vertmap)

    if (type == "count") {
      p <- p +
        geom_rect(aes(xmin = .data$ind.beg,
                      xmax = .data$ind.end,
                      ymin = 1, ymax = 2,
                      fill = .data$Region),
                  color = border.col) +
        labs(y = name, x = "Vertebral count")

      if (text) {
        # Add vertebral number on each colored square
        p <- p + geom_text(aes(x = .data$ind,
                               y = .data$text.y,
                               label = .data$vname))
      }
    }
    else if (type == "percent") {
      p <- p +
        geom_rect(aes(xmin = .data$pct.beg,
                      xmax = .data$pct.end,
                      ymin = 1, ymax = 2,
                      fill = .data$Region),
                  color = border.col) +
        labs(y = name, x = "% Vertebral count")

      if (text) {
        p <- p + geom_text(aes(x = (.data$pct.beg + .data$pct.end)/2,
                               y = .data$text.y,
                               label = .data$vname))
      }
    }
  }
  else {
    # Initialize plot with nonzero centrum lengths
    p <- ggplot(data = vertmap[vertmap$centraL > 0,])

    if (type == "count") {
      p <- p +
        geom_rect(aes(xmin = .data$L.beg,
                      xmax = .data$L.end,
                      ymin = 1, ymax = 2,
                      fill = .data$Region),
                  color = border.col) +
        labs(y = name, x = "Centrum position")

      if (text) {
        p <- p + geom_text(aes(x = (.data$L.beg + .data$L.end)/2,
                               y = .data$text.y,
                               label = .data$vname))
      }
    }
    else if (type == "percent") {
      p <- p +
        geom_rect(aes(xmin = .data$L.pct.beg,
                      xmax = .data$L.pct.end,
                      ymin = 1, ymax = 2,
                      fill = .data$Region),
                  color = border.col) +
        labs(y = name, x = "% Total centrum length")

      if (text) {
        p <- p + geom_text(aes(x = (.data$L.pct.beg + .data$L.pct.end)/2,
                               y = .data$text.y,
                               label = .data$vname))
      }
    }
  }

  if (!is.null(reg.lim)) {
    p <- p + geom_vline(xintercept = reg.lim,
                        color = lim.col, lwd = 1)
  }

  if (!is.null(bp.sd)) {
    # Add horizontal line for BPs standard deviation
    p <- p + geom_pointrange(data = bp.sd,
                             aes(xmin = .data$beg,
                                 xmax = .data$end,
                                 y = .data$y,
                                 x = .data$bp),
                             linewidth = 1, color = sd.col,
                             shape = "diamond", size = .8)
  }

  p <- p +
    scale_fill_manual("Legend", values = colors) +
    scale_x_continuous(expand = c(0, 0),
                       labels = {
                         if (type == "percent") scales::percent
                         else waiver()
                       }) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          axis.title.y = element_text(angle = 0, size = 10,
                                      vjust = 0.5, hjust = 1),
          legend.position = "none",
          plot.margin = unit(c(0, 10, 0, 0), "pt"))

  p
}
