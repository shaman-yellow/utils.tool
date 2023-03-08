# ==========================================================================
# Merge data tables based on continuous variables.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases align_merge
#'
#' @title Alignment merge of 'features' via m/z and RT
#'
#' @description Merge two data.frame of 'features' according to m/z and RT.
#' The formula is the same in MZmine2:
#' Score = (1 - rt.difference / rt.tolerance) * rt.weight +
#' (1 - mz.defference / mz.tolerance) * mz.weight
#'
#' @name align_merge
NULL
#> NULL

#' @export align_merge
#' @aliases align_merge
#' @description \code{align_merge}: ...
#' @rdname align_merge
align_merge <- 
  function(main, sub, id = ".features_id",
    mz.main = "mz", mz.sub = "mz",
    rt.main = "rt.min", rt.sub = "rt.min",
    mz.tol = .05, rt.tol = .3, mz.weight = 75, rt.weight = 25,
    unique = T
    ) {
    lst <- list(main = main, sub = sub)
    check_col <- function(col.x, col.y, lst) {
      if (col.x == col.y) {
        lst <<- coln_suffix(lst, col.x)
        cols <- paste0(col.x, c(".main", ".sub"))
      } else {
        cols <- c(col.x, col.y)
      }
      return(cols)
    }
    mz.cols <- check_col(mz.main, mz.sub, lst)
    rt.cols <- check_col(rt.main, rt.sub, lst)
    data <- tol_merge(lst[[1]], lst[[2]], main_col = mz.cols[1],
      sub_col = mz.cols[2], tol = mz.tol)
    data <- dplyr::mutate(data, .rt_diff = abs(.data[[rt.cols[1]]] - .data[[rt.cols[2]]]),
      .mz_diff = abs(.data[[mz.cols[1]]] - .data[[mz.cols[2]]]))
    data <- dplyr::filter(data, .rt_diff < !!rt.tol)
    data <- dplyr::mutate(data, .score = (1 - .rt_diff / !!rt.tol) * rt.weight +
      (1 - .mz_diff / !!mz.tol) * mz.weight)
    data <- dplyr::arrange(data, .data[[id]], dplyr::desc(.score))
    if (unique)
      data <- dplyr::distinct(data, .data[[id]], .keep_all = T)
    data <- dplyr::select(data, -.rt_diff, -.mz_diff, -.score)
    tibble::as_tibble(data)
  }

#' @export tol_merge
#' @aliases tol_merge
#' @description \code{tol_merge}: ...
#' @rdname align_merge
tol_merge <- 
  function(main,
    sub,
    main_col = "mz",
    sub_col = "mz",
    tol = 0.002,
    bin_size = 1
    ){
    if (main_col == sub_col) {
      new_name <- paste0(sub_col, ".sub")
      colnames(sub)[colnames(sub) == sub_col] <- new_name
      sub_col <- new_name
    }
    main$...seq <- 1:nrow(main)
    backup <- main
    ## to reduce computation, round numeric for limitation
    ## main
    main$...id <- round(main[[ main_col ]], bin_size)
    ## sub
    sub.x <- sub.y <- sub
    sub.x$...id <- round(sub.x[[ sub_col ]], bin_size)
    sub.y$...id <- sub.x$...id + ( 1 * 10^-bin_size )
    sub <- rbind(sub.x, sub.y)
    ## expand merge
    df <- merge(main, sub, by = "...id", all.x = T, allow.cartesian = T)
    df$...diff <- abs(df[[ main_col ]] - df[[ sub_col ]])
    df <- dplyr::filter(df, ...diff <= !!tol)
    ## get the non-merged
    backup <- backup[!backup$...seq %in% df$...seq, ]
    df <- dplyr::bind_rows(df, backup)
    ## remove the assist col
    dplyr::select(df, -...id, -...diff, -...seq)
  }

`:=` <- rlang::`:=`

coln_suffix <-
  function(lst, col,
    suffix = c(".main", ".sub",
      ifelse(length(lst) <= 2, character(0),
        paste0(".other", 1:(length(lst) - 2))))
    ) {
    lapply(1:length(lst),
      function(n) {
        dplyr::rename(lst[[n]], !!paste0(col, suffix[n]) := col)
      })
  }

