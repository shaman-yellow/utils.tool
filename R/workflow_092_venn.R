# ==========================================================================
# workflow of venn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_venn <- setClass("job_venn", 
  contains = c("job"),
  prototype = prototype(
    pg = "venn",
    info = c("..."),
    cite = "",
    method = "",
    tag = "venn",
    analysis = "Venn 交集"
    ))

job_venn <- function(..., lst = NULL)
{
  if (is.null(lst)) {
    object <- list(...)
  } else {
    object <- lst
  }
  if (all(vapply(object, is, logical(1), "feature"))) {
    snaps <- paste0("- ", vapply(object, snap, character(1)))
    snapAdd_onExit("x", "以下取交集：\n{bind(snaps, co = '\n')}")
    object <- lapply(object, function(x) unlist(x@.Data))
  }
  x <- .job_venn(object = object)
  return(x)
}

setMethod("step0", signature = c(x = "job_venn"),
  function(x){
    step_message("Prepare your data with function `job_venn`.")
  })

setMethod("step1", signature = c(x = "job_venn"),
  function(x, ...){
    step_message("Intersection.")
    p.venn <- new_venn(lst = object(x), ...)
    lab(p.venn) <- paste(sig(x), lab(p.venn))
    x@plots[[ 1L ]] <- namel(p.venn)
    x$.append_heading <- FALSE
    if (identical(parent.frame(1), .GlobalEnv)) {
      job_append_heading(
        x, heading = glue::glue(
          "交集: ", bind(names(object(x)), co = " + ")
        )
      )
    }
    x$.feature <- as_feature(p.venn$ins, x)
    return(x)
  })
