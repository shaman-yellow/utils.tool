# ==========================================================================
# class
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @exportClass graph
#'
#' @aliases graph
#'
#' @title A class built on 'grobs'
#'
#' @description ...
#'
#' @rdname graph-class
#'
#' @examples
#' \dontrun{
#' new('graph', ...)
#' }
graph <- setClass("graph", 
    contains = character(),
    representation = 
      representation(grob = "ANY",
        cvp = "ANY"
        ),
      prototype = NULL
  )

.grob_class <- c("grob", "frame", "gTree", "null",
  "text", "circle", "segments", "gtable",
  "curve", "polygon", "rastergrob")
lapply(.grob_class, setClass, where = topenv())
setClass("viewport")
setClassUnion("grob.obj", .grob_class)

.gg <- c("gg", "ggplot", "ggraph", "patchwork", "ggtree", "ggside")
setFakeClasses(.gg)
setClassUnion("gg.obj", .gg)

.class_unit <- c("unit", "simpleUnit", "unit_v2")
lapply(.class_unit, setClass, where = topenv())
setClassUnion("units", .class_unit)

# ==========================================================================
# color
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @export .default_color
.default_color <- ggsci::pal_npg()(9)

# ==========================================================================
# method
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases draw
#'
#' @title draw class of 'graph' and 'grobs'
#'
#' @description ...
#'
#' @name draw-methods
NULL
#> NULL

setGeneric("draw", 
  function(x, content)
    standardGeneric("draw"))

#' @exportMethod draw
setMethod("draw", signature = c(x = "graph", content = "grob.obj"),
  function(x, content){
    grid.draw(x@grob)
    pushViewport(x@cvp)
    grid.draw(content)
    upViewport(1)
  })

#' @exportMethod draw
setMethod("draw", signature = c(x = "graph"),
  function(x){
    grid.draw(x@grob)
  })

#' @exportMethod draw
setMethod("draw", signature = c(x = "grob.obj"),
  function(x){
    grid.draw(x)
  })

#' @exportMethod into
#' @title place 'grobs' into 'graph'
#' @description ...
#' @rdname into-methods
setGeneric("into", 
  function(x, content) standardGeneric("into"))

#' @exportMethod into
setMethod("into", signature = c(x = "graph", content = "grob.obj"),
  function(x, content){
    if (is.null(content$vp)) {
      content$vp <- x@cvp
    } else {
      content$vp <- vpStack(x@cvp, content$vp)
    }
    gTree(children = gList(x@grob, content))
  })

#' @exportMethod setvp
#' @title ...
#' @description ...
#' @rdname setvp-methods
setGeneric("setvp", 
  function(x, ...) standardGeneric("setvp"))

#' @exportMethod setvp
setMethod("setvp", signature = c(x = "ANY"),
  function(x, ...){
    viewport(grobX(x, 90), grobY(x, 0),
      grobWidth(x), grobHeight(x), ...)
  })

#' @exportMethod weight
setGeneric("weight", 
  function(x, sub) standardGeneric("weight"))
setMethod("weight", signature = c(x = "ANY", sub = "character"),
  function(x, sub){
    if (isS4(x)) {
      weight <-
        sapply(sub, simplify = FALSE, function(sub) {
          get_weight(slot(x, sub))
      })
    }
    as.list(sort(unlist(weight), decreasing = TRUE))
  })

#' @exportMethod as_grob
#' @title convert 'ggplot' object to 'grobs'
#' @description ...
#' @rdname as_grob-methods
setGeneric("as_grob",
  function(x, ...) standardGeneric("as_grob"))

#' @exportMethod as_grob
setMethod("as_grob", signature = c(x = "gg.obj"),
  function(x){
    ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
  })

.as_grob <- function(x) {
  ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
}

#' @exportMethod as_grob
setMethod("as_grob", signature = c(x = "patchwork"),
  function(x){
    patchwork::patchworkGrob(x)
  })

#' @export get_weight
get_weight <- function(x){
  if (isS4(x)) {
    n <- length(slotNames(x))
    if (n == 1) {
      if (is.list(slot(x, slotNames(x)))) {
        n <- n * length(slot(x, slotNames(x)))
      }
    }
    return(n)
  } else if (is.list(x)) {
    length(x)
  } else {
    length(x)
  }
}

#' @aliases frame
#'
#' @title draw in grid frame
#'
#' @description ...
#'
#' @name frame
NULL
#> NULL

#' @export layout_row
#' @aliases layout_row
#' @description \code{layout_row}: ...
#' @rdname frame
layout_row <- function(weight){
  grid.layout(length(weight), 1, heights = weight)
}

#' @export frame_row
#' @aliases frame_row
#' @description \code{frame_row}: ...
#' @rdname frame
frame_row <- function(weight, data, if.ex){
  do.call(frame_place, .fresh_param(list(type = "row")))
}

#' @export layout_col
#' @aliases layout_col
#' @description \code{layout_col}: ...
#' @rdname frame
layout_col <- function(weight){
  grid.layout(1, length(weight), widths = weight)
}

#' @export frame_col
#' @aliases frame_col
#' @description \code{frame_col}: ...
#' @rdname frame
frame_col <- function(weight, data, if.ex){
  do.call(frame_place, .fresh_param(list(type = "col")))
}

#' @exportMethod frame_place
#' @aliases frame_place
#' @description \code{frame_place}: ...
#' @rdname frame
setGeneric("frame_place", 
  function(weight, data, type, if.ex) 
    standardGeneric("frame_place"))

#' @exportMethod frame_place
setMethod("frame_place", signature = setMissing("frame_place",
    weight = "vector",
    data = "list",
    type = "character"),
  function(weight, data, type){
    fun <- paste0("layout_", type)
    layout <- match.fun(fun)(weight)
    frame <- frameGrob(layout = layout)
    data <- sapply(names(data), simplify = FALSE,
      function(name) {
        grob <- data[[ name ]]
        if (is(grob, "graph"))
          grob <- grob@grob
        gTree(children = gList(grob),
          vp = viewport(name = name))
      })
    names(weight) <- repSuffix(names(weight))
    for (i in seq_along(weight)) {
      i.name <- names(weight)[[ i ]]
      o.name <- gsub("\\.rep\\.[0-9]{1,}$", "", i.name)
      i.grob <- data[[ o.name ]]
      i.grob$vp$name <- i.name
      args <- list(frame, i.grob)
      args[[ type ]] <- i
      frame <- do.call(placeGrob, args)
    }
    frame
  })

setMethod("frame_place", signature = setMissing("frame_place",
    weight = "vector",
    data = "list",
    type = "character",
    if.ex = "logical"),
  function(weight, data, type, if.ex){
    main <- weight[!if.ex]
    sub <- weight[if.ex]
    funs <- list(frame_col, frame_row)
    ## which function
    which <- type == c("col", "row")
    ## ex
    ex <- paste0(names(sub), collapse = "__")
    data[[ ex ]] <- funs[!which][[1]](sub, data)
    ex.w <- list(sum(unlist(sub)))
    names(ex.w) <- ex
    weight <- c(main, ex.w)
    ## main
    funs[which][[1]](weight, data)
  })

#' @aliases setnull
#'
#' @title Set markers crossover viewports
#'
#' @description ...
#'
#' @name setnull
NULL
#> NULL

#' @export setnull
#' @aliases setnull
#' @description \code{setnull}: ...
#' @rdname setnull
setnull <- function(target, args, name = "null"){
  gPath <- grid.grep(gPath(target), vpPath = TRUE, grep = TRUE)
  vpPath <- attr(gPath, "vpPath")
  args <- .fresh_param2(list(name = name, vp = vpPath), args)
  do.call(nullGrob, args)
}

#' @export setnullvp
#' @aliases setnullvp
#' @description \code{setnullvp}: ...
#' @rdname setnull
setnullvp <- function(pattern, args, x, name = NULL, fix = TRUE, perl = FALSE){
  if (fix) pattern <- paste0("::", pattern, "$")
  vpPath <- gsub("ROOT::", "", grepPath(pattern, x = x, perl = perl)[1])
  vpPath <- vpPath(vpPath)
  args <- .fresh_param2(list(name = name, vp = vpPath), args)
  do.call(nullGrob, args)
}

#' @export ruler
ruler <- function(p1, p2){
  segmentsGrob(grobX(p1, 0), grobY(p1, 0),
    grobX(p2, 0), grobY(p2, 0))
}

#' @export grepPath
#' @aliases grepPath
#' @description \code{grepPath}: ...
#' @rdname setnull
grepPath <- function (pattern, x = NULL, grobs = TRUE, viewports = TRUE, perl = FALSE) {
    args <- list(x = x, grobs = grobs, viewports = viewports, print = FALSE)
    dl <- do.call(grid.ls, args)
    if (viewports) {
      keep <- dl$type == "vpListing" | dl$type == "grobListing" | 
        dl$type == "gTreeListing"
    } else {
      keep <- dl$type == "grobListing" | dl$type == "gTreeListing"
    }
    vpPaths <- dl$vpPath[keep]
    vpPaths[grepl(pattern, vpPaths, perl = perl)]
  }


#' @export sort_vpPaths
#' @aliases sort_vpPaths
#' @description \code{sort_vpPaths}: ...
#' @rdname setnull
sort_vpPaths <- function(vpPaths){
  nums <- vapply(vpPaths, FUN.VALUE = 0,
    function(ch) {length(grepRaw("::", ch, all = TRUE))})
  lapply(order(nums), function(n) vpPaths[[ n ]])
}

#' @export u
u <- function(n, unit){
  unit <- as.character(substitute(unit))
  unit(n, unit)
}

#' @export vptest
vptest <- function(r = .7, fill = "lightblue"){
  x11(width = 7, height = 7 * r,)
  pushViewport(viewport(, , .5, .5, gp = gpar(fill = fill)))
}

#' @export sym_chem
sym_chem <- function(smi){
  tmpsvg <- file.path(tempdir(), "tempsvg.svg")
  ChemmineOB::convertToImage("SMI", "SVG", source = smi, toFile = tmpsvg)
  svgtxt <- readLines(tmpsvg)
  svgtxt <- gsub("stroke-width=\"2.0", "stroke-width=\"4.0", svgtxt)
  writeLines(svgtxt, tmpsvg)
  rsvg::rsvg_svg(tmpsvg, tmpsvg)
  .rm_background(.cairosvg_to_grob(tmpsvg))
}

#' @aliases ggather
#'
#' @title a mutate of grid::gTree
#'
#' @description ...
#'
#' @name ggather
NULL
#> NULL

#' @export ggather
#' @aliases ggather
#' @description \code{ggather}: ...
#' @rdname ggather
ggather <- function(..., vp = NULL, gp = NULL){
  objs <- list(...)
  objs <- lapply(objs, function(obj) {
    if (is(obj, "graph"))
      obj@grob
    else
      obj
    })
  glist <- do.call(gList, objs)
  gTree(children = glist, gp = gp, vp = vp)
}

#' @export zo
#' @aliases zo
#' @description \code{zo}: ...
#' @rdname ggather
zo <- function(x, w = .9, h = .9) {
  ggather(x, vp = viewport(, , w, h))
}

# ==========================================================================
# arrow
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases arrow
#'
#' @title draw arrow
#'
#' @description ...
#'
#' @name arrow
NULL
#> NULL

#' @export .gpar_dashed_line
.gpar_dashed_line <- grid::gpar(fill = "black", lty = "dashed", lwd = unit(2, "line"))

#' @export .gpar_dotted_line
.gpar_dotted_line <- grid::gpar(fill = "black", lty = "dashed", lwd = unit(2, "line"))

#' @export parrow
#' @aliases parrow
#' @description \code{parrow}: ...
#' @rdname arrow
parrow <- function(n = 5, col, type = "dashed", lwd = u(1, line)){
  y <- seq(0, 1, length.out = n)[c(-1, -n)]
  segs <- segmentsGrob(rep(.1, n - 2), y,
    rep(.9, n - 2), y,
    gp = gpar(lwd = lwd, fill = col, col = col, lty = type)
  )
  arrs <- segmentsGrob(rep(.8, n - 2), y,
    rep(.9, n - 2), y,
    arrow = arrow(angle = 15, length = unit(.7, "line"), type = "closed"),
    gp = gpar(lwd = lwd, fill = col, col = col)
  )
  ggather(segs, arrs)
}

setClassUnion("maybe_p1p2", c("null", "list"))
setClassUnion("maybe_function", c("function", "logical"))
setClassUnion("maybe_list", c("list", "logical"))

#' @exportMethod garrow
#' @aliases garrow
#' @description \code{garrow}: ...
#' @rdname arrow
setGeneric("garrow", 
  function(p1, p2, args_line, args_arrow, fun_line, fun_arrow, city)
    standardGeneric("garrow"))

#' @exportMethod garrow
setMethod("garrow", signature = setMissing("garrow"),
  function(){
    args_line <- list(square = FALSE, ncp = 10, curvature = .3,
      gp = gpar(fill = "black"))
    args_arrow <- list(angle = 15, length = unit(.7, "line"), type = "closed")
    list(args_line = args_line,
      args_arrow = args_arrow,
      fun_line = match.fun(curveGrob),
      fun_arrow = match.fun(arrow),
      city = NULL
    )
  })

#' @exportMethod garrow
setMethod("garrow", signature = c(p1 = "maybe_p1p2", p2 = "maybe_p1p2"),
  function(p1, p2, args_line, args_arrow,
    fun_line, fun_arrow, city){
    args <- as.list(environment())
    default <- garrow()
    args$args_line <- .fresh_param2(default$args_line, args_line)
    args$args_arrow <- .fresh_param2(default$args_arrow, args_arrow)
    args <- .fresh_param2(default, args)
    reCallMethod("garrow", args)
  })

#' @exportMethod garrow
setMethod("garrow", signature = c(p1 = "null", p2 = "null",
    args_line = "list",
    args_arrow = "list",
    fun_line = "maybe_function",
    fun_arrow = "maybe_function",
    city = "maybe_list"),
  function(p1, p2, args_line, args_arrow,
    fun_line, fun_arrow, city){
    args <- as.list(environment())
    args$p1 <- list(x1 = grobX(p1, 0), y1 = grobY(p1, 0))
    args$p2 <- list(x2 = grobX(p2, 0), y2 = grobY(p2, 0))
    reCallMethod("garrow", args)
  })

#' @exportMethod garrow
setMethod("garrow", signature = c(p1 = "list", p2 = "list",
    args_line = "list",
    args_arrow = "list",
    fun_line = "maybe_function",
    fun_arrow = "maybe_function",
    city = "maybe_list"),
  function(p1, p2, args_line, args_arrow,
    fun_line, fun_arrow, city){
    if (!is.null(fun_arrow))
      args_line$arrow <- do.call(fun_arrow, args_arrow)
    if (!is.null(city)) {
      city <- .fresh_param2(garrow_city_args(), city)
      e <- segmentsGrob(p1$x1, p1$y1, p2$x2, p2$y2)
      if (city$axis == "x") {
        pm <- list(x = p1$x1 + city$shift, y = grobY(e, 0))
      } else if (city$axis == "y") {
        pm <- list(x = grobX(e, 90), y = p1$y1 + city$shift)
      } else {
        stop("city$axis != x & city$axis != y")
      }
      args_line <- .fresh_param2(args_line, city$args_line)
      args_line_mid <- args_line[names(args_line) != "arrow"]
      names(pm) <- c("x2", "y2")
      a1 <- do.call(fun_line, c(p1, pm, args_line_mid))
      names(pm) <- c("x1", "y1")
      if (city$rev) {
        args_line$curvature <- -(args_line$curvature)
      }
      a2 <- do.call(fun_line, c(pm, p2, args_line))
      ## as graph
      vp <- viewport(pm$x, pm$y, .1 * grobHeight(e), .1 * grobHeight(e))
      if (!is.null(city$grob_anno)) {
        if (!is.null(city$vp_anno)) {
          vp <- vpStack(vp, city$vp_anno)
        }
        anno <- ggather(city$grob_anno, vp = vp)
        return(ggather(a1, a2, anno))
      }
      return(graph(grob = ggather(a1, a2), cvp = vp))
    }
    args_line <- c(p1, p2, args_line)
    do.call(fun_line, args_line)
  })

#' @export garrow_city
#' @aliases garrow_city
#' @description \code{garrow_city}: ...
#' @rdname arrow
garrow_city <- function(p1, p2, up, left, shift, gp_line){
  shift <- abs(shift)
  curvature <- 1
  if (!up & left) {
    shift <- -shift
  } else if (!up & !left) {
    curvature <- -1
  } else if (up & left) {
    curvature <- -1
    shift <- -shift
  }
  args <- list(curvature = curvature,
    square = TRUE,
    ncp = 1)
  garrow(p1, p2, list(gp = gp_line), 
    city = list(args_line = args, shift = shift))
}

#' @export garrow_snake
#' @aliases garrow_snake
#' @description \code{garrow_snake}: ...
#' @rdname arrow
garrow_snake <- function(p1, p2, color, lwd = u(1, line), cur = 1){
  garrow(p1, p2, list(curvature = cur, square = TRUE, ncp = 1, inflect = TRUE,
      gp = gpar(lwd = lwd, col = color, fill = color)))
}

#' @export garrow_city_args
#' @aliases garrow_city_args
#' @description \code{garrow_city_args}: ...
#' @rdname arrow
garrow_city_args <- 
  function(shift = u(2, line), axis = "x", mid = .5,
    args_line = list(ncp = 1, curvature = 1, square = TRUE),
    grob_anno = NULL, vp_anno = NULL, rev = FALSE
    ){
    as.list(environment())
  }

#' @export sagnage
#' @aliases sagnage
#' @description \code{sagnage}: ...
#' @rdname arrow
sagnage <- function(grob, left = TRUE, l_gpar, borderF = 1.5, front_len = .1,
  vp_shift = u(1.5, line), ...){
  l_gpar <- .fresh_param2f(gpar(linejoin = "round", fill = "grey85",
      col = "transparent"), l_gpar)
  w <- grobWidth(grob) * borderF
  h <- grobHeight(grob) * borderF
  vp <- viewport(, , w * (1 + front_len), h)
  shift <- front_len / (front_len + 1)
  cw <- 1 / (front_len + 1)
  if (left) {
    poly <- polygonGrob(c(0, front_len, 1, 1, front_len), c(.5, 1, 1, 0, 0),
      vp = vp, gp = l_gpar)
    cvp <- viewport(.5 + shift, width = cw)
  } else {
    poly <- polygonGrob(c(0, 0, 1 - front_len, 1, 1 - front_len),
      c(0, 1, 1, .5, 0), vp = vp, gp = l_gpar)
    cvp <- viewport(.5 - shift, width = cw)
  }
  ggrob <- into(graph(grob = poly, cvp = vpStack(vp, cvp)), grob)
  if (left) {
    vp <- viewport(vp_shift, just = c("left", "centre"))
  } else {
    vp <- viewport(-vp_shift, just = c("left", "centre"))
  }
  c(list(grob_anno = ggrob, vp_anno = vp), list(...))
}

#' @export sagnage_shiny
#' @aliases sagnage_shiny
#' @description \code{sagnage_shiny}: ...
#' @rdname arrow
sagnage_shiny <- function(label, left, color){
  sagnage(gtext(label, gpar(col = "white", fontface = "plain")), 
    left, gpar(fill = color))$grob_anno
}

#' @export maparrow
#' @aliases maparrow
#' @description \code{maparrow}: ...
#' @rdname arrow
maparrow <-
  function(obj, data, pattern = list(), round_cur = .3,
    pos = list(r = list(x = 1), l = list(x = 0),
      t = list(y = 1), b = list(y = 0))
    ){
    .check_columns(data, c("from", "to", "group", "fun", "color", "left", "up",
        "shift"), "data")
    if (is.null(data$dup))
      data$dup <- duplicated(data$group)
    if (is.null(data$sag))
      data$sag <- paste0(substr(form(data$fun), 1, 3), ".")
    target <- unique(c(data$from, data$to))
    nulls <- sapply(target, simplify = FALSE,
      function(name) {
        pattern <- 
          if (!is.null(pattern[[ name ]]))
            pattern[[ name ]]
          else
            name
        sapply(names(pos), simplify = FALSE,
          function(p){
            setnullvp(pattern, pos[[p]], obj)
          })
      })
    bafs <- sapply(target, simplify = FALSE,
      function(name) {
        null <- nulls[[ name ]]
        null <- null[names(null) %in% c("l", "r")]
        lapply(null, function(null) {
          function(w = u(1, line), h = u(3, line)) {
            baf(grobX(null, 0), grobY(null, 0), w, h)
          }
          })
      })
    bafs <- unlist(bafs, FALSE)
    keep <- lapply(c("from", "to"),
      function(ch) {
        pos <- ifelse(data[[ "left" ]], "l", "r")
        paste0(data[[ ch ]], ".", pos)
      })
    bafs <- bafs[names(bafs) %in% unique(unlist(keep))]
    arr_city <- apply(data, 1, simplify = FALSE,
      function(v) {
        pos <- ifelse(as.logical(v[[ "left" ]]), "l", "r")
        garrow_city(nulls[[ v[["from"]] ]][[ pos ]],
          nulls[[ v[["to"]] ]][[ pos ]],
          as.logical(v[[ "up" ]]), as.logical(v[[ "left" ]]),
          u(as.numeric(v[["shift"]]), line),
          gpar(fill = v[["color"]], col = v[["color"]],
            lwd = u(1, line))
        )
      })
    arr_round <- apply(data, 1, simplify = FALSE,
      function(v) {
        pos <- ifelse(as.logical(v[[ "left" ]]), "l", "r")
        cur <- if (pos == "l") round_cur else -round_cur
        cur <- if (as.logical(v[[ "up" ]])) -cur else cur
        garrow(nulls[[ v[["from"]] ]][[ pos ]],
          nulls[[ v[["to"]] ]][[ pos ]],
          list(curvature = cur,
            gp = gpar(fill = v[["color"]], col = v[["color"]],
              lwd = u(1, line)))
        )
      })
    sags <- lapply(seq_along(arr_city),
      function(n) {
        if (data$dup[[n]]) return()
        left <- !data$left[[n]]
        up <- data$up[[n]]
        col <- data$color[[n]]
        gtext <- gtext(data$sag[[n]], gpar(col = "white", fontface = "plain"))
        sag <- sagnage(gtext, left, gpar(fill = col), 1.5)$grob_anno
        vp <- arr_city[[n]]@cvp
        x <- if (left) u(.5, npc) + u(2, line) else u(.5, npc) - u(2, line)
        y <- if (up) .5 + 5 else .5 -5
        vp <- vpStack(vp, viewport(x, y))
        ggather(sag, vp = vp)
      })
    sags <- sags[!vapply(sags, is.null, TRUE)]
    namel(nulls, bafs, arr_city, arr_round, sags, data)
  }

#' @export baf
baf <- function(x, y, width = u(1, line), height = u(3, line)) {
  rect <- grectn(bgp_args = gpar(lty = "solid"))@grob
  clip <- clipGrob(, , .7)
  ggather(clip, rect, vp = viewport(x, y, width, height))
}

# ==========================================================================
# text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases text
#'
#' @title a mutate of grid::textGrob
#'
#' @description ...
#'
#' @name text
NULL
#> NULL

#' @export .font
.font <- "Times"

#' @export .title_gp
.title_gp <- grid::gpar(col = "black", cex = 1, fontfamily = .font, fontface = "bold")

#' @export gtext
#' @aliases gtext
#' @description \code{gtext}: ...
#' @rdname text
gtext <- function(label, gp_arg, form = TRUE, ...){
  args <- list(...)
  args <- .fresh_param2(list(x = 0.5, y = 0.5), args)
  args$label <- if (form) form(label) else label
  args$gp <- .fresh_param2f(.title_gp, gp_arg)
  do.call(textGrob, args)
}

#' @export gtextp
#' @aliases gtextp
#' @description \code{gtextp}: ...
#' @rdname text
gtextp <- function(label, gp_arg, form = TRUE, ...){
  if (missing(gp_arg))
    gp_arg <- list()
  gp_arg$font <- c(plain = 1)
  gtext(label, gp_arg, form, ...)
}

#' @export gtext0
#' @aliases gtext0
#' @description \code{gtext0}: ...
#' @rdname text
gtext0 <- function(label, gp_arg, form = TRUE, ...) {
  gtext(label, gp_arg, form, x = .1, hjust = 0, ...)
}

#' @export form
#' @aliases form
#' @description \code{form}: ...
#' @rdname text
form <- function(label){
  Hmisc::capitalize(gsub("_", " ",  label))
}

#' @export gltext
#' @aliases gltext
#' @description \code{gltext}: ...
#' @rdname text
gltext <- function(label, gp_arg = list(), args = list(),
  l_gp = .gpar_dotted_line, flip = FALSE, borderF = 1.2){
  if (flip) args$rot <- 90
  title <- do.call(gtext, c(list(label, gp_arg), args))
  if (flip) {
    height <- grobHeight(title) * borderF
    seg <- list(.5, 0, .5, u(.5, npc) - height / 2)
    seg. <- list(.5, 1, .5, u(.5, npc) + height / 2)
  } else {
    width <- grobWidth(title) * borderF
    seg <- list(0, .5, u(.5, npc) - width / 2, .5)
    seg. <- list(1, .5, u(.5, npc) + width / 2, .5)
  }
  line <- do.call(segmentsGrob, c(seg, list(gp = l_gp)))
  line. <- do.call(segmentsGrob, c(seg., list(gp = l_gp)))
  ggather(title, line, line.)
}

# grid.draw(gtext("test", list(cex = 5)))

# ==========================================================================
# rect
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases rect
#'
#' @title a mutate of grid::rectGrob
#'
#' @description ...
#'
#' @name rect
NULL
#> NULL

#' @export .rect_gp
.rect_gp <- grid::gpar(fill = "transparent")
#' @export .rect.r
.rect.r <- grid::unit(0.3, "lines")
#' @export .vp.sep
.vp.sep <- grid::unit(0.25, "lines")

#' @export .grecti.vp.p
.grecti.vp.p <- list(width = unit(1, "npc") - .vp.sep,
  height = unit(1, "npc") - .vp.sep,
  clip = "on")

#' @export .grecti.vp
.grecti.vp <- do.call(grid::viewport, .grecti.vp.p)
#' @export .grecto.vp
.grecto.vp <- do.call(grid::viewport, .fresh_param2(.grecti.vp.p, list(clip = "inherit")))

#' @export grect
#' @aliases grect
#' @description \code{grect}: ...
#' @rdname rect
grect <- function(name, tfill = "#E18727FF", bfill = "white",
  t_args, tgp_args, t = roundrectGrob,
  b_args, bgp_args, b = roundrectGrob,
  order = c(1, 2), vp = NULL){
  t <- match.fun(t)
  t_args <- .fresh_param2(list(x = 0.5, y = 0.5,
      width = 0.5, height = 0.5,
      r = .rect.r),
    t_args)
  t_args$gp <- .fresh_param2f(gpar(fill = tfill), tgp_args)
  trect <- do.call(t, t_args)
  ## b
  b <- match.fun(b)
  b_args <- .fresh_param2(list(x = 0.5, y = 0.5,
      r = .rect.r),
    b_args)
  b_args$gp <- .fresh_param2f(gpar(fill = bfill), bgp_args)
  brect <- do.call(b, b_args)
  gTree(children = gList(trect, brect)[order], name = name, vp = vp)
}

#' @export grecti
#' @aliases grecti
#' @description \code{grecti}: ...
#' @rdname rect
grecti <- function(label, cex = 1, x = 0.5, y = 1.005,
  borderF = 2, just = c("center", "top"),
  tfill = "#E18727FF", vp = .grecti.vp,
  order = c(2, 1), cvp_clip = "on",
  cvp_fix = TRUE, ...){
  if (is.list(borderF)) {
    borderF <- borderF[[1]]
    xmax <- TRUE
  } else {
    xmax <- FALSE
  }
  if (is(label, "grob")) {
    gtext <- label
    label <- label$label
  } else {
    gtext <- gtext(label, x = x, y = y, gp_arg = list(cex = cex * borderF), just = just)
  }
  t_args <- list(x = grobX(gtext, 90), y = grobY(gtext, 0),
    width = grobWidth(gtext), height = grobHeight(gtext))
  if (xmax) {
    t_args$width <- u(1.2, npc)
  }
  lst <- list(t_args = t_args, tfill = tfill, order = order)
  args <- list(...)
  args <- .fresh_param2(lst, args)
  args$name <- label
  grect <- do.call(grect, args)
  gtext <- gtext(label, list(cex = cex, col = "white"), vp = setvp(gtext))
  grob <- gTree(children = gList(grect, gtext), name = label, vp = vp)
  g <- grect$children[[ order[2] ]]
  cvp_name <- paste0(form(label), "_content")
  if (cvp_fix) {
    cvp <- viewport(x = 0.5, y = 0,
      width = grobWidth(g),
      height = grobHeight(g) - t_args$height,
      just = c("center", "bottom"),
      clip = cvp_clip,
      name = cvp_name)
  } else {
    cvp <- setvp(g, clip = cvp_clip, name = cvp_name)
  }
  graph(grob = grob, cvp = cvp)
}

#' @export grecti2
#' @aliases grecti2
#' @description \code{grecti2}: ...
#' @rdname rect
grecti2 <- function(label, cex = 1, tfill = "#E18727FF",
  borderF = list(2), ...){
  tgp_args <- list(col = tfill)
  bgp_args <- list(col = tfill)
  args <- list(...)
  args <- c(namel(label, cex, borderF, tfill, tgp_args, bgp_args), args)
  do.call(grecti, args)
}

#' @export grecti3
#' @aliases grecti3
#' @description \code{grecti3}: ...
#' @rdname rect
grecti3 <- function(label, cex = 1, tfill = "#E18727FF",
  borderF = list(2), ...){
  tgp_args <- list(col = "black")
  bgp_args <- list(col = "transparent")
  args <- list(...)
  args <- c(namel(label, cex, borderF, tfill, tgp_args, bgp_args), args)
  do.call(grecti, args)
}

#' @export grecto
#' @aliases grecto
#' @description \code{grecto}: ...
#' @rdname rect
grecto <- function(label, cex = 1, x = 0.5, y = 0.995,
  borderF = 2, just = c("center", "bottom"),
  tfill = "#E18727FF", vp = .grecto.vp,
  order = c(1, 2), cvp_clip = "on", cvp_fix = FALSE, ...){
  args <- c(as.list(environment()), list(...))
  do.call(grecti, args)
}

#' @export grectn
#' @aliases grectn
#' @description \code{grectn}: ...
#' @rdname rect
grectn <- function(bfill = "white", b_args, bgp_args, b = roundrectGrob,
  cvp_clip = "inherit") {
  b <- match.fun(b)
  b_args <- .fresh_param2(list(x = 0.5, y = 0.5, r = .rect.r), b_args)
  b_args$gp <- .fresh_param2f(gpar(fill = bfill, lty = "dotted"), bgp_args)
  brect <- do.call(b, b_args)
  graph(grob = brect, cvp = setvp(brect, clip = cvp_clip))
}

#' @export grectn_frame
#' @aliases grectn_frame
#' @description \code{grectn_frame}: ...
#' @rdname rect
grectn_frame <- function(content, title, zo = TRUE){
  if (zo) content <- zo(content)
  content <- frame_row(c(title = .2, content = 1), namel(title, content))
  rect <- grectn(, , list(lty = "solid"))
  into(rect, content)
}

#' @export lst_grecti
#' @aliases lst_grecti
#' @description \code{lst_grecti}: ...
#' @rdname rect
lst_grecti <- function(names, pal, tar = "slot", fun = grecti, ...){
  sapply(names, simplify = FALSE,
    function(name){
      graph <- match.fun(fun)(form(name), , tfill = pal[[ tar ]], ...)
    })
}

#' @export grectN
#' @aliases grectN
#' @description \code{grectN}: ...
#' @rdname rect
grectN <- function(lab.1, lab.2, gp = gpar(fontface = "plain"),
  bfill = "white"){
  frame <- frame_row(list(lab.1 = 1, seg = .1, lab.2 = 1),
    list(lab.1 = gtext(lab.1, gp),
      seg = segmentsGrob(0, .5, 1, .5),
      lab.2 = gtext(lab.2, gp)))
  into(grectn(bfill, , list(lty = "solid")), frame)
}

#' @export grecta
#' @aliases grecta
#' @description \code{grecta}: ...
#' @rdname rect
grecta <- function(label, cex = 4) {
  grob <- gtext(
    label, list(cex = cex), form = FALSE,
    x = 0, y = u(1, npc),
    just = c("left", "top")
  )
  cvp <- viewport(
    grobWidth(grob), 0, u(1, npc) - grobWidth(grob), .95,
    just = c("left", "bottom")
  )
  graph(grob = grob, cvp = cvp)
}

#' @export gshiny
#' @aliases gshiny
#' @description \code{gshiny}: ...
#' @rdname rect
gshiny <- function(xn = 4, yn = 3,
  xps = seq(0, 1, , xn), yps = seq(0, 1, , yn),
  size = c(.15, .02),
  color = .default_color,
  rect = rectGrob(),
  vp = viewport(clip = "off"),
  cvp = viewport(, , .9, .9)
  )
{
  size <- seq(size[1], size[2], length.out = floor(max(c(length(xps), length(yps))) / 2) + 1)
  xsize <- sym_fill(xps, size)
  xpal <- sym_fill(xps, color)
  args <- list(c(xps, xps), c(rep(1, length(xps)), rep(0, length(xps))),
    c(xsize, xsize), c(xpal, xpal))
  fun <- function(n) {
    circleGrob(args[[1]][n], args[[2]][n], args[[3]][n],
      gp = gpar(fill = args[[4]][n], col = "transparent"))
  }
  xcir <- lapply(1:(2 * length(xps)), fun)
  rep <- xps[xps %in% c(0, 1)]
  rep <- yps[yps %in% rep]
  yps <- yps[!yps %in% rep]
  if (length(rep) > 0) {
    size <- size[-1]
    color <- color[-1]
  }
  ysize <- sym_fill(yps, size)
  ypal <- sym_fill(yps, color)
  args <- list(c(rep(0, length(yps)), rep(1, length(yps))), c(yps, yps),
    c(ysize, ysize), c(ypal, ypal))
  ycir <- lapply(1:(2 * length(yps)), fun)
  args <- c(list(rect), xcir, ycir, list(vp = vp))
  graph(grob = do.call(ggather, args), cvp = cvp)
}

#' @export sym_fill
sym_fill <- function(long, short){
  if (length(long) %% 2 == 0) {
    short <- short[1:(length(long) / 2)]
    res <- c(short, rev(short))
  } else {
    half <- (length(long) - 1) / 2
    res <- c(short[1:(half + 1)], rev(short[1:half]))
  }
  if (any(is.na(res)))
    stop( "any(is.na(res)) == TRUE" )
  else
    return(res)
}

# ==========================================================================
# layers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @export glayer
#' @aliases glayer
#' @description \code{glayer}: ...
#' @rdname rect
glayer <- 
  function(n = 5, to = .2, gp = gpar(fill = "white"), fun = rectGrob){
    grob <- fun(x = seq(0, to, length.out = n),
      y = seq(0, to, length.out = n),
      height = 1 - to,
      width = 1 - to,
      just = c("left", "bottom"),
      gp = gp
    )
    cvp <- viewport(to, to, 1 - to, 1 - to, just = c("left", "bottom"), clip = "on")
    graph(grob = grob, cvp = cvp)
  }

# ==========================================================================
# network
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases network
#'
#' @title Quickly draw random network diagrams
#'
#' @description ...
#'
#' @name network
NULL
#> NULL

#' @export fast_layout
#' @aliases fast_layout
#' @description \code{fast_layout}: ...
#' @rdname network
fast_layout <- function(edges, layout = "fr", nodes = NULL, ...){
  if (!is.null(edges)) {
    if (is(edges, "data.frame")) {
      graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
    } else if (is(edges, "igraph")) {
      graph <- edges
    }
  } else {
    edges <- data.frame(from = nodes[[1]][1], to = nodes[[1]][1])
    edges <- tibble::tibble(from = character(0), to = character(0))
    graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  }
  graph <- tidygraph::as_tbl_graph(graph)
  graph <- dplyr::mutate(graph, centrality_degree = tidygraph::centrality_degree(mode = 'all'),
    cent = centrality_degree)
  ggraph::create_layout(graph, layout, ...)
}

#' @export random_graph
#' @aliases random_graph
#' @description \code{random_graph}: ...
#' @rdname network
random_graph <- function(ids, n = length(ids), e = 4, layout = "fr") {
  df <- data.frame(id = ids, size = rnorm(n, .5, .2))
  edges <- data.frame(id1 = sample(ids, e), id2 = sample(ids, e),
    width = rnorm(e, .5, .2))
  fast_layout(edges, layout, df)
}

# ==========================================================================
# others
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @export sep_legend
sep_legend <- function(p, theme) {
  p.l <- MCnebula2:::.get_legend(p + theme)
  theme$legend.position <- "none"
  p.m <- p + theme
  p <- p + theme
  return(namel(p.l, p.m, p))
}

#' @aliases zoom_pdf
#'
#' @title Zoom in locally pdf to png
#'
#' @description ...
#'
#' @name zoom_pdf
NULL
#> NULL

#' @export zoom_pdf
#' @aliases zoom_pdf
#' @description \code{zoom_pdf}: ...
#' @rdname zoom_pdf
zoom_pdf <- function(file, position = c(.5, .5), size = c(.15, .1), page = 1, dpi = 2000,
  as.grob = TRUE)
{
  position[2] <- 1 - position[2]
  png <- pdftools::pdf_render_page(file, page = page, dpi = dpi)
  png <- png::readPNG(png::writePNG(png))
  wxh <- dim(png)[2:1]
  sel <- lapply(1:2,
    function(n) {
      center <- wxh[n] * position[n]
      long <- wxh[n] * size[n]
      start <- center - long / 2
      end <- center + long / 2
      round(start):round(end)
    })
  res <- png[sel[[2]], sel[[1]], ]
  if (as.grob) {
    res <- rasterGrob(res)
  }
  return(res)
}

# ==========================================================================
# shape
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

simulate_peaks <- function(all_range = list(1:30, 31:60, 61:100, 101:140),
  shift = rnorm(10, 2, 1))
{
  lst <- mapply(shift, seq_along(shift), SIMPLIFY = FALSE, FUN = function(shift, id){
    peak <- mapply(all_range, SIMPLIFY = FALSE,
      FUN = function(range){
        peak <- dnorm(range, median(range) + shift, rnorm(1, 5, 1.2)) *
          rnorm(1, 0.7, 0.15)
      })
    feature <- mapply(seq_along(all_range), lengths(all_range),
      FUN = function(seq, rep){
        rep(paste0("peak", seq), rep)
      })
    tibble::tibble(x = unlist(all_range), y = unlist(peak),
      sample = paste0("sample", id),
      peak = unlist(feature)
    )
  })
  data.table::rbindlist(lst)
}

# ==========================================================================
# get_ggsets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @export nebulae_as_grob
#' @aliases nebulae_as_grob
#' @title Convert Nebulae as 'grob' object
#' @description \code{nebulae_as_grob}: This will convert Child-Nebulae
#' as a 'grob' object. See package `grid` about 'grob' object.
#' @param x [mcnebula-class] object.
#' @rdname nebulae_as_grob
nebulae_as_grob <- function(x) {
  chAsGrob <- function(ch, x) {
    ggset <- modify_default_child(ch)
    as_grob(call_command(ggset))
  }
  sets <- lapply(ggset(child_nebulae(x)), chAsGrob, x = x)
  sets <- lapply(names(sets),
    function(name) {
      ggather(sets[[name]],
        vp = viewports(child_nebulae(x))[[name]])
    })
  sets_vp <- viewport(layout = grid_layout(child_nebulae(x)))
  sets <- do.call(ggather, c(sets, list(vp = sets_vp)))
  legendH <- MCnebula2:::.legend_hierarchy(child_nebulae(x), x)
  legendG <- MCnebula2:::.get_legend(
    call_command(modify_default_child(ggset(child_nebulae(x))[[1]], x))
  )
  ## integrate
  vis <- frame_row(list(sets = 5, legendH = .5), namel(sets, legendH))
  vis <- frame_col(list(vis = 4, legendG = 1), namel(vis, legendG))
  vis
}

setMethod("as_grob", signature = c(x = "expression"),
  function(x, envir = parent.frame(2)){
    ggplotify::base2grob(x, envir)
  })
