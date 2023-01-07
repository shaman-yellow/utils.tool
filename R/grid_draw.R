# ==========================================================================
# class
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @import grid
graph <- 
  setClass("graph", 
           contains = character(),
           representation = 
             representation(grob = "ANY",
                            cvp = "ANY"
                            ),
           prototype = NULL
           )

.grob_class <- c("grob", "frame", "gTree", "null",
                 "text", "circle", "segments", "gtable",
                 "curve", "polygon")
setOldClass(.grob_class)
setOldClass("viewport")
setClassUnion("grob.obj", .grob_class)

.gg <- c("gg", "ggplot", "ggraph")
setOldClass(.gg)
setClassUnion("gg.obj", .gg)

.class_unit <- c("unit", "simpleUnit", "unit_v2")
setOldClass(.class_unit)
setClassUnion("units", .class_unit)

# ==========================================================================
# color
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.default_color <- ggsci::pal_npg()(9)

# ==========================================================================
# method
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setGeneric("draw", 
           function(x, content)
             standardGeneric("draw"))

setMethod("draw", 
          signature = c(x = "graph", content = "grob.obj"),
          function(x, content){
            grid.draw(x@grob)
            pushViewport(x@cvp)
            grid.draw(content)
            upViewport(1)
          })

setMethod("draw", 
          signature = setMissing("draw",
                                 x = "graph"),
          function(x){
            grid.draw(x@grob)
          })

setMethod("draw", 
          signature = c(x = "grob.obj"),
          function(x){
            grid.draw(x)
          })

setGeneric("into", 
           function(x, content) standardGeneric("into"))

setMethod("into", 
          signature = c(x = "graph", content = "grob.obj"),
          function(x, content){
            if (is.null(content$vp)) {
              content$vp <- x@cvp
            } else {
              content$vp <- vpStack(x@cvp, content$vp)
            }
            gTree(children = gList(x@grob, content))
          })

setGeneric("setvp", 
           function(x, ...) standardGeneric("setvp"))

setMethod("setvp", 
          signature = c(x = "ANY"),
          function(x, ...){
            viewport(grobX(x, 90), grobY(x, 0),
                     grobWidth(x), grobHeight(x), ...)
          })

setGeneric("weight", 
           function(x, sub) standardGeneric("weight"))
setMethod("weight", 
          signature = c(x = "ANY", sub = "character"),
          function(x, sub){
            if (isS4(x)) {
              weight <-
                sapply(sub, simplify = F, function(sub) {
                         get_weight(slot(x, sub))
                     })
            }
            as.list(sort(unlist(weight), decreasing = T))
          })

setGeneric("as_grob",
           function(x) standardGeneric("as_grob"))
setMethod("as_grob", 
          signature = c(x = "gg.obj"),
          function(x){
            ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
          })

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

layout_row <- function(weight){
  grid.layout(length(weight), 1, heights = weight)
}

frame_row <- function(weight, data, if.ex){
  do.call(frame_place, .fresh_param(list(type = "row")))
}

layout_col <- function(weight){
  grid.layout(1, length(weight), widths = weight)
}

frame_col <- function(weight, data, if.ex){
    do.call(frame_place, .fresh_param(list(type = "col")))
}

setGeneric("frame_place", 
           function(weight, data, type, if.ex) 
             standardGeneric("frame_place"))
setMethod("frame_place", 
          signature = setMissing("frame_place",
                                 weight = "vector",
                                 data = "list",
                                 type = "character"),
          function(weight, data, type){
            fun <- paste0("layout_", type)
            layout <- match.fun(fun)(weight)
            frame <- frameGrob(layout = layout)
            data <- sapply(names(data), simplify = F,
                           function(name) {
                             grob <- data[[ name ]]
                             if (is(grob, "graph"))
                               grob <- grob@grob
                             gTree(children = gList(grob),
                                   vp = viewport(name = name))
                           })
            names(weight) <- repSuffix(names(weight))
            for (i in 1:length(weight)) {
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
setMethod("frame_place", 
          signature = setMissing("frame_place",
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

setnull <- function(target, args, name = "null"){
  gPath <- grid.grep(gPath(target), vpPath = T, grep = T)
  vpPath <- attr(gPath, "vpPath")
  args <- .fresh_param2(list(name = name, vp = vpPath), args)
  do.call(nullGrob, args)
}

setnullvp <- function(pattern, args, x, name = NULL, fix = T, perl = F){
  if (fix) pattern <- paste0("::", pattern, "$")
  vpPath <- gsub("ROOT::", "", grepPath(pattern, x = x, perl = perl)[1])
  vpPath <- vpPath(vpPath)
  args <- .fresh_param2(list(name = name, vp = vpPath), args)
  do.call(nullGrob, args)
}

ruler <- function(p1, p2){
  segmentsGrob(grobX(p1, 0), grobY(p1, 0),
               grobX(p2, 0), grobY(p2, 0))
}

## grid.grep
grepPath <- 
  function (pattern, x = NULL, grobs = T, viewports = T, perl = F) {
    args <- list(x = x, grobs = grobs, viewports = viewports, print = F)
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
    

## sort vpPaths
sort_vpPaths <- function(vpPaths){
  nums <- vapply(vpPaths, FUN.VALUE = 0,
                 function(ch) {length(grepRaw("::", ch, all = T))})
  lapply(order(nums), function(n) vpPaths[[ n ]])
}

u <- function(n, unit){
  unit <- as.character(substitute(unit))
  unit(n, unit)
}

vptest <- function(r = .7, fill = "lightblue"){
  x11(width = 7, height = 7 * r,)
  pushViewport(viewport(, , .5, .5, gp = gpar(fill = fill)))
}

sym_chem <- function(smi){
  tmpsvg <- paste0(tempdir(), "/tempsvg.svg")
  ChemmineOB::convertToImage("SMI", "SVG", source = smi, toFile = tmpsvg)
  svgtxt <- readLines(tmpsvg)
  svgtxt <- gsub("stroke-width=\"2.0", "stroke-width=\"4.0", svgtxt)
  writeLines(svgtxt, tmpsvg)
  rsvg::rsvg_svg(tmpsvg, tmpsvg)
  .rm_background(.cairosvg_to_grob(tmpsvg))
}

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

zo <- function(x, w = .9, h = .9) {
  ggather(x, vp = viewport(, , w, h))
}

# ==========================================================================
# arrow
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.gpar_dashed_line <- gpar(fill = "black", lty = "dashed", lwd = unit(2, "line"))
.gpar_dotted_line <- gpar(fill = "black", lty = "dashed", lwd = unit(2, "line"))

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
setClassUnion("maybe_function", c("NULL", "function"))
setClassUnion("maybe_list", c("NULL", "list"))

setGeneric("garrow", 
           function(p1, p2, args_line, args_arrow, fun_line, fun_arrow, city)
             standardGeneric("garrow"))
setMethod("garrow", 
          signature = setMissing("garrow"),
          function(){
            args_line <- list(square = F, ncp = 10, curvature = .3,
                              gp = gpar(fill = "black"))
            args_arrow <- list(angle = 15, length = unit(.7, "line"), type = "closed")
            list(args_line = args_line,
                 args_arrow = args_arrow,
                 fun_line = match.fun(curveGrob),
                 fun_arrow = match.fun(arrow),
                 city = NULL
            )
          })
setMethod("garrow", 
          signature = c(p1 = "maybe_p1p2", p2 = "maybe_p1p2"),
          function(p1, p2, args_line, args_arrow,
                   fun_line, fun_arrow, city){
            args <- as.list(environment())
            default <- garrow()
            args$args_line <- .fresh_param2(default$args_line, args_line)
            args$args_arrow <- .fresh_param2(default$args_arrow, args_arrow)
            args <- .fresh_param2(default, args)
            reCallMethod("garrow", args)
          })
setMethod("garrow", 
          signature = c(p1 = "null", p2 = "null",
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
setMethod("garrow", 
          signature = c(p1 = "list", p2 = "list",
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
               square = T,
               ncp = 1)
  garrow(p1, p2, list(gp = gp_line), 
         city = list(args_line = args, shift = shift))
}

garrow_snake <- function(p1, p2, color, lwd = u(1, line), cur = 1){
  garrow(p1, p2, list(curvature = cur, square = T, ncp = 1, inflect = T,
                      gp = gpar(lwd = lwd, col = color, fill = color)))
}

garrow_city_args <- 
  function(shift = u(2, line), axis = "x", mid = .5,
           args_line = list(ncp = 1, curvature = 1, square = T),
           grob_anno = NULL, vp_anno = NULL, rev = F
           ){
    as.list(environment())
  }

sagnage <- function(grob, left = T, l_gpar, borderF = 1.5, front_len = .1,
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

sagnage_shiny <- function(label, left, color){
  sagnage(gtext(label, gpar(col = "white", fontface = "plain")), 
          left, gpar(fill = color))$grob_anno
}

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
    nulls <- sapply(target, simplify = F,
                    function(name) {
                      pattern <- 
                        if (!is.null(pattern[[ name ]]))
                          pattern[[ name ]]
                        else
                          name
                      sapply(names(pos), simplify = F,
                             function(p){
                               setnullvp(pattern, pos[[p]], obj)
                             })
                    })
    bafs <- sapply(target, simplify = F,
                   function(name) {
                     null <- nulls[[ name ]]
                     null <- null[names(null) %in% c("l", "r")]
                     lapply(null, function(null) {
                              function(w = u(1, line), h = u(3, line)) {
                                baf(grobX(null, 0), grobY(null, 0), w, h)
                              }
                             })
                   })
    bafs <- unlist(bafs, F)
    keep <- lapply(c("from", "to"),
                   function(ch) {
                     pos <- ifelse(data[[ "left" ]], "l", "r")
                     paste0(data[[ ch ]], ".", pos)
                   })
    bafs <- bafs[names(bafs) %in% unique(unlist(keep))]
    arr_city <- apply(data, 1, simplify = F,
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
    arr_round <- apply(data, 1, simplify = F,
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
    sags <- lapply(1:length(arr_city),
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
    sags <- sags[!vapply(sags, is.null, T)]
    namel(nulls, bafs, arr_city, arr_round, sags, data)
  }

baf <- function(x, y, width = u(1, line), height = u(3, line)) {
  rect <- grectn(bgp_args = gpar(lty = "solid"))@grob
  clip <- clipGrob(, , .7)
  ggather(clip, rect, vp = viewport(x, y, width, height))
}

# ==========================================================================
# text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.font <- "Times"
.title_gp <- gpar(col = "black", cex = 1, fontfamily = .font, fontface = "bold")

gtext <- function(label, gp_arg, form = T, ...){
  args <- list(...)
  args <- .fresh_param2(list(x = 0.5, y = 0.5), args)
  args$label <- if (form) form(label) else label
  args$gp <- .fresh_param2f(.title_gp, gp_arg)
  do.call(textGrob, args)
}

gtextp <- function(label, gp_arg, form = T, ...){
  if (missing(gp_arg))
    gp_arg <- list()
  gp_arg$font <- c(plain = 1)
  gtext(label, gp_arg, form, ...)
}

gtext0 <- function(label, gp_arg, form = T, ...) {
  gtext(label, gp_arg, form, x = .1, hjust = 0, ...)
}

form <- function(label){
  Hmisc::capitalize(gsub("_", " ",  label))
}

gltext <- function(label, gp_arg = list(), args = list(),
                   l_gp = .gpar_dotted_line, flip = F, borderF = 1.2){
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

.rect_gp <- gpar(fill = "transparent")
.rect.r <- unit(0.3, "lines")
.vp.sep <- unit(0.25, "lines")

.grecti.vp.p <- list(width = unit(1, "npc") - .vp.sep,
                     height = unit(1, "npc") - .vp.sep,
                     clip = "on")
.grecti.vp <- do.call(viewport, .grecti.vp.p)
.grecto.vp <- do.call(viewport, .fresh_param2(.grecti.vp.p, list(clip = "inherit")))

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

# grid.draw(grect("test"))
grecti <- function(label, cex = 1, x = 0.5, y = 1.005,
                   borderF = 2, just = c("center", "top"),
                   tfill = "#E18727FF", vp = .grecti.vp,
                   order = c(2, 1), cvp_clip = "on",
                   cvp_fix = T, ...){
  if (is.list(borderF)) {
    borderF <- borderF[[1]]
    xmax <- T
  } else {
    xmax <- F
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

grecti2 <- function(label, cex = 1, tfill = "#E18727FF",
                    borderF = list(2), ...){
  tgp_args <- list(col = tfill)
  bgp_args <- list(col = tfill)
  args <- list(...)
  args <- c(namel(label, cex, borderF, tfill, tgp_args, bgp_args), args)
  do.call(grecti, args)
}

grecti3 <- function(label, cex = 1, tfill = "#E18727FF",
                    borderF = list(2), ...){
  tgp_args <- list(col = "black")
  bgp_args <- list(col = "transparent")
  args <- list(...)
  args <- c(namel(label, cex, borderF, tfill, tgp_args, bgp_args), args)
  do.call(grecti, args)
}

grecto <- function(label, cex = 1, x = 0.5, y = 0.995,
                   borderF = 2, just = c("center", "bottom"),
                   tfill = "#E18727FF", vp = .grecto.vp,
                   order = c(1, 2), cvp_clip = "on", cvp_fix = F, ...){
  args <- c(as.list(environment()), list(...))
  do.call(grecti, args)
}

grectn <- function(bfill = "white", b_args, bgp_args, b = roundrectGrob,
                   cvp_clip = "inherit") {
  b <- match.fun(b)
  b_args <- .fresh_param2(list(x = 0.5, y = 0.5, r = .rect.r), b_args)
  b_args$gp <- .fresh_param2f(gpar(fill = bfill, lty = "dotted"), bgp_args)
  brect <- do.call(b, b_args)
  graph(grob = brect, cvp = setvp(brect, clip = cvp_clip))
}

grectn_frame <- function(content, title, zo = T){
  if (zo) content <- zo(content)
  content <- frame_row(c(title = .2, content = 1), namel(title, content))
  rect <- grectn(, , list(lty = "solid"))
  into(rect, content)
}

lst_grecti <- function(names, pal, tar = "slot", fun = grecti, ...){
  sapply(names, simplify = F,
         function(name){
           graph <- match.fun(fun)(form(name), , tfill = pal[[ tar ]], ...)
         })
}

grectN <- function(lab.1, lab.2, gp = gpar(fontface = "plain"),
                   bfill = "white"){
  frame <- frame_row(list(lab.1 = 1, seg = .1, lab.2 = 1),
                     list(lab.1 = gtext(lab.1, gp),
                          seg = segmentsGrob(0, .5, 1, .5),
                          lab.2 = gtext(lab.2, gp)))
  into(grectn(bfill, , list(lty = "solid")), frame)
}

grecta <- function(label, cex = 4) {
  grob <- gtext(
    label, list(cex = cex), form = F,
    x = 0, y = u(1, npc),
    just = c("left", "top")
  )
  cvp <- viewport(
    grobWidth(grob), 0, u(1, npc) - grobWidth(grob), .95,
    just = c("left", "bottom")
  )
  graph(grob = grob, cvp = cvp)
}

gshiny <- function(xn = 4, yn = 3,
                   xps = seq(0, 1, , xn), yps = seq(0, 1, , yn),
                   size = c(.15, .02),
                   color = .default_color,
                   rect = rectGrob(),
                   vp = viewport(clip = "off"),
                   cvp = viewport(, , .9, .9)
                   ){
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

sym_fill <- function(long, short){
  if (length(long) %% 2 == 0) {
    short <- short[1:(length(long) / 2)]
    res <- c(short, rev(short))
  } else {
    half <- (length(long) - 1) / 2
    res <- c(short[1:(half + 1)], rev(short[1:half]))
  }
  if (any(is.na(res)))
    stop( "any(is.na(res)) == T" )
  else
    return(res)
}

# ==========================================================================
# get external grob
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.expathsvg <-
  system.file("extdata", "svg",
              package = gsub("^.*:", "",
                             environmentName(topenv(environment()))))
prefix <- c()

.check_external_svg <- function(){
  files <- list.files(.expathsvg, "\\.svg$", full.names = T)
  log.path <- paste0(.expathsvg, "/log")
  if (file.exists(log.path)) {
    log <- readLines(log.path)
    log <- log[vapply(log, file.exists, T)]
  } else {
    log <- c()
  }
  if (length(files) > 0) {
    new.log <-
      lapply(files,
             function(file) {
               if (!file %in% log) {
                 rsvg::rsvg_svg(file, file)
                 file
               }
             })
    new.log <- unlist(new.log)
    log <- c(log, new.log)
  }
  if (!is.null(log))
    writeLines(log, log.path)
}
.check_external_svg()

ex_grob <- function(name, fun = .cairosvg_to_grob){
  file <- paste0(.expathsvg, "/", name, ".svg")
  if (file.exists(file)) {
    fun(file)
  } else {
    stop("file.exsits(file) == F")
  }
}

ex_pic <- function(name){
  ex_grob(name, fun = grImport2::readPicture)
}

# ==========================================================================
# layers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

fast_layout <- function(edges, layout = "fr", nodes = NULL){
  graph <- igraph::graph_from_data_frame(edges, directed = T, vertices = nodes)
  graph <- tidygraph::as_tbl_graph(graph)
  ggraph::create_layout(graph, layout)
}

random_graph <- function(ids, n = 5, e = 4, layout = "fr") {
  df <- data.frame(id = ids, size = rnorm(n, .5, .2))
  edges <- data.frame(id1 = sample(ids, e), id2 = sample(ids, e),
                      width = rnorm(e, .5, .2))
  fast_layout(edges, layout, df)
}
