png_gather_two <- 
  function(
           png_file1,
           png_file2 = NA,
           width = 5000,
           height = 3500,
           internal = 0.06
           ){
    ## read png1
    png1 <- png::readPNG(png_file1)
    ratio.1 <- ncol(png1) / nrow(png1)
    ## read png2
    if(!is.na(png_file2)){
      png2 <- png::readPNG(png_file2)
      ratio.2 <- ncol(png2) / nrow(png2)
      ## filename fix
      fix <- "gather_"
      ## png postion shift
      position_shift <- "0"
    }else{
      ratio.2 <- NULL
      ## filename fix
      fix <- "ps_"
      ## png postion shift
      position_shift <- "(unit.width / 2 + internal / 2) / width.adjust"
    }
    ## ------------------ 
    max <- max(c(ratio.1, ratio.2))
    ## according to height, calculate needed width
    expect.width <- height * max / (1 - internal) * 2
    ## if width not enough
    if(width < expect.width){
      width <- expect.width
      width.adjust <- 1
    }else{
      width.adjust <- expect.width / width
    }
    ## save path and savename
    path <- get_path(png_file1)
    name <- get_filename(png_file1)
    png(paste0(path, "/", fix, name), width = width, height = height)
    plot.new()
    ## x, y, xend, yend
    unit.width <- (0.5 - internal / 2) * width.adjust
    rasterImage(png1,
                xleft = unit.width * (1 - ratio.1 / max) +
                  eval(parse(text = position_shift)),
                ybottom = 0,
                xright = unit.width + 
                  eval(parse(text = position_shift)),
                ytop = 1)
    if(!is.na(png_file2)){
      ## x, y, xend, yend
      rasterImage(png2,
                  xleft = unit.width + internal,
                  ybottom = 0,
                  xright = unit.width * (1 + ratio.2 / max) + internal, 
                  ytop = 1)
    }
    dev.off()
  }
get_path <- 
  function(
           path_str
           ){
    path <- stringr::str_extract(path_str, ".*(?=/)")
    return(path)
  }
get_filename <- 
  function(
           path_str
           ){
    filename <- stringr::str_extract(path_str, "(?<=/)[^/]*$")
    return(filename)
  }
## ---------------------------------------------------------------------- 
png_to_pdf <- 
  function(
           file
           ){
    file.pdf <- gsub("\\.png", "\\.pdf", file)
    pdf(file = file.pdf)
    png <- EBImage::readImage(file)
    EBImage::display(png, method = "raster", all = T)
    dev.off()
  }
