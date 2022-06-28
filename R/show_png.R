show_png <- 
  function(
           file,
           size = "500x"
           ){
    if(grepl("\\.svg$", file)){
      tofile <- sub("\\.svg$", ".png", file)
      if(!file.exists(tofile)){
        system(
               paste("cairosvg -i", file, "-d 300", "-o", tofile)
        )
      }
      file <- tofile
    }
    png <- magick::image_read(file)
    png <- magick::image_resize(png, size)
    grid::grid.raster(png)
  }

