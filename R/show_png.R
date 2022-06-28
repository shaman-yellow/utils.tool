show_png <- 
  function(
           file,
           path = "mcnebula_results",
           size = "500x"
           ){
    if(grepl("\\.svg$", file)){
      tofile <- sub("\\.svg$", ".png", file)
      if(!file.exists(tofile)){
        system(
               paste0("cairosvg -i ", path, "/", file, " -d 300", " -o ", path, "/", tofile)
        )
      }
      file <- tofile
    }
    png <- magick::image_read(file)
    png <- magick::image_resize(png, size)
    grid::grid.raster(png)
  }

