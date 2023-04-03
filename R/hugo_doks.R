
new_scene <- 
  function(scene, weight = rep(100, length(scene)), ex_weight = weight,
    path = "~/MCnebula2", suffix = "/content/en", tar = "docs", index_Draft = T){
    if (!is.vector(scene)) {
      stop("is.vector(scene) == F")
    }
    if (!is.character(names(scene))){
      stop( "is.character(names(scene)) == F" )
    }
    path <- paste0(path, suffix)
    path <- target_dir(path, tar)
    lapply(1:length(scene),
      function(n){
        name <- names(scene)[n]
        dir <- paste0(path, "/", name)
        if (!dir.exists(dir)) {
          dir.create(dir)
          if (index_Draft) {
            index <- paste0(dir, "/_index.Rmd")
            writeLines(index_Draft(name, weight[[n]]), index)
          }
        }
        file <- paste0(dir, "/", scene[[n]], ".Rmd")
        if (!file.exists(file)) {
          writeLines(Draft(scene[[n]], ex_weight[[n]], tar, name), file)
        }
      })
    message("Done")
  }

target_dir <- 
  function(path, tar){
    lst <- list.dirs(path, recursive = T)
    lst <- lst[grepl(paste0("(?<=/)", tar, "$"), lst, perl = T)][1]
    lst
  }

record_time <- function(){
  format(Sys.time(), "%Y %b %e %H:%M:%S | %a")
}

Draft <- function(title, weight = 100, tar, name){
  title <- gsub("_", " ", title)
  c("---",
    "contributors:\n- LiChuang Huang",
    paste0("title: ", "\"", Hmisc::capitalize(title), "\""),
    paste0("date: ", "\"", record_time(), "\""),
    paste0("lastmod: ", "\"", record_time(), "\""),
    "draft: false",
    "images: []",
    "menu:",
    strwrap(paste0(tar, ":"), indent = 2),
    paste0(strwrap("parent:", indent = 4), " \"", name, "\""),
    paste0("weight: ", weight),
    "toc: true",
    "---"
  )
}

index_Draft <- function(title, weight = 100){
  title <- gsub("_", " ", title)
  c("---",
    paste0("title: ", "\"", Hmisc::capitalize(title), "\""),
    paste0("date: ", "\"", format(Sys.time(), record_time()), "\""),
    paste0("lastmod: ", "\"", format(Sys.time(), record_time()), "\""),
    "draft: false",
    "images: []",
    paste0("weight: ", weight),
    "toc: true",
    "---"
  )
}

target_file <- 
  function(path, tar){
    lst <- list.files(path, recursive = T, full.names = T)
    lst <- lst[grepl(paste0("(?<=/)", tar, "$"), lst, perl = T)][1]
    lst
  }

setGeneric("set_home", 
  function(x, ...) standardGeneric("set_home"))
setMethod("set_home", 
  signature = setMissing("set_home"),
  function(){
    function(path = "~/MCnebula2/config", tar = "params.toml") {
      target_file(path, tar)
    }
  })

setMethod("set_home", 
  signature = setMissing("set_home",
    x = "vector"),
  function(x, ...){
    path <- set_home()(...)
    lines <- readLines(path)
    for (i in names(x)) {
      lines <- repl_huto(i, x[[i]], lines)
    }
    writeLines(lines, path)
  })

repl_huto <-
  function(key, content, lines,
    left = "\"", right = "\"", link = " = "){
    pattern <- paste0("^\\s{0,}", key, "\\s{0,}(?![a-z|A-Z|0-9|_|.])")
    n <- grep(pattern, lines, perl = T)
    pattern <- paste0("(?<=", link, left, ")", "[^\"]{1,}", "(?=", right, ")")
    lines[n] <- gsub(pattern, content, lines[n], perl = T)
    lines
  }
repl_yaml <- function(key, content, lines){
  repl_huto(key, content, lines, link = ": |: \"",
    left = "", right = "\"|$"
  )
}

setGeneric("set_index", 
  function(x, tar, ...) standardGeneric("set_index"))
setMethod("set_index", 
  signature = setMissing("set_index"),
  function(){
    function(tar = "en/_index.Rmd", path = "~/MCnebula2/content/en") {
      target_file(path, tar)
    }
  })

setMethod("set_index", 
  signature = setMissing("set_index",
    x = "vector",
    tar = "character"),
  function(x, tar, ...){
    path <- set_index()(tar, ...)
    lines <- readLines(path)
    time <- grep("^lastmod|^date", lines)
    for (i in time) {
      if (!grepl(": \"", lines[i]))
        lines[i] <- sub(": ", ": \"", lines[i])
      if (!grepl("\"$", lines[i]))
        lines[i] <- sub("$", "\"", lines[i])
    }
    for (i in names(x)) {
      lines <- repl_yaml(i, x[[i]], lines)
    }
    writeLines(lines, path)
  })

inclu.fig.ht <- function(src, to, caption = "This is a figure",
  width = "100%", height = NULL,
  rel.path = "~/MCnebula2/content/en")
{
  if (!file.exists(paste0(rel.path, to))) {
    if (!file.exists(src)) {
      stop("file.exists(", src, ") == F")
    }
    if (grepl("\\.pdf$", src)) {
      system(paste0("pdf2svg ", src, " ", rel.path, to, " 1"))
    } else {
      file.copy(src, paste0(rel.path, to))
    }
  }
  draft <- c("<figure>", "", "", "</figure>")
  img <- paste0("<center>", "<img src=\"", to, "\"", ">", "</center>")
  cap <- paste0("<center>", "<figcaption>", caption, "</figcaption>", "</center>")
  draft[2:3] <- c(img, cap)
  writeLines(draft)
}
