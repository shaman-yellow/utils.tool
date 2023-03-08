# ==========================================================================
# Used to shift a series of functions (from files) from one package to another.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

new_package.fromFiles <- function(pkg.path, files, path = NULL,
  depends = c("ggplot2", "grid"),
  exclude = c("MSnbase"),
  ...) {
  imports <- parent_packs(files, path)
  imports <- imports[ !imports %in% exclude ]
  new_package(pkg.path, imports, depends, ...)
  if (!is.null(path))
    files <- paste0(path, "/", files)
  file.copy(files, paste0(pkg.path, "/R"), T)
}

new_package <- function(path, imports = NULL, depends = NULL, extdata = T,
  fields = .new_package_fields(),
  gitignore_templ = .gitignore_templ()
  ) {
  pre.wd <- getwd()
  if (!file.exists(path)) {
    usethis::create_package(path, fields)
    usethis::use_mit_license()
    dir.create(paste0(path, "/inst/extdata"), recursive = T)
    file.copy(gitignore_templ, ".")
    writeLines(c(paste0("# ", stringr::str_extract(path, "[^/]*$")),
        "", "Under preparation...", ""), "README.md")
  } else {
    setwd(path)
  }
  cat("Now:", getwd(), "\n")
  if (!is.null(depends)) {
    lapply(depends, usethis::use_package, type = "depends")
  }
  if (!is.null(imports)) {
    lapply(imports, usethis::use_package)
  }
  setwd(pre.wd)
}

parent_packs <- function(files, path = NULL){
  if (!is.null(path))
    files <- paste0(path, "/", files)
  names <- lapply(files,
    function(file){
      file <- readLines(file)
      names <- grep_operater(file)
      names
    })
  unique(unlist(names))
}

grep_operater <- function(txt){
  packs <- stringr::str_extract(txt, "[a-zA-Z][a-zA-Z0-9._]{0,}(?=::)")
  packs <- packs[!vapply(packs, is.na, logical(1))]
  packs <- packs[vapply(packs, requireNamespace, logical(1), quietly = T)]
}

match_function <- function(txt){
  funs <- stringr::str_extract(txt, "[.a-zA-Z][a-zA-Z0-9._]{0,}(?=\\()")
  funs <- funs[!vapply(funs, is.na, logical(1))]
  funs <- unique(funs)
  parent <- lapply(funs, findFunction)
  parent <- parent[vapply(parent, function(p) if (length(p) > 0) T else F, logical(1))]
  parent <- lapply(parent, function(envs) {
    name <- lapply(envs,
      function(env) {
        attr(env, "name")
      })
    gsub("^.*:", "", unlist(name))
    })
  unique(unlist(parent))
}

.new_package_fields <- function() {
  author <- c(person(given = "LiChuang", family = "Huang",
      email = "shaman.yellow@foxmail.com",
      role = c("aut"),
      comment = c(ORCID = "0000-0002-5445-1988")),
    person(given = "Gang", family = "Cao",
      role = c("cre")))
  list(`Authors@R` = author, Author = author,
    Maintainer = "LiChuang Huang <shaman.yellow@foxmail.com>"
  )
}

.gitignore_templ <- function() {
  "~/MCnebula2/.gitignore"
}
