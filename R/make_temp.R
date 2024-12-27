# ==========================================================================
# Backing up temporary files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

backup.tmp <- function(tmp, backup_path = "~/tmp_backup/R") {
  if (!file.exists(backup_path))
    dir.create(backup_path, recursive = TRUE)
  if (!file.exists(tmp))
    stop("file.exists(tmp) == FALSE")
  file.copy(tmp, backup_path, TRUE, TRUE)
}
