# ==========================================================================
# Backing up temporary files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

backup.tmp <- function(tmp, backup_path = "~/tmp_backup/R") {
  if (!file.exists(backup_path))
    dir.create(backup_path, recursive = T)
  if (!file.exists(tmp))
    stop("file.exists(tmp) == F")
  file.copy(tmp, backup_path, T, T)
}
