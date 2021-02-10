remove_white <- function(x) gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)

to_Bytes <- function(value) {
  num <- as.numeric(value[1])
  # Avoid case-sensitive
  units <- tolower(value[2])
  lookup <- list(
    "kb" = "kilobytes",
    "mb" = "megabytes",
    "gb" = "gigabytes",
    "tb" = "terabytes"
  )

  if (units %in% lookup) {
    power <- which(units == lookup)
    result <- num * 1024^power
  }
  else if (units %in% names(lookup)) {
    power <- which(units == names(lookup))
    result <- num * 1024^power
  }
  else {
    result <- num
  }

  return(result)
}

system_ram <- function(os) {
  if (length(grep("^linux", os))) {
    cmd <- "awk '/MemTotal/ {print $2}' /proc/meminfo"
    ram <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    ram <- as.numeric(ram) * 1024
  }
  else if (length(grep("^darwin", os))) {
    ram <- system("system_profiler -detailLevel mini | grep \"  Memory:\"",
      intern = TRUE,
      ignore.stderr = TRUE
    )[1]
    ram <- remove_white(ram)
    ram <- to_Bytes(unlist(strsplit(ram, " "))[2:3])
  }
  else if (length(grep("^solaris", os))) {
    cmd <- "prtconf | grep Memory"
    ram <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    ram <- remove_white(ram)
    ram <- to_Bytes(unlist(strsplit(ram, " "))[3:4])
  }
  else {
    ram <- system("wmic MemoryChip get Capacity", intern = TRUE)[-1]
    ram <- remove_white(ram)
    ram <- ram[nchar(ram) > 0]
    ram <- sum(as.numeric(ram))
  }
  return(as.double(ram))
}

ramSize <- function() {
  os <- R.version$os
  ram <- try(system_ram(os), silent = TRUE)
  if (class(ram) == "try-error") {
    ram <- 2048 * 1024 * 1024
  }

  return(ram[1])
}

ldet <- function(m) spam::determinant(m, logarithm = TRUE)$modulus

DIST <- fields::rdist
SQLdbList <- filehashSQLite::dbList
log <- base::log
diag.spam <- spam::diag.spam
rlimit <- ramSize()
