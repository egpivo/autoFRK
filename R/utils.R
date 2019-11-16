system_ram <- function(os) {
  remove_white <- function(x) gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
  
  to_Bytes <- function(value) {
    num <- as.numeric(value[1])
    units <- value[2]
    power <- match(units, c("kB", "MB", "GB", "TB"))
    
    if (!is.na(power)) return(num * 1024^power)
    
    power <- match(units, c("Kilobytes",
                            "Megabytes",
                            "Gigabytes", 
                            "Terabytes"))
    
    if (!is.na(power)) 
      return(num * 1024^power)
    else 
      return(num)
  }
  
  if (length(grep("^linux", os))) {
    cmd <- "awk '/MemTotal/ {print $2}' /proc/meminfo"
    ram <- system(cmd, intern = TRUE)
    ram <- as.numeric(ram) * 1024
  }
  else if (length(grep("^darwin", os))) {
    ram <- system("system_profiler -detailLevel mini | grep \"  Memory:\"", 
                  intern = TRUE)[1]
    ram <- remove_white(ram)
    ram <- to_Bytes(unlist(strsplit(ram, " "))[2:3])
  }
  else if (length(grep("^solaris", os))) {
    cmd <- "prtconf | grep Memory"
    ram <- system(cmd, intern = TRUE)
    ram <- remove_white(ram)
    ram <- to_Bytes(unlist(strsplit(ram, " "))[3:4])
  }
  else {
    ram <- system("wmic MemoryChip get Capacity", intern = TRUE)[-1]
    ram <- remove_white(ram)
    ram <- ram[nchar(ram) > 0]
    sum(as.numeric(ram))
  }
  as.double(ram)
}

ramSize <- function() {
  os <- R.version$os
  ram <- try(system_ram(os), silent = TRUE)
  if (class(ram) == "try-error") 
    ram = 2048 * 1024 * 1024
  
  return(ram[1])
}

ldet <- function(m) spam::determinant(m, logarithm = TRUE)$modulus

DIST <- fields::rdist
SQLdbList <- filehashSQLite::dbList
log <- base::log
diag.spam <- spam::diag.spam

