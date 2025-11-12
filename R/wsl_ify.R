wsl_ify <- function(file_path){
  if(grepl("WSL2", Sys.info()[2]) & substr(file_path, 1, 3) == "C:/"){
    return(gsub("C:/", "/mnt/c/", file_path))
  }
  if(grepl("Windows", Sys.info()[1]) & substr(file_path, 1, 4) == "/mnt"){
    return(gsub("/mnt/c/", "C:/",  file_path))
  }
  return(file_path)
}
