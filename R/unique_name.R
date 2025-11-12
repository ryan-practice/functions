unique_name <- function(base_name){
  if(file.exists(base_name)){
    extension <- substr(base_name, max(str_locate_all(base_name, "\\.")[[1]]), nchar(base_name))
    copy_suffix <- 1
    copy_suffix <- sprintf("%02d", copy_suffix)
    base_name <- gsub(extension, paste0("_", copy_suffix, extension), base_name)
    while(file.exists(base_name)){
      copy_suffix <- as.numeric(copy_suffix) + 1
      copy_suffix <- sprintf("%02d", copy_suffix)
      base_name <- gsub(paste0("_\\d{2}",extension), paste0("_", copy_suffix, extension), base_name)
    }
  }
  return(base_name)
}
