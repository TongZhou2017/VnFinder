#' Usearch wrapper
#'
#' @description Wrapper function for usearch. Check host local first, then check
#' docker container environment. Should setup the option parameter by the guide
#' in document.
#' @param command a character string specifying the usearch command, which start
#' with 'usearch'. In docker environment, single quote is defualt for command
#' transform. Double quote should be avoid to be used in command.
#' @return stdout screen from usearch
#' @export
#' @examples
#' options(usearch.docker.path = "docker ps -aqf 'name=usearch' | \\
#'                                xargs -I {} docker exec -i {} /bin/sh -c ")
#' options(usearch.path = '/bin/')
#' run_usearch(command = "usearch -fastx_get_sample_names /mnt/Octocarol.fa \\
#'                                -output /mnt/Octocarol_samples.txt")
run_usearch <- function(command){
  if(is.null(getOption("usearch.path"))){
    stop("Please setup usearch environment parameters by 'Setup' section in document")
  }
  if(file.exists(paste0(getOption("usearch.path"),"usearch"))){
    system(command = paste0(getOption("usearch.path"),command),ignore.stdout = T)
  }else{
    str_quote = "'"
    if(grepl("\"",command)){
      str_quote = "'"
    }
    if(grepl("'",command)){
      str_quote = '"'
    }
    system(command = paste0(getOption("usearch.docker.path"),str_quote,command,str_quote),ignore.stdout = T)
  }
}
