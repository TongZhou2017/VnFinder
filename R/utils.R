#' Gene region scroller
#'
#' @description Select gene region by a vector of index
#' @param index a vector of numbers specifying the gene regions. 1-18S_v1,
#' 2-18S_v2, 3-18S_v3, 4-18S_v4, 5-18S_v5, 6-18S_v6, 7-18S_v7, 8-18S_v8,
#' 9-18S_v9, 10-ITS_ITS1, 11-5.8S_5.8S, 12-ITS_ITS2.
#' @return a vector of characters specifying the pasted string of gene and
#' region name
#' @importFrom dplyr case_when
#' @export
#' @examples
#' # single region
#' range_scroller(4)
#' # a range of regions
#' range_scroller(3:6)
range_scroller <- function(index){
  gene <- c()
  region <- c()
  gene_region <- c()
  for (i in 1:length(index)) {
    gene[i] <- dplyr::case_when(index[i] %in% 1:9 ~ "18S",
                                index[i] %in% c(10,12) ~ "ITS",
                                index[i] == 11 ~ "5.8S")
    region[i] <- dplyr::case_when(index[i] %in% 1:9 ~ paste0("v",index[i]),
                                  index[i] == 10 ~ "ITS1",
                                  index[i] == 11 ~ "5.8S",
                                  index[i] == 12 ~ "ITS2")
    gene_region[i] <- paste0(gene[i],"_",region[i])
  }
  return(gene_region)
}

#' Unique fasta file
#'
#' @description Unique fasta file by label and sequence
#' @param file a character specifying the fasta file path
#' @param dir a character specifying the output dir
#' @return No return value, only output files
#' @export
fa_unique <- function(file,dir){
  seq <- itol.toolkit::fa_read(file)
  file_name=itol.toolkit::file_get_name(file,with_ext = TRUE,keep_dir = FALSE)
  seq <- seq %>% distinct(seq_name,sequence)
  output = paste0(dir,"/",file_name)
  fa_write(seq,output)
}
