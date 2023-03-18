#' Primer search
#'
#' @description Using usearch -search_pcr2 function to extract hypervariable
#' region based on primer alignment method.
#' @param keywords a vector of characters specifying the query words in any
#' order. But the letter case should be considered, such like using 18S not 18s,
#' using v4 not V4. The short words for direction are format strict, using fwd
#' for forward primer and rev for reverse primer.
#' @param target a vector of characters specifying the query result target
#' information. The character should use the builtin_primers list elements name,
#' such as gene, region, set, name, id, direction, specificity, seq,
#' start_yeast, reference. If the length of vector is 1, the result will be a
#' character. If the length of vector is more than 1, the result will be a list.
#' @param ignore a cector of characters specifying the ignored query range. The
#' character should use the builtin_primers list elements name, such as gene,
#' region, set, name, id, direction, specificity, seq, start_yeast, reference.
#' @return a character or a vector of characters specifying the query target
#' information
#' @export
#' @examples
#' primer_v4_f = primer_search(keywords = c("18S","v4","fwd"))
primer_search <- function(keywords,target="seq",ignore=NULL){
  query_item <- names(builtin_primers[[1]])
  query_item <- query_item[which(!query_item %in% ignore)]
  query_keywords <- data.frame(matrix(data = FALSE,
                                      nrow = length(keywords),
                                      ncol = length(builtin_primers)))
  for (i in 1:length(keywords)) {
    query_keyword <- data.frame(matrix(data = FALSE,
                                       nrow = length(query_item),
                                       ncol = length(builtin_primers)))
    for (j in 1:length(query_item)) {
      query_keyword[j,] <- sapply(builtin_primers,
                                  "[[",
                                  query_item[j]) == keywords[i]
    }
    query_keywords[i,] <- apply(query_keyword, 2, any)
  }
  query_indexs <- which(colSums(query_keywords) == max(colSums(query_keywords)))
  query_results <- builtin_primers[query_indexs]
  query_targets <- list()
  for (k in 1:length(query_indexs)) {
    query_targets[[k]] <- builtin_primers[[query_indexs[k]]][target]
  }
  if(length(query_indexs) == 1 && length(target == 1)){
    query_targets <- unlist(query_targets)
  }
  return(query_targets)
}

#' Usearch search_pcr2
#'
#' @description Using usearch -search_pcr2 function to extract gene region by
#' primer alignment method.
#' @param template a character specifying the template sequence file in fasta or
#' fastq format.
#' @param gene_region a vector of characters specifying the gene name and region
#' name in a pasted string. This parameter usually defined by the range_scroller
#' function.
#' @param fwd a character specifying the forward primer set name.
#' @param rev a character specifying the reverse primer set name.
#' @return usearch output files in docker volumn
#' @importFrom itol.toolkit file_get_name
#' @export
usearch_pcr2 <- function(template, gene_region,fwd="default",rev="default"){
  template_name <- itol.toolkit::file_get_name(template,with_ext = F,keep_dir = F)
  list_gene_region <- strsplit(gene_region,"_")
  n <- length(list_gene_region)
  primer_f <- primer_search(keywords = c(list_gene_region[[1]][1],list_gene_region[[1]][2],"fwd",fwd))
  primer_r <- primer_search(keywords = c(list_gene_region[[n]][1],list_gene_region[[n]][2],"rev",rev))
  if(n == 1){
    command <- paste0("usearch -search_pcr2 ", template, " -fwdprimer ", primer_f, " -revprimer ", primer_r, " -strand both -fastaout /mnt/", template_name, "_", list_gene_region[[1]][1],"_",list_gene_region[[1]][2], ".fa")
  }else{
    command <- paste0("usearch -search_pcr2 ", template, " -fwdprimer ", primer_f, " -revprimer ", primer_r, " -strand both -fastaout /mnt/", template_name, "_", list_gene_region[[1]][1],"_",list_gene_region[[1]][2], "_to_", list_gene_region[[n]][1], "_",list_gene_region[[n]][2], ".fa")
  }
  run_usearch(command = command)
}

#' Theoretical length of amplicon
#'
#' @description Query theoretical length of amplicon based on yeast location of
#' primer.
#' @param index a vector of numbers specifying the gene regions. 1-18S_v1,
#' 2-18S_v2, 3-18S_v3, 4-18S_v4, 5-18S_v5, 6-18S_v6, 7-18S_v7, 8-18S_v8,
#' 9-18S_v9, 10-ITS_ITS1, 11-5.8S_5.8S, 12-ITS_ITS2.
#' @return a number specifying the theoretical length of amplicon
#' @export
#' @examples
#' theor_length(1:9)
theor_length <- function(index){
  gene_region <- range_scroller(index)
  list_gene_region <- strsplit(gene_region,"_")
  n <- length(list_gene_region)
  as.numeric(primer_search(c(list_gene_region[[n]],"rev"),"start_yeast") -
               primer_search(c(list_gene_region[[1]],"fwd"),"start_yeast") +
               nchar(primer_search(c(list_gene_region[[n]],"rev"),"seq")))
}
