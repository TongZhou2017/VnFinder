#' Call HMM model
#'
#' @description call hmm model in package
#' @param gene a character specifying the model target gene name, such as 18S,
#' 5.8S, ITS.
#' @param region a character specifying the model target region of gene, such as
#' v1-9 for 18S, 5.8S for 5.8S, ITS1 and ITS2 for ITS.
#' @param strand a character specifying the model target strand of gene region,
#' such as forward(fwd) for all 18S regions, ITS1 and 5.8S, type 1-3 for ITS2.
#' @return a character specifying the HMM model path
#' @export
#' @examples
#' model_call()
model_call <- function(gene="18S",region="v4",strand="fwd"){
  #gene = stringr::str_to_upper(gene)
  #region = stringr::str_to_upper(region)
  #strand = stringr::str_to_upper(strand)
  dir = system.file("extdata",package = "VnFinder")
  database <- data.frame(gene=c(rep("18S",18),"ITS","5.8S","ITS","ITS","ITS","ITS"),region=c(paste0("v",rep(c(1:9),each=2)),"ITS1","5.8S","ITS2","ITS2","ITS2","ITS2"),strand=c(rep(c("fwd","rev"),9),rep("fwd",3),c("type1","type2","type3")),model=c("18S_v1_fwd","18S_v1_rev","18S_v2_fwd","18S_v2_rev","18S_v3_fwd","18S_v3_rev","18S_v4_fwd","18S_v4_rev","18S_v5_fwd","18S_v5_rev","18S_v6_fwd","18S_v6_rev","18S_v7_fwd","18S_v7_rev","18S_v8_fwd","18S_v8_rev","18S_v9_fwd","18S_v9_rev","ITS_ITS1","5.8S","ITS_ITS2","ITS_ITS2_type1","ITS_ITS2_type2","ITS_ITS2_type3"),stringsAsFactors = FALSE)
  m <- database[which(database$gene == gene & database$region == region & database$strand == strand),]$model
  model <- paste0(dir,"/",m)
  return(model)
}

#' Extract Region by HMM
#'
#' @description Extract target region of gene by HMM model. The output is fasta
#' format sequences and txt format table.
#' @param seq a character specifying the template sequence file path.
#' @param gene_region a vector of characters specifying the gene name and region
#' name in a pasted string. This parameter usually defined by the range_scroller
#' function.
#' @param extend a logical specifying the extract method is or not extend to the
#' ending.
#' @return No return value, only output files
#' @importFrom itol.toolkit file_get_dir
#' @importFrom itol.toolkit file_get_name
#' @importFrom itol.toolkit fa_read
#' @importFrom itol.toolkit fa_write
#' @import dplyr
#' @importFrom stringr str_sub
#' @importFrom utils read.table
#' @export
region_extract_hmm <- function(seq,gene_region,extend=FALSE){
  list_gene_region <- strsplit(gene_region,"_")
  list_gene_region <- unlist(list_gene_region[[1]])
  gene = list_gene_region[1]
  region = list_gene_region[2]
  strands <- "fwd"
  if(list_gene_region[2] == "ITS2"){
    strands <- c("fwd","type1","type2","type3")
    models <- c(model_call(gene,region,strands[1]),
                model_call(gene,region,strands[2]),
                model_call(gene,region,strands[3]),
                model_call(gene,region,strands[4]))
  }else{
    models = model_call(gene,region,strands)
  }
  dir=itol.toolkit::file_get_dir(seq)
  seq_name=itol.toolkit::file_get_name(seq,with_ext = FALSE)
  for (model in models) {
    model_name=itol.toolkit::file_get_name(model)
    output = paste0(dir,"/",seq_name,"_by_",model_name,".tab")
    output_fa = paste0(dir,"/",seq_name,"_by_",model_name,".fa")
    # state file check point  # nhmmer
    if(file.exists(output)){
      message("nhmmer result already exist")
    }else{
      message("nhmmer result not exist")
      system(paste0("nhmmer --noali --tblout ",output," ",model," ",seq),ignore.stdout = T)
    }
    nhmmer <- read.table(output,stringsAsFactors = FALSE)
    nhmmer <- nhmmer %>% distinct(V1,.keep_all = T)
    names(nhmmer) <- c("seq_name", "target_accession", "query_name", "query_accession", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to","sq_len", "strand", "evalue", "score", "bias", "description")
    nhmmer <- nhmmer[,1:16]
    fa <- itol.toolkit::fa_read(seq)
    #dim(fa)
    cat(length(intersect(fa$seq_name,nhmmer$seq_name)),"/",nrow(fa))
    fa_nhmmer <- dplyr::inner_join(fa,nhmmer,by="seq_name")
    if(extend){
      fa_nhmmer_seq <- fa_nhmmer %>% rowwise() %>% mutate(seq = stringr::str_sub(sequence,min(env_from,env_to),-1)) %>% select(seq_name,seq)
    }else{
      fa_nhmmer_seq <- fa_nhmmer %>% rowwise() %>% mutate(seq = stringr::str_sub(sequence,min(env_from,env_to),max(env_from,env_to))) %>% select(seq_name,seq)
    }

    names(fa_nhmmer_seq)[2] <- "sequence"

    # write
    itol.toolkit::fa_write(fa_nhmmer_seq,output_fa)
  }
  if(length(models) == 4){
    type_0 <- paste0(seq_name,"_by_ITS_ITS2.fa")
    type_1 <- paste0(seq_name,"_by_ITS_ITS2_type1.fa")
    type_2 <- paste0(seq_name,"_by_ITS_ITS2_type2.fa")
    type_3 <- paste0(seq_name,"_by_ITS_ITS2_type3.fa")
    seq_0 <- itol.toolkit::fa_read(paste0(dir,"/",type_0))
    seq_1 <- itol.toolkit::fa_read(paste0(dir,"/",type_1))
    seq_2 <- itol.toolkit::fa_read(paste0(dir,"/",type_2))
    seq_3 <- itol.toolkit::fa_read(paste0(dir,"/",type_3))
    seq <- dplyr::full_join(seq_0,seq_1,by="seq_name",suffix = c(".type_0",".type_1")) %>% dplyr::full_join(seq_2,by="seq_name") %>% dplyr::full_join(seq_3,by="seq_name",suffix = c(".type_2",".type_3"))
    seq$sequence <- apply(seq[,2:ncol(seq)],1, function(x) x[which.max(nchar(x))])
    seq <- seq %>% select(seq_name,sequence)
    output_seq = paste0(dir,"/",seq_name,"_by_ITS_ITS2_types.fa")
    itol.toolkit::fa_write(seq,output_seq)
  }
}
