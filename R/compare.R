#' Evaluate Correction
#'
#' @description Using NMDS and procrustes analysis to evaluate the region effect
#' .
#' @param otutab a character specifying the otutab file or a data.frame of
#' otutab. The first column should be ID. The first row should be header.
#' @param set1 a character specifying the first dataset to compare.
#' @param set2 a character specifying the second dataset to compare.
#' @return a list contains p-value, mean, procrustes plot, bar plot
#' @importFrom data.table fread
#' @import dplyr
#' @import vegan
#' @import ggplot2
#' @importFrom raster pointDistance
#' @importFrom cluster pam
#' @export
evalue_change <- function(otutab,set1,set2){
  tab <- data.table::fread(otutab)
  tab_raw <- tab %>% select(!ends_with("HMM")) %>% select(-ID)
  tab_raw <- as.data.frame(t(tab_raw))
  colnames(tab_raw) <- tab$ID
  tab_raw <- tab_raw[grep(paste0("^",set1,'|',set2),rownames(tab_raw)),]
  tab_hmm <- tab %>% select(ends_with("HMM"))
  tab_hmm <- as.data.frame(t(tab_hmm))
  colnames(tab_hmm) <- tab$ID
  tab_hmm <- tab_hmm[grep(paste0("^",set1,'|',set2),rownames(tab_hmm)),]
  rename_tab <- data.table::fread("/Users/zhoutong/Downloads/docker_volumn/mee/script/rename.sh",header = F)
  rename_tab <- rename_tab %>% mutate(old_name = stringr::str_remove(V1,"/.*"))
  rename_tab_2 <- data.table::fread("/Users/zhoutong/Downloads/docker_volumn/mee/script/rename.txt",header = F)
  names(rename_tab_2) <- c("old_name","new_name")
  rename_tab <- left_join(rename_tab,rename_tab_2,by="old_name")
  rename_tab$sample_id <- c(rep(1:12,each=2),rep(1:73,each=2),rep(1:23,each=2),rep(1:33,each=2),rep(1:32,each=2),rep(1:17,each=2),rep(1:44,each=2),rep(1:119,each=2),rep(1:60,each=2),rep(1:149,each=2),rep(1:16,each=2))
  rename_tab$type <- rep(c("RAW","HMM"),1156/2)
  rename_tab <- rename_tab %>% mutate(label = paste0(new_name,"_",sample_id,"_",type))
  group_raw <- rename_tab %>% filter(type == "RAW", new_name %in% c(set1,set2))
  #nmds
  dist_raw <- vegdist(tab_raw,index="bray")
  dist_hmm <- vegdist(tab_hmm,index="bray")
  mds_raw <- suppressMessages(metaMDS(dist_raw))
  mds_hmm <- suppressMessages(metaMDS(dist_hmm))
  score_raw <- as.data.frame(scores(mds_raw))
  score_hmm <- as.data.frame(scores(mds_hmm))
  mds_data_raw <- as.data.frame(mds_raw$points)
  mds_data_hmm <- as.data.frame(mds_hmm$points)
  mds_data_raw$SampleID <- rownames(mds_data_raw)
  mds_data_hmm$SampleID <- rownames(mds_data_hmm)
  mds_data_raw <- dplyr::left_join(mds_data_raw, rename_tab,by=c("SampleID" = "label"))
  mds_data_hmm <- dplyr::left_join(mds_data_hmm, rename_tab,by=c("SampleID" = "label"))
  proc <- procrustes(X = score_raw, Y = score_hmm, symmetric = TRUE)
  result <- data.frame(residuals = residuals(proc),group_raw)
  Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X),result)
  X <- data.frame(proc$rotation)
  p <- ggplot(Y) +
    geom_vline(xintercept = pam(data.frame(Y$X1, Y$X2),1)$medoids[1], color = 'gray', linetype = 2, size = 0.3) +
    geom_hline(yintercept = pam(data.frame(Y$X1, Y$X2),1)$medoids[2], color = 'gray', linetype = 2, size = 0.3) +
    geom_segment(aes(x = NMDS1, y = NMDS2, xend = X1, yend = X2), arrow = arrow(length = unit(0.2, 'cm')), color = 'black',alpha=0.3, size = 0.3) +
    geom_point(aes(X1, X2, color=new_name),alpha=0.7, size=2,shape = 1) +
    geom_point(aes(NMDS1, NMDS2, color=new_name),alpha=0.7, size=2, shape = 16) +
    geom_point(aes(x=pam(data.frame(Y$X1, Y$X2),1)$medoids[1],y=pam(data.frame(Y$X1, Y$X2),1)$medoids[2]),colour="red",shape=10,size=2)+
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.key = element_rect(fill = 'transparent')) +
    labs(x = 'NMDS1', y = 'NMDS2', color = '')
  Y$change = NA
  for (i in 1:nrow(Y)) {
    Y$change[i] <- pointDistance(c(pam(data.frame(Y$X1, Y$X2),1)$medoids[1],pam(data.frame(Y$X1, Y$X2),1)$medoids[2]), c(Y$X1[i],Y$X2[i]), lonlat=FALSE) - pointDistance(c(pam(data.frame(Y$X1, Y$X2),1)$medoids[1],pam(data.frame(Y$X1, Y$X2),1)$medoids[2]), c(Y$NMDS1[i],Y$NMDS2[i]), lonlat=FALSE)
  }
  t <- t.test(Y$change,mu=0)
  return(list(t$estimate,t$p.value,p,sort(Y$change),proc))
}

