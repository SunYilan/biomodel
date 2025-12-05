dir.create("../9.PTR_network")
setwd("./9.PTR_network")

#---------------------------------
rm(list = ls())

library(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(multipleROC)
library(mlr3viz)
library(future)

rm(list = ls())

sp_diff = read.csv("../4.masslin2_bacdiff/bac_diff_1229.csv",row.names = 1, header = T)
sp_cor = sp_diff

kegg_diff =  read.csv("../5.func_diff/func_diff_1229.csv",row.names = 1, header = T)
kegg_cor = kegg_diff

meta_diff =  read.csv("../8.PTR_diff/metab_diff_data.csv", header = T)

meta_diff = meta_diff[,1:203] %>% 
  group_by(feature) %>%
  slice_head(n = 1) %>%
  column_to_rownames("feature")

meta_diff[] <- lapply(meta_diff, function(x) as.numeric(x))
meta_diff_t <- t(meta_diff)
meta_cor <- as.data.frame(meta_diff_t)

meta_cor =  log(meta_cor+1)

meta_cor = meta_cor[rownames(sp_cor),]
kegg_cor = kegg_cor[rownames(sp_cor),]

cor_spearman=function(microbiome,metabolites,name, method = "spearman"){
  microbiome=microbiome[intersect(row.names(microbiome),row.names(metabolites)),]
  metabolites=metabolites[row.names(microbiome),]
  pvals = matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
  cors = matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
  for(i in 1:ncol(microbiome)){
    for(j in 1:ncol(metabolites)){
      a=microbiome[is.na(microbiome[,i])==F&is.na(metabolites[,j])==F,i]
      b=metabolites[is.na(microbiome[,i])==F&is.na(metabolites[,j])==F,j]
      if(length(a)>2&length(b)>2){
        cor = cor.test(a,b,method = method)
        pvals[j,i] = cor$p.value
        cors[j,i] = cor$estimate
      }
    }
  }
  qvals = matrix(p.adjust(pvals,method = "BH"),ncol = ncol(pvals))
  colnames(qvals) = colnames(microbiome)
  rownames(qvals) = colnames(metabolites)
  colnames(cors) = colnames(microbiome)
  rownames(cors) = colnames(metabolites)
  ind = which(pvals <=1,arr.ind = T)
  association = data.frame(metabolites = colnames(metabolites)[ind[,1]],microbiome = colnames(microbiome)[ind[,2]],Cor = cors[ind],Pval = pvals[ind],Qval = qvals[ind],stringsAsFactors = F)
  association$name=paste(association$metabolites,association$microbiome,sep = "_with_")
  colnames(association)[1:2]=name[c(2,1)]
  qvals[is.na(qvals)]=1
  cors=cors[row.names(qvals),colnames(qvals)]
  pvals[pvals<0.05]="."
  pvals[qvals<0.05]="*"
  pvals[pvals>0.05]=NA
  qvals[qvals<0.05]="*"
  qvals[qvals>0.05]=NA
  return(list(association,cors,qvals,pvals))
}


iv_fct = sp_cor

mv_fct = kegg_cor

dv_fct = meta_cor

iv_fct_name <- sub("_cor","",deparse(substitute(sp_cor)))
mv_fct_name <- sub("_cor","",deparse(substitute(kegg_cor)))
dv_fct_name <- sub("_cor","",deparse(substitute(meta_cor)))

cor_spearman(iv_fct,dv_fct,c(iv_fct_name ,dv_fct_name), "spearman") -> ccor_heatdat

ccor_heatdat[[1]] %>% 
  dplyr::select(meta,sp,Cor) %>%
  pivot_wider(meta,names_from = sp,values_from = Cor) %>%
  column_to_rownames("meta") -> plot_heatmap

ccor_heatdat[[1]] %>% 
  mutate(sign = case_when(
    Pval < 0.001 ~ "***",
    Pval < 0.01 ~ "**",
    Pval <= 0.05 ~ "*",
    Pval > 0.05 ~ ""
  )) %>%
  dplyr::select(meta,sp,sign) %>%
  pivot_wider(meta,names_from = sp,values_from = sign) %>%
  column_to_rownames("meta") -> plot_heatmap_text

pheatmap::pheatmap(plot_heatmap,
                   cellwidth = 10,
                   cellheight = 10,
                   display_numbers = plot_heatmap_text,
                   angle_col = 45,
                   filename = "UIA_sp_meta_cor_heatmap.pdf")
#---
cor_spearman(iv_fct,mv_fct,c(iv_fct_name ,mv_fct_name), "spearman") -> ccor1

cor_spearman(mv_fct,dv_fct,c(mv_fct_name,dv_fct_name),"spearman") -> ccor2

cor_spearman(iv_fct,dv_fct,c(iv_fct_name,dv_fct_name),"spearman") -> ccor3

ccor1[[1]] %>% filter(Pval < 0.05) -> ccor1_test
ccor2[[1]] %>% filter(Pval < 0.05) -> ccor2_test
ccor3[[1]] %>% filter(Pval < 0.05) -> ccor3_test

write.csv(ccor1_test,paste0(iv_fct_name,"-",mv_fct_name,"_spearman_Cor_p005.csv"))
write.csv(ccor2_test,paste0(mv_fct_name,"-",dv_fct_name,"_spearman_Cor_p005.csv"))

ccor1_test[,c(2,1,3,4)] %>% 
  dplyr::inner_join(ccor2_test[,c(2,1,3,4)],by = mv_fct_name) -> cor_pair_data

cor_pair_data %>% drop_na() -> cor_pair_data



cor_pair_data %>% mutate(updown.x = case_when(
  Cor.x >0 ~ "Postive",
  Cor.x <0 ~ "Negtive"
)) %>% mutate(updown.y = case_when(
  Cor.y >0 ~ "Postive",
  Cor.y <0 ~ "Negtive"
))  -> cor_pair_data


cor_pair_data %>% 
  filter(updown.x == "Postive" & updown.y == "Postive") %>%
  filter(Pval.x < 0.05 & Pval.y < 0.05) %>%
  mutate(genus = sub("_.*","",sp)) %>%
  mutate(keggsp = sub(".*__","",kegg)) %>%
  mutate(kegggenus = sub("_.*","",keggsp)) %>%
  filter(kegggenus == genus) %>%
  filter(Cor.x > 0.3 | Cor.y > 0.3) %>%
  filter(sp == keggsp) %>% 
  select(-c(genus,keggsp,kegggenus))-> cor_pair_data


colnames(cor_pair_data) = c(iv_fct_name,
                            mv_fct_name,
                            paste0("cor_",iv_fct_name,"_",mv_fct_name),
                            paste0("p_",iv_fct_name,"_",mv_fct_name),
                            dv_fct_name,
                            paste0("cor_",mv_fct_name,"_",dv_fct_name),
                            paste0("p_",mv_fct_name,"_",dv_fct_name),
                            paste0("updown_",iv_fct_name,"_",mv_fct_name),
                            paste0("updown_",mv_fct_name,"_",dv_fct_name)
)

write.csv(cor_pair_data,"cor_pair_data.csv")


colnames(cor_pair_data)

data_tmp1 = cor_pair_data[,c("sp","kegg")]
colnames(data_tmp1) = c("source","target")
data_tmp2 = cor_pair_data[,c("kegg","meta")]
colnames(data_tmp2) = c("source","target")

netdata_edge = rbind(data_tmp1,data_tmp2) %>% unique()


data_tmp1 = data.frame(cor_pair_data[["sp"]],"sp")
data_tmp2 = data.frame(cor_pair_data[["kegg"]],"kegg")
data_tmp3 = data.frame(cor_pair_data[["meta"]],"meta")

colnames(data_tmp1) = c("node","class")
colnames(data_tmp2) = c("node","class")
colnames(data_tmp3) = c("node","class")

netdata_node = rbind(data_tmp1,
                     data_tmp2,
                     data_tmp3) %>% unique()

write.csv(netdata_edge,"netdata_edge.csv")
write.csv(netdata_node,"netdata_node.csv")

data_tmp1 = cor_pair_data[,c("sp","kegg")]

table(data_tmp1$sp)


#000000000000000000
cor_pair_data %>% 
  group_by(
    !!sym(iv_fct_name),
    !!sym(mv_fct_name),
    !!sym(dv_fct_name),
    !!sym(paste0("updown_",iv_fct_name,"_",mv_fct_name)),
    !!sym(paste0("updown_",mv_fct_name,"_",dv_fct_name))) %>% 
  count() %>% filter(sp == "Fusobacterium_nucleatum")-> cor_pair_data_plot




cor_pair_data_plot$sp = sub("\\.t__.*","",cor_pair_data_plot$sp)
cor_pair_data_plot$sp = gsub(".*\\.","",cor_pair_data_plot$sp)
cor_pair_data_plot$kegg = gsub("_\\..*","",cor_pair_data_plot$kegg)
#cor_pair_data_plot = cor_pair_data_plot %>% filter(!grepl("s__",sp))

#cor_pair_data_plot %>% filter(sp ==  "Escherichia_coli") -> cor_pair_data_plot

write_tsv(cor_pair_data_plot,"chongjitu.data.txt")

dev.off()
require(ggalluvial)


colpal = readRDS("../0.colpal/colpal.rds")
ggplot(data = cor_pair_data_plot,
       aes(axis1 = !!sym(iv_fct_name), axis2 = !!sym(mv_fct_name), axis3 = !!sym(dv_fct_name),
           y = n)) +
  scale_x_discrete(limits = c(iv_fct_name, mv_fct_name, dv_fct_name), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = !!sym(mv_fct_name)), show.legend = F) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("#ff757a","#f2acbf",
                               "#efc7c7","#ffa394","#ffeace","#cbdeb8",
                               "#92c8ac","#3e658c","#c9687f","#4d316c","#f99417",
                               "#407d54","#89cdd5"))+
  theme_minimal() +
  ggtitle("")

ggsave("chonjitu1_RIA.pdf", width = 12, height = 6)

