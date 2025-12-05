setwd("./6.GCMS_diff")

require(tidyverse)
rm(list = ls())

# get metadata---------
alldata = read.csv("../1.data_tidy/OSCC_M.sum.1212.csv")

# get gcmsdata --------
gdata =  read.csv("../0.Rawdata/4.GCMS/OSCC_vs_HC_info_v1.csv",row.names = 1)

intersect(alldata$Sample_ID,colnames(gdata)[7:228])
setdiff(alldata$Sample_ID,colnames(gdata)[7:228])
setdiff(colnames(gdata)[7:228],alldata$Sample_ID)




gmeta = gdata[,1:6]

pie_data = as.data.frame(table(gmeta$Class.I)) %>% 
  arrange(rev(Freq))
pie_data$Var1 = factor(pie_data$Var1,levels = pie_data$Var1)
pie_data %>% mutate(
  ratio = round(Freq / sum(Freq) * 100, 1),
  Compounds = paste0(Var1," (",Freq,")")
) -> pie_data


palette =  readRDS("../0.colpal/colpal.rds")

ggplot(pie_data, aes(x = "", y = Freq, fill = Compounds)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Classification of Metabolites in Human Saliva by GC-MS") +
  scale_fill_manual(values = c("#c9687f", "#92c8ac", "#ff757a", "#f2acbf", "#efc7c7", "#ffa394",
                               "#ffeace", "#cbdeb8", "#92c8ac", "#3e658c", "#89cdd5"))+
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = paste0(round(Freq / sum(Freq) * 100, 1), "%")), 
            position = position_stack(vjust = 0.5))
ggsave("class_I_pie.pdf",width = 8,height = 8)


# mname = read.csv("../0.Rawdata/5.PTR/metab_name.csv")
# 
# 
# mname %>% left_join(gmeta,by = c("compound" = "Compounds")) %>%
#   select(compound,`Class.I`)


meta = read.csv("../1.data_tidy/meta_info_fix.csv")

diffset =  gdata[,7:228] %>% 
  rownames_to_column("feature") %>% 
  sjmisc::rotate_df(cn =T,rn = "Sample_ID") %>%
  column_to_rownames("Sample_ID")
diifmeta = meta %>%  
  filter(!Category == "outset") %>%
  column_to_rownames("Sample_ID")

otu_data2 = diffset[rownames(diifmeta),]
phe = diifmeta

otu_exp1 = otu_data2
set.seed(666)
source("../addin.r")


meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")

s_id = attributes(arg_distance)$Labels
results=NULL

for(i in colnames(phe)[-c(1:3)]){
  adonis_re <- tryCatch(
    {
      adonis3(arg_distance ~ phe[s_id, ][[i]], permutations = 1000, na.action = na.omit)
    },
    error = function(e) {
      message(paste("Error occurred in adonis3 function for column", i, ":", e))
      return(NULL)
    }
  )
  
  if (!is.null(adonis_re)) {
    summary(adonis_re)
    rownames(adonis_re)[1] = i
    results <- rbind(results, adonis_re[1, ])
  }
}


#what's mean of R2 of adonis
#https://www.researchgate.net/post/When_do_we_reject_null_hypothesis_of_not_difference_between_groups_in_MRPP_and_Adonis

adonis_data = results

adonis_data$R2_2 = round(adonis_data$R2, digits = 4)

adonis_data$R2_2 = paste0(adonis_data$R2_2*100,"%")

adonis_data %>% rownames_to_column("feature") %>% arrange(R2) -> adonis_data

adonis_data$feature = factor(adonis_data$feature , levels = unique(adonis_data$feature))

adonis_data %>%
  mutate(pif = case_when(
    `Pr(>F)`>0.05 ~ "ns",
    `Pr(>F)`<0.001 ~ "***",
    `Pr(>F)`<0.01 ~ "**",
    `Pr(>F)`<= 0.05 ~ "*",
  )) -> adonis_data

write.csv(adonis_data,"merge_adonis_clin_data_R2_result.csv")

# plot adonis------------

adonis_data = read.csv("merge_adonis_clin_data_R2_result.csv",row.names = 1)

nid = c("Stage",
        "Cancer",
        "Group",
        "Smoking_frequency",
        "Smoking_Quantity",
        "Areca",
        "BMI",
        "Age",
        "Alcohol_consumption",
        "Drinking_frequency",
        "Weight",
        "Height",
        "Inheritance")
nid = adonis_data %>% arrange(R2) %>% pull(feature)


adonis_data  = adonis_data %>% filter(feature %in% nid)

adonis_data$feature = factor(adonis_data$feature , levels = rev(nid))

adonis_data %>% 
  mutate(test = case_when(
    `Pr..F.`>0.05 ~ "ns",
    TRUE ~ "signif"
  )) -> adonis_data1

ggplot(adonis_data1 , aes(x = R2, y = feature))+
  geom_col(aes(fill = test))+
  geom_text(aes(label = pif,hjust =0),size = 3,color = "red")+
  scale_fill_manual(values = c("ns" = "#e5e5d8", "signif" = "#D1BCDC")) +
  coord_cartesian(xlim = c(0,0.5))+
  labs(y = "Feature")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y  = element_blank())

ggsave("GCMS_adonis_clin_data_R2_result.pdf", width = 5, height = 3)
