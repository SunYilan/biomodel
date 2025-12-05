#dir.create("10.factor")
setwd("10.factor")


require(vegan)
require(tidyverse)

rm(list = ls())


alldata = read.csv("../1.data_tidy/metagenome_exp.csv")
meta = read.csv("../1.data_tidy/meta_info.csv")

diffset =  alldata %>% 
  filter(!Category == "outset") %>%
  column_to_rownames("Sample_ID")

diifmeta = meta %>%  
  filter(!Category == "outset") %>%
  column_to_rownames("Sample_ID")

otu_data2 = diffset[,c(2:425)]

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

# nid = c("Stage",
#         "Cancer",
#         "Group",
#         "Smoking_frequency",
#         "Smoking_Quantity",
#         "Areca",
#         "BMI",
#         "Age",
#         "Alcohol_consumption",
#         "Drinking_frequency",
#         "Weight",
#         "Height",
#         "Inheritance")

nid = adonis_data %>% arrange(R2) %>% pull(feature)

adonis_data  = adonis_data %>% filter(feature %in% nid)

adonis_data$feature = factor(adonis_data$feature , levels = nid)

adonis_data %>% 
  mutate(test = case_when(
    `Pr..F.`>0.05 ~ "ns",
    TRUE ~ "signif"
  )) -> adonis_data1

ggplot(adonis_data1 , aes(x = R2, y = feature))+
  geom_col(aes(fill = test))+
  geom_text(aes(label = pif,hjust =0),size = 3,color = "red")+
  scale_fill_manual(values = c("ns" = "#e5e5d8", "signif" = "#D1BCDC")) +
  coord_cartesian(xlim = c(0,0.2))+
  labs(y = "Feature")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y  = element_blank())

ggsave("metagenome_adonis_clin_data_R2_result.pdf", width = 5, height = 8)

# OR value---------
rm(list = ls())
setwd("10.factor/")
meta = read.csv("../1.data_tidy/meta_info.csv")
require(tidyverse)

meta$Sex = ifelse(meta$Sex == "Female",0,1)
meta$Inheritance = ifelse(meta$Inheritance == "N",0,1)
meta$Areca = ifelse(meta$Areca == "N",0,1)
table(meta$Stage)
meta %>% mutate(Cancer = case_when(
  Cancer == "Post-Chemotherapy/Post-Treatment/Post-Neoadjuvant" ~ 2,
  Cancer == "Primary" ~ 1,
  Cancer == "Untreated Relapse" ~ 3,
  is.na(Cancer) ~ 0
)) %>% mutate(Stage = case_when(
  is.na(Stage) ~ 0,
  Stage == "I" ~ 1,
  Stage == "rI" ~ 1.5,
  Stage == "II" ~ 2,
  Stage == "rII" ~ 2.5,
  Stage == "III" ~ 3,
  Stage == "rIII" ~ 3.5,
  Stage == "IVA" ~ 4,
  Stage == "rIVA" ~ 4.5,
  Stage == "IVB" ~ 5,
  Stage == "rIVB" ~ 5.5,
  Stage == "IVC" ~ 6,
)) -> meta
meta$Group <- ifelse(meta$Group == "HC", 0, 1)


# get 
diffmeta = meta %>%  
  filter(!Category == "outset") %>%
  column_to_rownames("Sample_ID")

df = diffmeta[,-c(1,2)]
df <- hablar::retype(df)

# Normalization data------------

cols_to_scale <- colnames(df)[-c(12)]
convert_to_scale <- function(x){
  z = (x - min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))
}
df <- df %>%
  mutate(across(all_of(cols_to_scale), convert_to_scale))

require(ggcor)
quickcor(df[,-c(13,14)],cor.test = TRUE, use = "complete.obs") +
  geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
  geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 2.5) +
  geom_abline(slope = -1, intercept = 34)+
  scale_fill_gradient2(low = "#41ab5d",mid = "#f1fb67",high = "#e4007c") -> cor_data

ggsave(plot = cor_data,"mutivar_cor_spearman.pdf",width = 8,height = 6)

cdata = cor_data[["data"]]

cdata$p.adj = p.adjust(cdata$p.value,method = "BH")

cdata %>% 
  filter(`.col.names` == "Group") %>% 
  #filter(p.value < 0.05) -> lk
  mutate(sig = case_when(
    p.adj > 0.05 ~ "No",
    p.adj <= 0.05 ~ "FDR < 0.05"
  )) %>% 
  mutate(updown = case_when(
    r > 0 ~ "Risk",
    r < 0 ~ "Protective"
  )) %>% filter(`.row.names` != "Group")-> cdata