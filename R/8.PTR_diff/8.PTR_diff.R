setwd("./8.PTR_diff/")
rm(list = ls())
require(tidyverse)

# get metadata---------
rm(list = ls())
ptr2 =  read.csv("../1.data_tidy//old/meta_info_fix_0120.csv")
alldata = read.csv("../1.data_tidy/old/meta_info_fix_0120.csv")

intersect(alldata$Sample_ID,ptr2$Sample_ID)

setdiff(alldata$Sample_ID,ptr2$Sample_ID)

ptr3 = ptr2 %>% filter(Sample_ID %in% alldata$Sample_ID)

alldata %>% 
  filter(!Category == "outset") %>%
  pull(Sample_ID) -> nid


ptr2 %>% 
  filter(Sample_ID %in% nid) -> pexp

meta = alldata

meta %>% filter(Sample_ID %in% nid) -> meta

# adonis -------

otu_data2 = pexp %>%  column_to_rownames("Sample_ID")
phe = meta %>% column_to_rownames("Sample_ID")

otu_exp1 = otu_data2
set.seed(666)
source("../addin.r")

meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf,na.rm = T) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray",na.rm =T)

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

ggsave("PTR_adonis_clin_data_R2_result.pdf", width = 5, height = 3)

pexp1 = pexp %>% sjmisc::rotate_df(cn = T, rn = "feature")
pexp = pexp1 %>% 
  column_to_rownames("feature")
pexp[,meta$Sample_ID] -> pexp

pexp


pexp[apply(pexp>0,1,sum) > 6,] -> pexp_f

calculate_fold_change_and_pvalue <- function(data, group) {
  group1 <- data[, group == "HC"]
  group2 <- data[, group == "OSCC"]
  fold_changes <- log2(rowMeans(group2) / rowMeans(group1))
  p_values <- apply(data, 1, function(x) {
    wilcox.test(x[group == "HC"], x[group == "OSCC"])$p.value
  })
  return(list(fold_changes = fold_changes, p_values = p_values))
}
result <- calculate_fold_change_and_pvalue(pexp_f, meta$Group)

fold_changes <- result$fold_changes
p_values <- result$p_values

# OPLAS-DA analysis-----------

X <- t(pexp_f)  

require(ropls)
opls_model <- opls(X, meta$Group, predI = 5, orthoI = 5)

vip_values <- opls_model@vipVn


data <- as.data.frame(opls_model@scoreMN)
data$o1  <- opls_model@orthoScoreMN[,1]
data$group = meta$Group
data$samples = rownames(X)

x_lab <- opls_model@modelDF[1, "R2X"] * 100

colpal=readRDS("../0.colpal/colpal.rds")

# colpal = c(c("OSCC" = "#c9687f", 
#              "HC" = "#92c8ac"),colpal)
# saveRDS(colpal,"../0.colpal/colpal.rds")

ggplot(data,aes(x=p1,y=o1,color=group))+
  theme_bw()+
  geom_point(size=3)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color="red")+
  geom_hline(yintercept = 0,lty="dashed",color="red")+
  annotate("text", x = -Inf, y = -Inf,
           label = "R2Y(cum) = 0.749\nQ2(cum) = 0.464\nRMSEE = 0.255\n", 
           hjust = -0.1, vjust = -0.1)+
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("to1"))+
  stat_ellipse(data=data,
               geom = "polygon",level = 0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = colpal) +
  scale_fill_manual(values = colpal)+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())+
  ggtitle("Scores (OPLS-DA)")
ggsave("oplsda-dimplot.pdf",width = 6,height = 4)

gexp = pexp_f
gexp$log2fc = fold_changes
gexp$pvalue = p_values
gexp$vip = vip_values

gmeta =  read.csv("../0.Rawdata/5.PTR/metab_name.csv")
gmeta$feature = paste0("X",gmeta$feature)

allmeta = gexp %>% 
  rownames_to_column("feature") %>% 
  left_join(gmeta, by = "feature")

write.csv(allmeta, "incohort_all_metab_exp.csv")

gplot = allmeta %>% select(feature,
                           log2fc,
                           vip,
                           pvalue,
                           pvalue,
                           compound)

gplot$Compounds = paste0(gplot$feature,"_",gplot$compound)

gplot_diff = gplot %>% 
  filter(vip >= 1) %>%
  filter(pvalue < 0.05) %>%
  filter(abs(log2fc) > 1)

colnames(gplot_diff)

gplot_diff %>% arrange(log2fc) -> gplot_diff
gplot_diff$Compounds = factor(gplot_diff$Compounds, 
                              levels = gplot_diff$Compounds)
ggplot(gplot_diff, aes(x = log2fc,y = Compounds,fill = log2fc))+
  geom_col(width = 0.8)+
  scale_fill_gradient(low = "#307dc2",high = "#e61e4e")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("metab_diff_col.pdf",width = 6,height = 4)

allmeta %>% filter(feature %in% gplot_diff$feature) -> gplot_data

write.csv(gplot_data, "metab_diff_data.csv")

#meta_name

gmeta =  read.csv("../0.Rawdata/5.PTR/metab_name.csv")
gmeta$feature = paste0("X",gmeta$feature)

allmeta = pexp %>% 
  rownames_to_column("feature") %>% 
  inner_join(gmeta, by = "feature")
write.csv(allmeta, "metab_diff_data_before.csv")


setwd("../8.PTR_diff/")
rm(list = ls())
require(tidyverse)
gmeta = read.csv("PTR_calss.csv")
colnames(gmeta) = c("Metabolite","Class")

pie_data = as.data.frame(table(gmeta$Class)) %>% 
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
  labs(title = "Classification of Metabolites in Human Saliva by PTR-VOCs") +
  scale_fill_manual(values = c("#c9687f", "#92c8ac", "#ff757a", "#f2acbf", "#efc7c7", "#ffa394",
                               "#ffeace", "#cbdeb8", "#92c8ac", "#3e658c", "#89cdd5"))+
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = paste0(round(Freq / sum(Freq) * 100, 1), "%")), 
            position = position_stack(vjust = 0.5))
ggsave("class_I_pie.pdf",width = 8,height = 8)

pexp = read.csv("incohort_all_metab_exp.csv")



pexp[apply(pexp>0,1,sum) > 6,] -> pexp_f

calculate_fold_change_and_pvalue <- function(data, group) {
  group1 <- data[, group == "HC"]
  group2 <- data[, group == "OSCC"]
  fold_changes <- log2(rowMeans(group2) / rowMeans(group1))
  p_values <- apply(data, 1, function(x) {
    wilcox.test(x[group == "HC"], x[group == "OSCC"])$p.value
  })
  return(list(fold_changes = fold_changes, p_values = p_values))
}

result <- calculate_fold_change_and_pvalue(pexp_f, meta$Group)

fold_changes <- result$fold_changes
p_values <- result$p_values