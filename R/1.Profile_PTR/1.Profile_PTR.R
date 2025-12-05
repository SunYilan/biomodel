xbox::chdir("./1.Profile_PTR/")

rm(list = ls())

require(patchwork)
library(vegan)
require(tidyverse)

alldata = read.csv("../1.data_tidy/VOC_PTR_exp.csv")
meta = read.csv("../1.data_tidy/meta_info.csv")


diffset =  alldata %>% 
  column_to_rownames("Sample_ID")


input_metadata = read.csv("../1.data_tidy/meta_info.csv")

input_data = diffset[,c(1:390)]

temp = input_data[,-1]/apply(input_data[,-1],1,sum) 

input_data = cbind(rownames(input_data),temp)

colnames(input_data)[1] = "Sample_ID"

top10 <- names(rev(sort(apply(input_data[,2:390],2,sum))))[1:11]

# 95 is contain

top10  = top10[-8]


top10_data = input_data[,c("Sample_ID",top10)] %>% 
  left_join(meta %>% select(Sample_ID,Group),by = "Sample_ID")

top10_data$Other = 1-apply(top10_data[2:11],1, sum) # 有些加起来大于1？
top10_data$Other[top10_data$Other < 0] = 0
# change names
colnames(top10_data)[ncol(top10_data)] = "Others"

dodge_pdata = top10_data %>% 
  pivot_longer(-c(Sample_ID,Group), values_to = "Abundance", names_to = "Bac")

top10_data %>%  
  arrange(Group) -> temp1

dodge_pdata$Sample_ID = factor(dodge_pdata$Sample_ID ,
                               levels = temp1$Sample_ID)

dodge_pdata$Bac = factor(dodge_pdata$Bac, levels = c(top10,"Others"))

palette =  readRDS("../0.colpal/colpal.rds")

palette = c(
"#89cdd5",
"#c9687f",
"#92c8ac",
"#ff757a",
"#f2acbf",
"#efc7c7",
"#ffa394",
 "#ffeace",
"#cbdeb8",
"#3e658c",
"#4d316c",
"#A3CF9A",
"#987AB8"
)


#saveRDS(c(palette,c("Robert_2024_M" = "#E0B794")),"../../colpal.rds")

ggplot(dodge_pdata,aes(x = Sample_ID,y  = Abundance, fill =  Bac))+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = palette) +
  labs(y = "relative abundance", x = "Time (Days)")+
  ggthemes::theme_clean()+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank(),
        axis.text.x = element_blank())

ggsave("all_cohort_top10_bac_lineplot.pdf",width = 7,height = 4)

input_metadata = read.csv("../1.data_tidy/meta_info.csv")
input_metadata = input_metadata %>% 
  column_to_rownames("Sample_ID") %>% 
  select(Calculus,Group,Exercise,Education_level)

input_metadata = input_metadata[temp1$Sample_ID,] %>% rownames_to_column("Sample_ID")

input_metadata$Sample_ID = factor(input_metadata$Sample_ID ,
                               levels = temp1$Sample_ID)
ggplot(input_metadata,aes(x = Sample_ID, y = 1,fill = as.character(Calculus)))+
  geom_tile() ->p1
ggplot(input_metadata,aes(x = Sample_ID, y = 1,fill = as.character(Group)))+
  geom_tile() ->p2
require(patchwork)
p1/p2
ggsave("all_cohort_top10_bac_lineplot_a1.pdf",width = 12,height = 3)


ggplot(input_metadata,aes(x = Sample_ID, y = 1,fill = as.character(Education_level)))+
  geom_tile()

require(gghalves)
ggplot(data = dodge_pdata %>% 
          filter(!is.na(Group))%>% 
          filter(Bac != "Others")
        ) +
  geom_half_violin(aes(x = Abundance, y = Bac, fill = Group),
               width = 1,coef = 0,side = "l",nudge = 0.01) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif",
                             size = 3,
                             label.y = 0.1) +
  labs(x = "Method", y = "1-Distance")+
  scale_fill_manual(values = c("#ff8bbf","#a0a5e1"))+
  scale_color_manual(values = c("#ff8bbf","#a0a5e1"))+
  facet_wrap(".~method",ncol = 1,scales = "free")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("beta_distance_4method_bk.pdf",width = 5,height = 8)



# sankey plot------------

dodge_pdata_pc1c2 <- dodge_pdata %>%
  group_by(Group,Bac) %>%
  summarise(Abundance = mean(Abundance))

library(ggalluvial)

dodge_pdata_pc1c2 %>%
  arrange(Bac,Group) %>% 
  group_by(Bac) %>% 
  mutate(Abundance_diff = Abundance - lag(Abundance)) %>%
  drop_na() %>% 
  ungroup() %>% 
  mutate(updown = case_when(
    Abundance_diff > 0 ~ "up",
    Abundance_diff <= 0 ~ "down"
  ),
  Abundance_cumsum = cumsum(Abundance),
  Abundance = 1- Abundance_cumsum)-> diff_text



ggplot(dodge_pdata_pc1c2, aes(x = Group,
                              y  = Abundance, 
                              fill =  Bac,
                              stratum = Bac, alluvium = Bac))+
  geom_bar(stat = "identity", position = "stack",width = 0.5) +
  geom_text(data = diff_text,aes(x = 2.5,label = Bac, color = updown))+
  geom_stratum(width = 0.5, color='white')+
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "sine")+
  scale_fill_manual(values = palette) +
  scale_color_manual(values = c("#A3CF9A","#987AB8")) +
  labs(y = "relative abundance", x = "")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("source_bac_sankey.pdf",width = 5,height = 4)

# a diversity------------
source("../addin.r")

data_exp = input_data 

count_data = round(data_exp*10000)
alpha <- alpha_diversity(count_data)

alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp1

a_data = plot_a_data_temp1 %>% left_join(input_metadata, by = "Sample_ID")

write.csv(a_data,"a_diversity_value.csv")

ggplot(a_data %>% select(Shannon,Simpson,Group) %>% hablar::retype(),
       aes(x = Group, y = Simpson,fill = Group))+
  geom_boxplot()+
  ggpubr::stat_compare_means()+
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(y = "a Diversity", x = "Time (Days)")+
  ggthemes::theme_few()+
  ggtitle("Simpson")-> p1
ggplot(a_data %>% select(Shannon,Simpson,Group) %>% hablar::retype(),
       aes(x = Group, y = Shannon,fill = Group))+
  geom_boxplot()+
  ggpubr::stat_compare_means()+
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(y = "a Diversity", x = "Time (Days)")+
  ggthemes::theme_few()+
  ggtitle("Shannon")-> p2

p1+p2
ggsave("a_diversity_boxplot.pdf",width = 7,height = 3)

# b diversity-----------
meio.mdf <- as.matrix.data.frame(input_data)

meio.mdf[meio.mdf < 0.00001] = 0.00001

require(vegan)

meio.bray <- vegdist(meio.mdf, method = "bray")
pcoa.meio.bray <- cmdscale(meio.bray, k = 2, eig = T)

pcoa.meio.bray.plotting <- as.data.frame(pcoa.meio.bray$points)
colnames(pcoa.meio.bray.plotting) <- c("axis_1", "axis_2")

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pc1_r = pcoa.meio.bray$eig[1]/(sum(pcoa.meio.bray$eig))
pc2_r = pcoa.meio.bray$eig[2]/(sum(pcoa.meio.bray$eig))


merge_sum = input_data %>% 
  rownames_to_column("Sample_ID") %>% 
  left_join(input_metadata, by = "Sample_ID")
# Age_DOL
pcoa.meio.bray.plotting$Sample_ID <- merge_sum$Sample_ID
pcoa.meio.bray.plotting$Group <- merge_sum$Group

write.csv(pcoa.meio.bray.plotting,"pcao_meio,bray.data.csv")

# anova p值检验
source("../addin.r")
anova_data = data.frame()

for(i in c("Group")){
  anova_result = adonis3(as.formula(paste("meio.bray","~",i)),
                         data = pcoa.meio.bray.plotting,
                         permutations = 999,
                         na.action=na.omit)
  anova_temp = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
  anova_data = rbind(anova_data,anova_temp)
}
write.csv(anova_data,"bdiversity_adonis.csv")

colpal = palette

# 直接分组-----------------
pcoa.meio.bray.plotting %>% group_by(Group) %>%
  summarise(pc1 = mean(axis_1),
            pc2 = mean(axis_2)) -> g_p

ggplot(pcoa.meio.bray.plotting,
       aes(x = axis_1, y = axis_2, colour = Group)) +
  geom_point(size = 3) +
  ggplot2::stat_ellipse(geom = "polygon", 
                        mapping = ggplot2::aes(color = Group),fill = "white", type = "t", 
                        level = 0.95, linetype = 2, alpha = 0)+
  #geom_point(data = g_p,aes(x =  pc1, y = pc2),size = 2)+
  ggthemes::theme_few() +
  scale_color_manual(values = colpal)+
  xlab(paste0("PCoA 1 (",round(pc1_r,3)*100,")%")) +
  ylab(paste0("PCoA 2 (",round(pc2_r,3)*100,")%")) +
  annotate(geom = 'text', label = paste0('ANOVA P:',"0.001", "\n R2:",0.032*100,"%"), x = -Inf, y = -Inf, hjust = -0.1,vjust = -1)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("pcoa.pdf",width = 6,height = 4)
