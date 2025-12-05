setwd("../4.masslin2_bacdiff")
rm(list  = ls())

#-------
library(tidyverse)

alldata = read.csv("../1.data_tidy/metagenome_exp.csv",row.names = 1)

diffset =  alldata %>% 
  filter(!Category == "outset")


input_metadata = diffset[,c(426:439)]
input_data = diffset[,c(2:425)]

library(Maaslin2)
set.seed(20241213)

fit_data_sepsis = Maaslin2(
  input_data, input_metadata,
  min_prevalence = 0.001,
  normalization = "NONE",
  transform = "NONE",
  output = 'maaslin2_diff_mix_OSCC_vs_HC', 
  fixed_effects = c("Group"),
  reference = c('Group,HC')
)

shapdata = read.csv("../2.model/shapley_data.csv")
colnames(shapdata) = c("feature","shape_value")

ms2 = read_tsv("maaslin2_diff_mix_OSCC_vs_HC//all_results.tsv") %>%
  filter(metadata == "Group")

ms2$qval = p.adjust(ms2$pval, method = "BH")

ms2 = ms2 %>% left_join(shapdata, by= c("feature"))

data = ms2
data$neg_log10_pval <- -log10(data$pval)


ms2_q = data %>% 
  filter(pval < 0.05) %>%
  mutate(updown = case_when(
    coef > 0 ~ "OSCC",
    coef < 0 ~ "HC",
    coef == 0 ~ "None"
  ))

ggplot(data, aes(x = coef, y = neg_log10_pval,size = shape_value)) +
  geom_point(color = "#d9d9d8", size = 3) +
  geom_point(data = ms2_q ,aes(color = updown)) +
  scale_color_manual(values = c("OSCC" = "#c9687f", 
                                "HC" = "#92c8ac",
                                "None" = "#d9d9d8")) +
  ggrepel::geom_text_repel(aes(label = ifelse(pval < 0.05 & abs(coef) > 0.002, as.character(feature), "")), 
                           hjust = 0, vjust = 0, nudge_x = 0.0005) +
  geom_vline(xintercept = c(0.002,-0.002))+
  coord_cartesian(xlim = c(-0.04,0.04))+
  scale_x_continuous(breaks = c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03),
                     labels = c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03))+
  ggthemes::theme_clean(base_size = 30)+
  theme(panel.grid.major.y  = element_blank())+
  labs(x = "coef", y = "Negtive log10(Q value)", color = "FDR < 0.05")
ggsave("maaslin2_diff.pdf",width = 12,height = 8)
ggsave("maaslin2_diff_big.pdf",width = 24,height = 16)

ggplot(data %>% 
         filter(pval >= 0.05), aes(x = coef, y = neg_log10_pval,size = shape_value)) +
  geom_point(color = "#d9d9d8", size = 3) +
  geom_point(data = ms2_q ,aes(color = updown)) +
  scale_color_manual(values = c("OSCC" = "#c9687f", 
                                "HC" = "#92c8ac",
                                "None" = "#d9d9d8")) +
  geom_vline(xintercept = c(0.002,-0.002))+
  coord_cartesian(xlim = c(-0.02,0.02))+
  theme_bw()+
  theme(panel.grid.major.y  = element_blank())+
  labs(x = "coef", y = "Negtive log10(Q value)", color = "FDR < 0.05")

ggsave("maaslin2_diff_sim.pdf",width = 5,height = 4)

shape_top50 = ms2 %>% slice_max(n = 50, order_by = shape_value)

fdrdata = ms2 %>% filter(pval < 0.05) %>%  filter(abs(coef) > 0.002)

pdf("massline_shape_venn.pdf",width = 4,height = 4)
re = xbox::vennplot(shape_top50$feature, fdrdata$feature)
dev.off()

ms2 %>% 
  filter(feature %in% re[["id"]][["xy"]]) %>%
  write.csv("massline_shape_venn.csv")


setwd("../4.masslin2_bacdiff/")

ms2 = read.csv("massline_shape_venn.csv",row.names = 1)
ms2 %>% mutate(updown = case_when(
  coef > 0 ~ "OSCC",
  coef < 0 ~ "HC",
  coef == 0 ~ "None"
)) -> ms2
ms2 %>% arrange(coef) -> ms2
ms2$feature = factor(ms2$feature,levels = ms2$feature)
ggplot(data = ms2,
       aes(y = feature, x = coef))+
  geom_linerange(aes(color = updown,xmin = 0,xmax = coef),linewidth = 1)+
  geom_point(aes(color = updown,size = shape_value),alpha = 1)+
  scale_color_manual(values = c("OSCC" = "#c9687f", 
                                "HC" = "#92c8ac")) +
  ggthemes::theme_clean(base_size = 15)+
  theme(panel.grid.major.y  = element_blank())+
  labs(x = "coef", y = "", fill = "Group")
ggsave("intersect52_diff.pdf",width = 8,height = 6)

nid = ms2  %>% pull(feature)

input_data[,colnames(input_data) %in% nid] -> outdata

write.csv(outdata,"bac_diff_1229.csv")
