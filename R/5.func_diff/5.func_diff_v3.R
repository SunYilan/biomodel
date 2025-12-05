setwd("../5.func_diff/")
rm(list  = ls())
#-------
library(tidyverse)

alldata = read.csv("../1.data_tidy/OSCC_P.sum.1212.csv",row.names = 1)

diffset =  alldata %>% 
  filter(!Category == "outset") %>%
  column_to_rownames("Sample_ID")


input_metadata = diffset[,c(19232:19245)]
input_data = diffset[,c(2:19231)]

library(Maaslin2)
set.seed(20241213)

fit_data_sepsis = Maaslin2(
  input_data, input_metadata,
  min_prevalence = 0.001,
  normalization = "NONE",
  transform = "NONE",
  output = 'maaslin2_diff_mix_OSCC_vs_HC', 
  fixed_effects = c("Group","Sex", 
                    "Age", "BMI"),
  reference = c('Group,HC')
)


ms2 = read_tsv("maaslin2_diff_mix_OSCC_vs_HC//all_results.xls") %>%
  filter(metadata == "Group")
ms2$qval = p.adjust(ms2$pval, method = "BH")

data = ms2
data$neg_log10_qval <- -log10(data$qval)


ms2_q = data %>% 
  dplyr::filter(qval < 0.05) %>%
  mutate(updown = case_when(
    coef > 0 ~ "OSCC",
    coef < 0 ~ "HC",
    coef == 0 ~ "None"
  ))

ms2_q %>% arrange(coef) -> ms2_q


ms2_q$feature = factor(ms2_q$feature,levels = ms2_q$feature)

ggplot(ms2_q, aes(x = coef, y = feature,fill = updown)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("OSCC" = "#c9687f", 
                               "HC" = "#92c8ac")) +
  coord_cartesian(xlim = c(-40,45))+
  ggthemes::theme_clean(base_size = 15)+
  theme(panel.grid.major.y  = element_blank())+
  labs(x = "coef", y = "", fill = "Group")
ggsave("maaslin2_diff.pdf",width = 16,height = 24)


f_diff = input_data[,colnames(input_data) %in% ms2_q$feature]

write.csv(f_diff,"func_diff_1229.csv")
