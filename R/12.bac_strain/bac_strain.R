require(tidyverse)
setwd("../15.bac_strain/")
rm(list = ls())

data = read_tsv("fn_strain.txt")
colnames(data)[1] = "Sample_ID"
meta = read.csv("../1.data_tidy/meta_info.csv")

data_plot = data %>% 
  left_join(meta %>% select(Sample_ID,Group),by = "Sample_ID") %>%
  group_by(Sample_ID) %>%
  slice_head(n = 1) %>%
  drop_na()

data_plot %>% 
  group_by(Group) %>% 
  select(starts_with("Fn")) %>%
  summarise(across(.cols = everything(), .fns = ~mean(., na.rm = TRUE))) %>%
  ungroup() -> data_p
data_p  = data_p %>% pivot_longer(-Group,names_to = "type",values_to = "exp")

ggplot(data_p,aes(x = Group, y = exp, fill = type))+
  geom_col(position = "stack",width = 0.5)+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  labs(x = "", y = "Mean expression in \nFusobacterium nucleatum subspecies")
ggsave("fn_sp_diff.pdf",width = 4,height = 5)

data_p2 = data_plot %>% 
  pivot_longer(-c(Sample_ID,Group),names_to = "type",values_to = "exp")

ggplot(data_p2,aes(x = Group, y = exp, fill = type))+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  scale_fill_brewer(palette = "Set1")+
  coord_cartesian(ylim = c(0,100000))+
  facet_wrap(".~type",scales = "free",nrow = 1)+
  theme_bw()+
  labs(x = "", y = "Mean expression in \nFusobacterium nucleatum subspecies")
ggsave("fn_sp_diff2.pdf",width =8,height = 2.3)

data_p2 %>% group_by(type) %>%
  rstatix::wilcox_test(exp ~ Group) %>%
write_csv("fn_sp_diff2.csv")
