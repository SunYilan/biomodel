require(tidyverse)

data = read.csv("data.levels_exp.csv",row.names = 1)

data[is.na(data)] = 0

data_diff = data %>% 
  rownames_to_column("feature") %>% 
  pivot_longer(-feature,names_to = "sample",values_to = "exp") %>% 
  mutate(group = case_when(
    grepl("OSCC",sample) ~ "OSCC",
    grepl("PO",sample) ~ "OSCC",
    grepl("HC",sample) ~ "HC",
    grepl("WK",sample) ~ "HC"
  ))

require(rstatix)

data_diff %>% group_by(feature) %>% 
  rstatix::t_test(exp~group) ->data_test

ggplot(data_diff,aes(x = group,y = exp,fill = group))+
  geom_boxplot()+
  ggpubr::stat_compare_means(vjust = 1)+
  facet_wrap(".~feature",scales = "free_y")+
  scale_fill_manual(values = c("#A3CF9A","#987AB8")) +
  ggthemes::theme_few()
ggsave("data.levels_exp.pdf",width = 14,height = 10)
