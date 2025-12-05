dir.create("14.biaoqu")
setwd("14.biaoqu/")


require(tidyverse)
require(ggpmisc)
data = read_csv("./area_peak.csv")

ggplot(data,aes(y = ppb, x = `49`))+
  geom_point(color = "#f72585")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  labs(x = "Peak Area (a.u.)", y = "Analyte Concentration (ppb)")+
  theme_bw()
ggsave("biaoqu.pdf",width =4 ,height = 3.5)


data = read.csv("../1.data_tidy/VOC_PTR_exp.csv") %>% select(Sample_ID,X49)
meta = read.csv("../1.data_tidy/meta_info.csv")

data_plot = data %>% left_join(meta %>% select(Sample_ID,Group),by = "Sample_ID") %>% drop_na()

data_plot = data_plot %>% mutate(abs_f = X49*0.0164+1.45)
colpal = readRDS("../0.colpal/colpal.rds")
ggplot(data_plot,aes(x = Group, y = abs_f,fill = Group))+
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means()+
  coord_cartesian(ylim = c(0,100))+
  scale_fill_manual(values = colpal)+
  labs(x = "", y = "Methanethiol absolute concentration")+
  theme_bw()

ggsave("meth_boxplpt.pdf",width =4 ,height = 3.5)






