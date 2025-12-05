require(tidyverse)

shap_v = read.csv("../2.model/shap_values_0.csv",row.names = 1)
shap_v_abs = abs(shap_v)

tops = data.frame(apply(shap_v_abs,2,mean))

#------------
bacteria <- c("Prevotella_denticola", 
              "Rothia_dentocariosa", 
              "Lancefieldella_rimae", 
              "Aggregatibacter_sp_oral_taxon_458", 
              "Kingella_oralis", 
              "Haemophilus_sputorum", 
              "Actinomyces_dentalis", 
              "Capnocytophaga_granulosa", 
              "Prevotella_baroniae", 
              "Rothia_aeria")

fdata = read.csv("../1.data_tidy/oc_bac_pathway_species.csv",row.names = 1)

sfdata = fdata %>% filter(Species %in% bacteria)

sfdata = sfdata[apply(sfdata[,4:11] > 0,1,sum) > 3,]

sfdata$Pathway = sub(".*: ","",sfdata$Pathway)

sfdata$exp = apply(sfdata[,4:11],1,sum)

require(ggalluvial)
ggplot(sfdata, aes(axis1 = Pathway, axis3 = Species,axis2 = Genus, y = exp)) +
  geom_alluvium(fill = "#d9d9d8",color = "black") +
  geom_stratum(aes(fill = Genus),width = 0.6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(title = "Sankey Diagram using ggalluvial", x = "", y = "Normalized Weight")
ggsave("sankey_plot.pdf",width = 12,height = 8)
