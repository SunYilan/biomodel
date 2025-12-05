setwd("../2.model/")
require(tidyverse)
rm(list = ls())
# split data-----------------
alldata = read.csv("../1.data_tidy/metagenome_exp.csv",row.names = 1)

outset =  alldata %>% filter(Category == "outset")
set.seed("20241212")
alldata %>% 
  filter(Category == "onset") %>%
  sample_frac(0.7) -> trainset
alldata %>%
  filter(Category == "onset") %>%
  filter(!Sample_ID %in% (trainset %>% pull(Sample_ID))) -> vaildset

write.csv(outset,"outset.csv")
write.csv(trainset,"trainset.csv")
write.csv(vaildset,"vaildset.csv")

# model select ---------------
model_re = read.csv("model_compare_results.csv",row.names = 1)

model_re = model_re %>% 
  arrange(AUC)
model_re$Model = factor(model_re$Model,levels = model_re$Model)

ggplot(model_re %>% slice_max(n = 10, order_by = AUC),aes(y = Model,x = AUC)) +
  geom_col(aes(fill = AUC),width = 0.6) +
  scale_fill_gradient(low = "#307dc2",high = "#e61e4e")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  coord_cartesian(xlim = c(0.6,1))+
  labs(y = "Trained Model", x = "AUC Score")

ggsave("model_select.pdf",height = 3,width = 6)


# shape value ---------------

shapley_results = read.csv("shap_values_0.csv",row.names = 1)

shapley_results_abs_sum = apply(abs(shapley_results),2,mean)
data.frame(shapley_results_abs_sum) %>%
  rename("shapley_value" = "shapley_results_abs_sum")  %>%
  rownames_to_column("bac") %>%
  slice_max(order_by = shapley_value, n = 10) -> shape_plot

shape_plot %>%
  add_row(shapley_value = 0.03,bac = "aa") -> shape_plot1

shape_plot1$bac = factor(shape_plot1$bac,
                         levels = c(rev(c("x",as.character(shape_plot$bac))),"aa")
                         )
colpal = readRDS("../0.colpal/colpal.rds")
ggplot(shape_plot1,
       aes(y = bac, x = shapley_value))+
  geom_col(aes(fill = bac))+
  scale_fill_manual(values = colpal)+
  coord_polar()+
  scale_x_continuous(breaks = seq(0,0.03,0.01),
                     labels = seq(0,0.03,0.01))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("shape_plot.data_s1.pdf",width = 10,height = 3)

write.csv(data.frame(shapley_results_abs_sum),"shapley_data.csv")

#build model--------------
library(caret)
library(randomForest)
trainset = trainset[,c(1:424,439)]
vaildset = vaildset[,c(1:424,439)]
outset = outset[,c(1:424,439)]

trainset$Group <- as.factor(trainset$Group)
vaildset$Group <- as.factor(vaildset$Group)
outset$Group <- as.factor(outset$Group)

print(run_model)
best_model <- run_model$finalModel
write.csv(data.frame(run_model$results),"tarin_model_result.csv")

saveRDS(run_model,"run_RF_model.rds")

library(plotROC)
library(pROC)

auc_values <- run_model$pred %>%
  group_by(Resample) %>%  # 按折分组
  summarise(auc = as.numeric(roc(obs, OSCC)$auc)) %>% 
  mutate(auc_label = paste0("Fold ", 
                            gsub("Resample", "", Resample), ": AUC = ", 
                            round(auc, 3)))
    
    
    
run_model$pred$obs <- factor(run_model$pred$obs, levels = c("HC", "OSCC"))

ggplot() +
  geom_roc(data = run_model$pred, 
           aes(d = obs, m = OSCC, color = Resample),
           n.cuts = 0,show.legend = FALSE) +  # 不显示图例（用文字代替）
  style_roc() +
  geom_text(
    data = auc_values,
    aes(x = 0.6, 
        y = 0.2 - 0.05 * as.numeric(factor(Resample)),
        color = Resample),  # 调整标签位置
        label = auc_values$auc_label,
        hjust = 0,
        size = 3
    ) +
      ggtitle("10-Fold CV ROC Curves with AUC") +
      scale_color_brewer(palette = "Set1")+
  theme(legend.position = "none") -> p3

predictions <- predict(run_model, newdata = vaildset[,-425],type = "prob")

library(pROC)
library(multipleROC)

roc_curve <- roc(vaildset$Group, predictions[, "HC"])
auc_value <- round(auc(roc_curve),3)

ggroc(roc_curve,
      legacy.axes = T,
      color="#ffa394",
      linewidth=1)+
  annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1,
           linetype="dashed",color="#186F65")+
  annotate("text",x=0.8,y=0.3,label=paste0("AUC = ",auc_value))+
  labs(x="1-Specificity",y="Sensitivity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Model Vaild Set(30%)") -> p1

predictions <- predict(run_model, newdata = outset[,-425],type = "prob")

roc_curve <- roc(outset$Group, predictions[, "HC"])
auc_value <- round(auc(roc_curve),3)

ggroc(roc_curve,
      legacy.axes = T,
      color="#cbdeb8",
      linewidth=1)+
  annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1,
           linetype="dashed",color="#186F65")+
  annotate("text",x=0.8,y=0.3,label=paste0("AUC = ",auc_value))+
  labs(x="1-Specificity",y="Sensitivity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("External Cohort") -> p2
require(patchwork)
p3+p2+p1
ggsave("Model_test_ROC.pdf",width = 12,height = 4)


data17 = read.csv("../20.yanzheng_cohort/OSCC_M.sum.1212.csv",row.names = 1)
data17_p = data17 %>% 
  column_to_rownames("Sample") %>% 
  select(-c(Sample_ID,Group))

nid = setdiff(colnames(outset[,-425]),colnames(data17_p))
for (feature in nid) {
  data17_p[[feature]] <- 0
}


predictions <- predict(run_model, newdata = data17_p,type = "prob")

roc_curve <- roc(data17$Group, predictions[, "HC"])
auc_value <- round(auc(roc_curve),3)

ggroc(roc_curve,
      legacy.axes = T,
      color="#cbdeb8",
      linewidth=1)+
  annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1,
           linetype="dashed",color="#186F65")+
  annotate("text",x=0.8,y=0.3,label=paste0("AUC = ",auc_value))+
  labs(x="1-Specificity",y="Sensitivity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("External Cohort") -> p2

ggsave("data17_roc.pdf",width = 3,height = 2.8)

#----------------


predictions1 <- predict(run_model, newdata = outset[,-425])
predictions2 <- predict(run_model, newdata = vaildset[,-425])
predictions3 = run_model[["finalModel"]]$predicted


predd1 = data.frame(rowid = rownames(outset),obs=outset$Group,pred = predictions1)
predd2 = data.frame(rowid = rownames(vaildset),obs=vaildset$Group,pred = predictions2)
predd3 = data.frame(rowid = rownames(trainset),obs=trainset$Group,pred = predictions3)

predd = rbind(predd1,predd2,predd3)

write.csv(predd,"best_predict_all_sample.csv")





#分stage来做----------
setwd("12.model_all/")
require(tidyverse)

meta = read.csv("../1.data_tidy/meta_info.csv")
pred = read.csv("best_predict_all_sample.csv")

meta %>% select(Sample_ID,Stage,Cancer) %>% 
  replace_na(list(Stage = "Control")) %>% 
  replace_na(list(Cancer = "Control")) -> meta_s

pred %>% left_join(meta_s,by = c("rowid" = "Sample_ID")) %>% 
  replace_na(list(Stage = "Control")) %>% 
  replace_na(list(Cancer = "Control")) -> pred

pred %>% mutate(ident = case_when(
  obs == pred ~ 1,
  obs != pred ~ 0
)) %>% group_by(Stage) %>%
  summarise(a = sum(ident),
            b = n()) %>%
  mutate(ratio = a/b*100) %>% 
  mutate(label = paste0(a,"/",b," (",round(ratio,2),"%)")) %>%  
  arrange(ratio) -> pred_plot_stage

pred_plot_stage$Stage = factor(pred_plot_stage$Stage,levels = pred_plot_stage$Stage)
ggplot(pred_plot_stage,aes(y = Stage,x = ratio))+
  geom_linerange(aes(xmin = 0,xmax=  ratio),
                 size = 1,alpha = 0.8,color = "#f99417") +
  geom_point(color = "#c9687f",size = 4)+
  geom_text(aes(x = ratio/2,label = label),vjust =-0.5)+
  labs(y = "Cancer Stage", x = "Prediction Accuracy (%)")+
  ggthemes::theme_clean() -> p4

pred %>% mutate(ident = case_when(
  obs == pred ~ 1,
  obs != pred ~ 0
)) %>% group_by(Cancer) %>%
  summarise(a = sum(ident),
            b = n()) %>%
  mutate(ratio = a/b*100) %>% 
  mutate(label = paste0(a,"/",b," (",round(ratio,2),"%)")) %>%  
  arrange(ratio) -> pred_plot_stage

pred_plot_stage$Cancer = factor(pred_plot_stage$Cancer,levels = pred_plot_stage$Cancer)
ggplot(pred_plot_stage,aes(y = Cancer,x = ratio))+
  geom_linerange(aes(xmin = 0,xmax=  ratio),
                 size = 1,alpha = 0.8,color = "#f99417") +
  geom_point(color = "#c9687f",size = 4)+
  geom_text(aes(x = ratio/2,label = label),vjust =-0.5)+
  labs(y = "Cancer Type", x = "Prediction Accuracy (%)")+
  ggthemes::theme_clean() -> p5

p4/p5+plot_layout(heights = c(3.5,1))

ggsave("Prediction Accuracy_stage.pdf",width =6,height = 5)