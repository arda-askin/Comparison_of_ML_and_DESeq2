library(h2o)
h2o.init(nthreads = -1) 

dl_n <- h2o.loadModel("C:\\TCGA_SKCM_model\\all_models\\ann_valid_fi.RDS\\dl_grid_no_model_4")
dl_s <- h2o.loadModel("C:\\TCGA_SKCM_model\\all_models\\noise_ifn_ann_5.rds\\dl_grid_5_model_4")
dl_m <- h2o.loadModel("C:\\TCGA_SKCM_model\\all_models\\noise_ifn_ann_0.5.rds\\dl_grid_0.5_model_4")
dl_l <- h2o.loadModel("C:\\TCGA_SKCM_model\\all_models\\noise_ifn_ann_0.05.rds\\dl_grid_0.05_model_4")
dl_xl <- h2o.loadModel("C:\\TCGA_SKCM_model\\all_models\\noise_ifn_ann_0.005.rds\\dl_grid_0.005_model_1")




Importance_func <- function(
    no,
    small ,
    medium,
    large,
    xlarge,
){
  
  # browser()
  
  models <- list(
    no = no,
    small = small,
    medium = medium,
    large = large,
    extra_large = xlarge
  )
  
  
  importance_results <- data.frame(
    Model = character(),
    Important_Gene_Number = numeric(),
    Ifng_Gene_Number = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  
  for (model_index in 1 : length(names(models))) {
    
    model <- models[[model_index]]
    
    
    # Check if the model is an H2O model
    if (class(model) == "H2OBinomialModel") {
      
      
      
      modelimp <- h2o.varimp(model)
      modelimp_d <- as.data.frame(modelimp)
      modelimp_d <- filter(modelimp_d, !scaled_mportance == 0)
      modelimp_d$gene <- rownames(modelimp_d) 
      
      gene_importance_score <- modelimp_d$scaled_mportance
      gene_name <- modelimp_d$variable
      condition <- model_index
      
    } else {
      
      modelImp <- varImp(model, scale = T)
      modelImp_d <- data.frame(modelImp$importance)
      modelImp_d <- modelImp_d %>% filter(!Overall == 0)
      modelImp_d$gene <- rownames(modelImp_d)
      
      
      important_gene_number <- nrow(modelImp_d)
      ifng_gene_number <- length(intersect(modelImp_d$gene, gene_set))
      
    }
    
    # Store the results in the data frame
    importance_results <- rbind(importance_results, data.frame(
      gene_importance_score = gene_importance_score,
      gene_name = gene_name,
      condition = paste0("condition", condition)
    ))
    
  }
  
  importance_results
  
  
  }




dl_imp <- Importance_func(dl_n,
                          dl_s,
                          dl_m,
                          dl_l,
                          dl_xl,
                          infg_gene_set)




library(tidyr)

# Pivot to wide format
wide_varimp <- dl_imp %>%
  pivot_wider(names_from = condition, values_from = gene_importance_score)


# Example data frame: Variables x Conditions
importance_scores <- data.frame(
  Variable = paste0("Variable", 1:100),
  Condition1 = runif(100, 0, 1),
  Condition2 = runif(100, 0, 1),
  Condition3 = runif(100, 0, 1)
)


importance_scores <- wide_varimp

importance_scores$Diff_Cond1_Cond2 <- importance_scores$condition2 - importance_scores$condition1
importance_scores$Diff_Cond1_Cond3 <- importance_scores$condition3 - importance_scores$condition1
importance_scores$Diff_Cond1_Cond4 <- importance_scores$condition4 - importance_scores$condition1
importance_scores$Diff_Cond1_Cond5 <- importance_scores$condition5 - importance_scores$condition1

importance_scores$Pct_Change_Cond1_Cond2 <- (importance_scores$Diff_Cond1_Cond2 / importance_scores$condition1)*100 
importance_scores$Pct_Change_Cond1_Cond3 <- (importance_scores$Diff_Cond1_Cond3 / importance_scores$condition1) *100 
importance_scores$Pct_Change_Cond1_Cond4 <- (importance_scores$Diff_Cond1_Cond4 / importance_scores$condition1) *100 
importance_scores$Pct_Change_Cond1_Cond5 <- (importance_scores$Diff_Cond1_Cond5 / importance_scores$condition1) *100 



importance_scores <- filter(importance_scores, importance_scores$gene_name %in% intersect(importance_scores$gene_name, infg_gene_set))

library(RColorBrewer)

# Set a color palette with 12 distinguishable colors
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3", 
            "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", 
            "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
            "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3")


min_value = -100
max_value = +200

library(ggplot2)

# Melt the data for ggplot
library(reshape2)

importance_scores <- importance_scores[,-c(2:10)]

melted_data <- melt(importance_scores, id.vars = "gene_name")

# Plot importance scores for each condition
ggplot(melted_data, aes(x = gene_name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank()  
  )+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values = colors)

sub1 <- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond2" )

ggplot(sub1, aes(x = gene_name, y = abs(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank()  
  )
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values = "#7570B3")+
  ylim(c(min_value,max_value))

sub2 <- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond3" )

ggplot(sub2, aes(x = gene_name, y = abs(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_blank())+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values = "#A65628")+
  ylim(c(min_value,max_value))


sub3<- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond4" )

ggplot(sub3, aes(x = gene_name, y = abs(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
theme(axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values =  "#1B9E77")+
  ylim(c(min_value,max_value))


sub4<- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond5" )

ggplot(sub4, aes(x = gene_name, y = abs(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme( axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values =  "#377EB8")+
  ylim(c(min_value,max_value))



sub5 <- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond2", value >0 )

ggplot(sub5, aes(x = gene_name, y = abs(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
theme( axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values = "#7570B3")

sub6 <- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond3", value >0 )

ggplot(sub6, aes(x = gene_name, y = abs(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
theme( axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values = "#A65628")


sub7 <- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond2", value <0 )

ggplot(sub7, aes(x = gene_name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme( axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values = "#7570B3")


sub8 <- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond3", value <0 )

ggplot(sub8, aes(x = gene_name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
theme(axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values = "#A65628")

sub9<- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond4", value >0)

ggplot(sub9, aes(x = gene_name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme(axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values =  "#1B9E77")+
  ylim(c(min_value,max_value))

sub10<- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond4", value < 0)

ggplot(sub10, aes(x = gene_name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme(axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values =  "#1B9E77")

sub11<- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond5", value >0)

ggplot(sub11, aes(x = gene_name, y = abs(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme(axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values =  "#377EB8")+
  ylim(c(min_value,max_value))

sub11<- melted_data %>% filter(variable == "Pct_Change_Cond1_Cond5", value <0)

ggplot(sub11, aes(x = gene_name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +  
  theme(axis.text.x = element_blank(),  
axis.ticks.x = element_blank()  
)+
  labs(y = "Importance Score", fill = "Condition") +
  scale_fill_manual(values =  "#377EB8")


  