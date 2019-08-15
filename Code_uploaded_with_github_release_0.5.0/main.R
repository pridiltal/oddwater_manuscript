## ---- load
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(plotly)
library(lubridate)
library(GGally)
library(zoo)
library(oddwater)
library(DDoutlier)
library(kableExtra)
library(patchwork)
library(microbenchmark)

## Load user-defined functions
source("src/functions.R")



###########################################################
#### Analysis of data from Sandy Creek
###########################################################

## ---- filter_outliers_from_anormalies

data("data_sandy_anom")
## Change the format and the class of the existing variable
data_sandy_anom$Timestamp <- lubridate::dmy_hm(data_sandy_anom$Timestamp)
# Turbidity
data_sandy_anom$out_type_Tur <- ifelse(data_sandy_anom$type_Tur %in% c("B", "C", "E", "H", "L"), "0", as.character(data_sandy_anom$type_Tur))
# Conductivity
data_sandy_anom$out_type_Cond <- ifelse(data_sandy_anom$type_Cond %in% c("B", "C", "E", "H", "L"), "0", as.character(data_sandy_anom$type_Cond))
# Level
data_sandy_anom$out_type_Level <- ifelse(data_sandy_anom$type_Level %in% c("B", "C", "E", "H", "L"), "0", as.character(data_sandy_anom$type_Level))
# creat a new dataset including only outliers
data_sandy_out <- data_sandy_anom[, c("Timestamp", "Level", "Cond", "Tur", "out_type_Tur", "out_type_Cond", "out_type_Level")]
data_sandy_out$out_label_Tur <- ifelse((data_sandy_out$out_type_Tur == "0"), 0, 1)
data_sandy_out$out_label_Cond <- ifelse((data_sandy_out$out_type_Cond == "0"), 0, 1)
data_sandy_out$out_label_Level <- ifelse((data_sandy_out$out_type_Level == "0"), 0, 1)

data_sandy_out[, c(8:10)] <- lapply(
  data_sandy_out[, c(8:10)],
  function(x) {
    x <- as.factor(x)
    levels(x) <- c("typical", "outlier")
    x
  }
)

out1 <- data_sandy_out[which(data_sandy_out$out_label_Tur == "outlier"), 1]
out2 <- data_sandy_out[which(data_sandy_out$out_label_Cond == "outlier"), 1]
out3 <- data_sandy_out[which(data_sandy_out$out_label_Level == "outlier"), 1]
out <- unique(c(out1, out2, out3))

index <- which(data_sandy_out$Timestamp %in% out)

# neighbours
n1 <- which(data_sandy_out$Timestamp %in% out1) - 1
n1b <- which(data_sandy_out$Timestamp %in% out1) + 1
n2 <- which(data_sandy_out$Timestamp %in% out2) - 1
n2b <- which(data_sandy_out$Timestamp %in% out2) + 1
n3 <- which(data_sandy_out$Timestamp %in% out3) - 1
n3b <- which(data_sandy_out$Timestamp %in% out3) + 1

data_sandy_out$out_label_Tur <- as.character(data_sandy_out$out_label_Tur)
data_sandy_out$out_label_Tur[c(n1, n1b)] <- "neighbour"
data_sandy_out$out_label_Tur <- as.factor(data_sandy_out$out_label_Tur)

data_sandy_out$out_label_Cond <- as.character(data_sandy_out$out_label_Cond)
data_sandy_out$out_label_Cond[c(n2, n2b)] <- "neighbour"
data_sandy_out$out_label_Cond <- as.factor(data_sandy_out$out_label_Cond)

data_sandy_out$out_label_Level <- as.character(data_sandy_out$out_label_Level)
data_sandy_out$out_label_Level[c(n3, n3b)] <- "neighbour"
data_sandy_out$out_label_Level <- as.factor(data_sandy_out$out_label_Level)

data_sandy_out$label_all <- ifelse((data_sandy_out$out_label_Tur == "typical" &
                                      data_sandy_out$out_label_Cond == "typical" &
                                      data_sandy_out$out_label_Level == "typical"),
                                   "typical", "outlier"
)

data_sandy_out$label_all[unique(c(n1, n1b, n2, n2b, n3, n3b))] <- "neighbour"

colnames(data_sandy_out)[2:4] <- c("Level", "Conductivity", "Turbidity")








# Figure 2: Data preparation
## ---- transform_sandy
data_sandy_out$flag_TC <- ifelse((data_sandy_out$Cond <= 0 |
                                    data_sandy_out$Tur <= 0),
                                 "outlier", "typical"
)

data_sandy_out$flag_TCL <- ifelse((data_sandy_out$Cond <= 0 |
                                     data_sandy_out$Tur <= 0 |
                                     data_sandy_out$Level < 0), "outlier", "typical")
data_sandy_out <- drop_na(data_sandy_out)
sandy_trans <- oddwater::transform_data(data_sandy_out[, 1:4])
sandy_full <- left_join(sandy_trans, data_sandy_out[, -c(2:4)], by = "Timestamp")

# Figure 2
## ---- transformType
p1 <- plot_trans(x = "Turbidity", y = "Conductivity", sandy_full) + ggtitle("(a)")
p2 <- plot_trans(x = "log_Turbidity", y = "log_Conductivity", sandy_full) + ggtitle("(b)")
p3 <- plot_trans(x = "diff_log_Turbidity", y = "diff_log_Conductivity", sandy_full) + ggtitle("(c)")
p4 <- plot_trans(x = "der_log_bound_Turbidity", y = "der_log_bound_Conductivity", sandy_full) + ggtitle("(d)")
p5 <- plot_trans(x = "neg_der_log_bound_Turbidity", y = "pos_der_log_bound_Conductivity", sandy_full) + ggtitle("(e)")
p6 <- plot_trans(x = "rc_Turbidity", y = "rc_Conductivity", sandy_full) + ggtitle("(f)")
p7 <- plot_trans(x = "rdiff_Turbidity", y = "rdiff_Conductivity", sandy_full) + ggtitle("(g)")
p8 <- plot_trans(x = "rdifflog_Turbidity", y = "rdifflog_Conductivity", sandy_full) + ggtitle("(h)")

p <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
                       nrow = 2, ncol = 4, legend = "bottom",
                       common.legend = TRUE
)
ggsave("./fig/transform_type.png")






# Figure 3
## ---- Visualise_outlier_pairs_original_data
data_sandy_out$label_all <- ifelse((data_sandy_out$out_label_Tur == "typical" &
  data_sandy_out$out_label_Cond == "typical" &
  data_sandy_out$out_label_Level == "typical"),
"typical", "outlier"
)

data_sandy_out$label_all[unique(c(n1, n1b, n2, n2b, n3, n3b))] <- "neighbour"
colnames(data_sandy_out)[2:4] <- c("Level", "Conductivity", "Turbidity")

# (Top panel)
p1 <- plot_trans(x = "Turbidity", y = "Conductivity", data = data_sandy_out ) + ggtitle("(a)")
p2 <- plot_trans(x = "Turbidity", y = "Level", data = data_sandy_out ) + ggtitle("(b)")
p3 <- plot_trans(x = "Conductivity", y = "Level", data = data_sandy_out ) + ggtitle("(c)")


data <- sandy_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity",
  "neg_der_log_bound_Level", "label_all"
)]
colnames(data)[2:4] <- c("Trans_Tur", "Trans_Cond", "Trans_level")

# Bottom panel
p4 <- plot_trans(x = "Trans_Tur", y = "Trans_Cond", data = data ) + ggtitle("(d)")
p5 <- plot_trans(x = "Trans_Tur", y = "Trans_level", data = data ) + ggtitle("(e)")
p6 <- plot_trans(x = "Trans_Cond", y = "Trans_level", data = data ) + ggtitle("(f)")

p <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6,
                       nrow = 2, ncol = 3, legend = "bottom",
                       common.legend = TRUE
)

ggsave("./fig/Visualise_outlier_pairs_trans_data.png", plot = p)








# Figure 4
## ---- trans_demo_TCL
Tur <- data_sandy_out %>%
  select(Timestamp, Turbidity, out_label_Tur) %>%
  mutate(Var = "Turbidity")
Cond <- data_sandy_out %>%
  select(Timestamp, Conductivity, out_label_Cond) %>%
  mutate(Var = "Conductivity")
Level <- data_sandy_out %>%
  select(Timestamp, Level, out_label_Level) %>%
  mutate(Var = "Level")
colnames(Tur) <- colnames(Cond) <- colnames(Level) <- c("Time", "Value", "Type", "Var")


Trans_Tur <- sandy_full %>%
  select(Timestamp, neg_der_log_bound_Turbidity, out_label_Tur) %>%
  mutate(Var = "Trans_Tur")
Trans_Cond <- sandy_full %>%
  select(Timestamp, pos_der_log_bound_Conductivity, out_label_Tur) %>%
  mutate(Var = "Trans_cond")
Trans_Level <- sandy_full %>%
  select(Timestamp, neg_der_log_bound_Level, out_label_Tur) %>%
  mutate(Var = "Trans_Level")
colnames(Trans_Tur) <- colnames(Trans_Cond) <- colnames(Trans_Level) <- c("Time", "Value", "Type", "Var")


data <- bind_rows(Tur, Trans_Tur, Cond,Trans_Cond, Level, Trans_Level) %>%
  mutate(Type = as.factor(Type))
data$Var_o <- factor(data$Var, levels = c("Turbidity","Trans_Tur", "Conductivity", 
                                          "Trans_cond" , "Level", "Trans_Level" ))

out_data <- data %>% filter(Type %in% c("neighbour", "outlier"))

p <- ggplot(data, aes(
  x = Time, y = Value, color = Type,
  shape = Type, size = Type, alpha = Type
)) +
  geom_point() +
  scale_colour_manual(values = c(
    "typical" = "black", "outlier" = "#d95f02", "neighbour" = "#1b9e77"
  )) +
  scale_size_manual(values = c("outlier" = 4.5, "typical" = 1, "neighbour" = 3)) +
  scale_shape_manual(values = c("outlier" = 17, "typical" = 19, "neighbour" = 15)) +
  scale_alpha_manual(values = c("outlier" = 1, "typical" = 0.5, "neighbour" = 0.8)) +
  geom_point(data = out_data, aes(x = Time, y = Value)) +
  facet_grid(rows = vars(Var_o), scales = "free") +
  labs(x = "", y = "") +
  theme(legend.position = "bottom", legend.title = element_blank())

p <- ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = TRUE, legend = "bottom")
ggsave("./fig/trans_demo_TCL.png")






## ---- timegapSandy

p <- ggplot(sandy_full, aes(Timestamp, time)) +
  geom_point() +
  geom_line() +
  xlab("Time gap (in minutes)")
print(p)






# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: First  Derivative
# Variables: Turbidity, conductivity
# Site :  Sandy Creek
## ---- derivative_TC_sandy

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier"),
"outlier", "typical"
)

p1 <- plot_original(
  data = sandy_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3"))
p2 <- plot_original(
  data = sandy_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3"))

test_data <- sandy_full[, c(
  "Timestamp", "der_log_bound_Turbidity",
  "der_log_bound_Conductivity"
)] %>% drop_na()

p3 <- performance_plot(sandy_full, test_data,
  method = "NN-HD",
  flag = "flag_TC", title = "NN-HD", variable = "out_score"
)
p4 <- performance_plot(sandy_full, test_data,
  method = "KNN-AGG", flag = "flag_TC",
  title = "KNN-AGG", variable = "out_score"
)
p5 <- performance_plot(sandy_full, test_data,
  method = "KNN-SUM", flag = "flag_TC",
  title = "KNN-SUM", variable = "out_score"
)
p6 <- performance_plot(sandy_full, test_data,
  method = "LOF", flag = "flag_TC",
  variable = "out_score"
)
p7 <- performance_plot(sandy_full, test_data,
  method = "COF", flag = "flag_TC",
  variable = "out_score"
)
p8 <- performance_plot(sandy_full, test_data,
  method = "INFLO", flag = "flag_TC",
  variable = "out_score"
)
p9 <- performance_plot(sandy_full, test_data,
  method = "LDOF", flag = "flag_TC",
  variable = "out_score"
)
p10 <- performance_plot(sandy_full, test_data,
  method = "RKOF", flag = "flag_TC",
  variable = "out_score"
) +
  theme(legend.position = "bottom", legend.title = element_blank())

pD <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 10)
ggsave("./fig/derivative_TC_sandy.png", height = 13, width = 12)







# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: First  Derivative
# Variables: Turbidity, conductivity, Level
# Site :  Sandy Creek
## ---- derivative_TCL_sandy
sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier" |
  data_sandy_out$out_label_Level == "outlier"),
"outlier", "typical"
)

p1 <- plot_original(
  data = sandy_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3"))
p2 <- plot_original(
  data = sandy_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3"))
p11 <- plot_original(
  data = sandy_full, y = "Level",
  colour = "out_label_Level", y_label = "Level"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3"))

test_data <- sandy_full[, c(
  "Timestamp", "der_log_bound_Turbidity",
  "der_log_bound_Conductivity", "der_log_bound_Level"
)] %>% drop_na()

p3 <- performance_plot(sandy_full, test_data,
  method = "NN-HD",
  flag = "flag_TC", title = "NN-HD", variable = "out_score"
)
p4 <- performance_plot(sandy_full, test_data,
  method = "KNN-AGG", flag = "flag_TCL",
  title = "KNN-AGG", variable = "out_score"
)
p5 <- performance_plot(sandy_full, test_data,
  method = "KNN-SUM", flag = "flag_TCL",
  title = "KNN-SUM", variable = "out_score"
)
p6 <- performance_plot(sandy_full, test_data,
  method = "LOF", flag = "flag_TCL",
  variable = "out_score"
)
p7 <- performance_plot(sandy_full, test_data,
  method = "COF", flag = "flag_TCL",
  variable = "out_score"
)
p8 <- performance_plot(sandy_full, test_data,
  method = "INFLO", flag = "flag_TCL",
  variable = "out_score"
)
p9 <- performance_plot(sandy_full, test_data,
  method = "LDOF", flag = "flag_TCL",
  variable = "out_score"
)
p10 <- performance_plot(sandy_full, test_data,
  method = "RKOF", flag = "flag_TCL",
  variable = "out_score"
) +
  theme(legend.position = "bottom", legend.title = element_blank())

pC <- ggpubr::ggarrange(p1, p2, p11, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 11)

ggsave("./fig/derivative_TCL_sandy.png", height = 13, width = 12)






# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: One-sided First  Derivative
# Variables: Turbidity, conductivity
# Site:  Sandy Creek
## ---- one_sided_derivative_TC_sandy
sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier"),
"outlier", "typical"
)

p1 <- plot_original(
  data = sandy_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c(
    "typical" = "black", "outlier" = "red",
    "neighbour" = "lightsalmon3"
  ))
p2 <- plot_original(
  data = sandy_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c(
    "typical" = "black", "outlier" = "red",
    "neighbour" = "lightsalmon3"
  ))

test_data <- sandy_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity"
)] %>% drop_na()

p3 <- performance_plot(sandy_full, test_data,
  method = "NN-HD",
  flag = "flag_TC", title = "NN-HD", variable = "out_score"
)
p4 <- performance_plot(sandy_full, test_data,
  method = "KNN-AGG", flag = "flag_TC",
  title = "KNN-AGG", variable = "out_score"
)
p5 <- performance_plot(sandy_full, test_data,
  method = "KNN-SUM", flag = "flag_TC",
  title = "KNN-SUM", variable = "out_score"
)
p6 <- performance_plot(sandy_full, test_data,
  method = "LOF", flag = "flag_TC",
  variable = "out_score"
)
p7 <- performance_plot(sandy_full, test_data,
  method = "COF", flag = "flag_TC",
  variable = "out_score"
)
p8 <- performance_plot(sandy_full, test_data,
  method = "INFLO", flag = "flag_TC",
  variable = "out_score"
)
p9 <- performance_plot(sandy_full, test_data,
  method = "LDOF", flag = "flag_TC",
  variable = "out_score"
)
p10 <- performance_plot(sandy_full, test_data,
  method = "RKOF", flag = "flag_TC",
  variable = "out_score"
) +
  theme(legend.position = "bottom", legend.title = element_blank())

pB <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 10)
ggsave("./fig/one_sided_derivative_TC_sandy.png", height = 13, width = 12)






# Figure 6
# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: One-sided First  Derivative
# Variables: Turbidity, conductivity, Level
# Site :  Sandy Creek
## ---- one_sided_derivative_TCL_sandy

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier" |
  data_sandy_out$out_label_Level == "outlier"),
"outlier", "typical"
)

p1 <- plot_original(
  data = sandy_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3")) +
  ylab("(i) Turbidity") + theme(text = element_text(size = 11))

p2 <- plot_original(
  data = sandy_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3")) +
  ylab("(ii) Conductivity") + theme(text = element_text(size = 11))

p11 <- plot_original(
  data = sandy_full, y = "Level",
  colour = "out_label_Level", y_label = "Level"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red", "neighbour" = "lightsalmon3")) +
  ylab("(iii) Level") + theme(text = element_text(size = 11))

test_data <- sandy_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
)] %>% drop_na()

p3 <- performance_plot(sandy_full, test_data,
  method = "NN-HD",
  flag = "flag_TCL", title = "NN-HD", variable = "out_score"
) + ylab("(c) NN-HD") + theme(text = element_text(size = 11))
p4 <- performance_plot(sandy_full, test_data,
  method = "KNN-AGG", flag = "flag_TCL",
  title = "KNN-AGG", variable = "out_score"
) + ylab("(b) KNN-AGG") + theme(text = element_text(size = 11))
p5 <- performance_plot(sandy_full, test_data,
  method = "KNN-SUM", flag = "flag_TCL",
  title = "KNN-SUM", variable = "out_score"
) + ylab("(a) KNN-SUM") + theme(text = element_text(size = 11))
p6 <- performance_plot(sandy_full, test_data,
  method = "LOF", flag = "flag_TCL",
  variable = "out_score"
) + ylab("(d) LOF") + theme(text = element_text(size = 11))
p7 <- performance_plot(sandy_full, test_data,
  method = "COF", flag = "flag_TCL",
  variable = "out_score"
) + ylab("(e) COF") + theme(text = element_text(size = 11))
p8 <- performance_plot(sandy_full, test_data,
  method = "INFLO", flag = "flag_TCL",
  variable = "out_score"
) + ylab("(f) INFLO") + theme(text = element_text(size = 11))
p9 <- performance_plot(sandy_full, test_data,
  method = "LDOF", flag = "flag_TCL",
  variable = "out_score"
) + ylab("(g) LDOF") + theme(text = element_text(size = 11))
p10 <- performance_plot(sandy_full, test_data,
  method = "RKOF", flag = "flag_TCL",
  variable = "out_score"
) + ylab("(h) RKOF") + theme(text = element_text(size = 11))

p10 <- p10 + theme(legend.position = "bottom", legend.title = element_blank())

pA <- ggpubr::ggarrange(p1, p2, p11, p5, p4, p3, p6, p7, p8, p9, p10, nrow = 11)
ggsave("./fig/one_sided_derivative_TCL_sandy.png", height = 13, width = 12)







# Table 2
## ---- sandyTable

## Sandy Creek Derivative TCL
results_table <- NULL
sandy_full <- as_tibble(sandy_full)

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier" |
  data_sandy_out$out_label_Level == "outlier"),
"outlier", "typical"
)

test_data <- sandy_full[, c(
  "Timestamp", "der_log_bound_Turbidity",
  "der_log_bound_Conductivity", "der_log_bound_Level"
)] %>%
  drop_na()

r1 <- apply_ADmethod(sandy_full, test_data, method = "NN-HD", flag = "flag_TCL")
r2 <- apply_ADmethod(sandy_full, test_data, method = "KNN-AGG", flag = "flag_TCL")
r3 <- apply_ADmethod(sandy_full, test_data, method = "KNN-SUM", flag = "flag_TCL")
r4 <- apply_ADmethod(sandy_full, test_data, method = "INFLO", flag = "flag_TCL")
r5 <- apply_ADmethod(sandy_full, test_data, method = "COF", flag = "flag_TCL")
r6 <- apply_ADmethod(sandy_full, test_data, method = "LDOF", flag = "flag_TCL")
r7 <- apply_ADmethod(sandy_full, test_data, method = "LOF", flag = "flag_TCL")
r8 <- apply_ADmethod(sandy_full, test_data, method = "RKOF", flag = "flag_TCL")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C-L", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("First Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)



## Sandy Creek Derivative TC

sandy_full <- as_tibble(sandy_full)

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier"),
"outlier", "typical"
)


test_data <- sandy_full[, c("Timestamp", "der_log_bound_Turbidity", "der_log_bound_Conductivity")] %>% drop_na()

r1 <- apply_ADmethod(sandy_full, test_data, method = "NN-HD", flag = "flag_TC")
r2 <- apply_ADmethod(sandy_full, test_data, method = "KNN-AGG", flag = "flag_TC")
r3 <- apply_ADmethod(sandy_full, test_data, method = "KNN-SUM", flag = "flag_TC")
r4 <- apply_ADmethod(sandy_full, test_data, method = "INFLO", flag = "flag_TC")
r5 <- apply_ADmethod(sandy_full, test_data, method = "COF", flag = "flag_TC")
r6 <- apply_ADmethod(sandy_full, test_data, method = "LDOF", flag = "flag_TC")
r7 <- apply_ADmethod(sandy_full, test_data, method = "LOF", flag = "flag_TC")
r8 <- apply_ADmethod(sandy_full, test_data, method = "RKOF", flag = "flag_TC")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("First Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)


## Sandy Creek One sided Derivative TCL

sandy_full <- as_tibble(sandy_full)

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier" |
  data_sandy_out$out_label_Level == "outlier"),
"outlier", "typical"
)


test_data <- sandy_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
)] %>%
  drop_na()

r1 <- apply_ADmethod(sandy_full, test_data, method = "NN-HD", flag = "flag_TCL")
r2 <- apply_ADmethod(sandy_full, test_data, method = "KNN-AGG", flag = "flag_TCL")
r3 <- apply_ADmethod(sandy_full, test_data, method = "KNN-SUM", flag = "flag_TCL")
r4 <- apply_ADmethod(sandy_full, test_data, method = "INFLO", flag = "flag_TCL")
r5 <- apply_ADmethod(sandy_full, test_data, method = "COF", flag = "flag_TCL")
r6 <- apply_ADmethod(sandy_full, test_data, method = "LDOF", flag = "flag_TCL")
r7 <- apply_ADmethod(sandy_full, test_data, method = "LOF", flag = "flag_TCL")
r8 <- apply_ADmethod(sandy_full, test_data, method = "RKOF", flag = "flag_TCL")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C-L", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("One sided Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)


## Sandy Creek One sided Derivative TC

sandy_full <- as_tibble(sandy_full)

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier"),
"outlier", "typical"
)

test_data <- sandy_full[, c("Timestamp", "neg_der_log_bound_Turbidity", "pos_der_log_bound_Conductivity")] %>% drop_na()

r1 <- apply_ADmethod(sandy_full, test_data, method = "NN-HD", flag = "flag_TC")
r2 <- apply_ADmethod(sandy_full, test_data, method = "KNN-AGG", flag = "flag_TC")
r3 <- apply_ADmethod(sandy_full, test_data, method = "KNN-SUM", flag = "flag_TC")
r4 <- apply_ADmethod(sandy_full, test_data, method = "INFLO", flag = "flag_TC")
r5 <- apply_ADmethod(sandy_full, test_data, method = "COF", flag = "flag_TC")
r6 <- apply_ADmethod(sandy_full, test_data, method = "LDOF", flag = "flag_TC")
r7 <- apply_ADmethod(sandy_full, test_data, method = "LOF", flag = "flag_TC")
r8 <- apply_ADmethod(sandy_full, test_data, method = "RKOF", flag = "flag_TC")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("One sided Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)


## Sandy Creek Original TCL

sandy_full <- as_tibble(sandy_full)

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier" |
  data_sandy_out$out_label_Level == "outlier"),
"outlier", "typical"
)


test_data <- sandy_full[, c(
  "Timestamp", "Turbidity", "Conductivity", "Level"
)] %>%
  drop_na()

r1 <- apply_ADmethod(sandy_full, test_data, method = "NN-HD", flag = "flag_TCL")
r2 <- apply_ADmethod(sandy_full, test_data, method = "KNN-AGG", flag = "flag_TCL")
r3 <- apply_ADmethod(sandy_full, test_data, method = "KNN-SUM", flag = "flag_TCL")
r4 <- apply_ADmethod(sandy_full, test_data, method = "INFLO", flag = "flag_TCL")
r5 <- apply_ADmethod(sandy_full, test_data, method = "COF", flag = "flag_TCL")
r6 <- apply_ADmethod(sandy_full, test_data, method = "LDOF", flag = "flag_TCL")
r7 <- apply_ADmethod(sandy_full, test_data, method = "LOF", flag = "flag_TCL")
r8 <- apply_ADmethod(sandy_full, test_data, method = "RKOF", flag = "flag_TCL")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C-L", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("Original series", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)


## Sandy Creek Original TC

sandy_full <- as_tibble(sandy_full)

sandy_full$out_label_Tur <- ifelse(sandy_full$out_label_Tur == "outlier", "outlier", "typical")
sandy_full$out_label_Cond <- ifelse(sandy_full$out_label_Cond == "outlier", "outlier", "typical")
sandy_full$out_label_Level <- ifelse(sandy_full$out_label_Level == "outlier", "outlier", "typical")

sandy_full$label_all <- ifelse((data_sandy_out$out_label_Tur == "outlier" |
  data_sandy_out$out_label_Cond == "outlier"),
"outlier", "typical"
)

test_data <- sandy_full[, c(
  "Timestamp", "Turbidity", "Conductivity"
)] %>%
  drop_na()

r1 <- apply_ADmethod(sandy_full, test_data, method = "NN-HD", flag = "flag_TC")
r2 <- apply_ADmethod(sandy_full, test_data, method = "KNN-AGG", flag = "flag_TC")
r3 <- apply_ADmethod(sandy_full, test_data, method = "KNN-SUM", flag = "flag_TC")
r4 <- apply_ADmethod(sandy_full, test_data, method = "INFLO", flag = "flag_TC")
r5 <- apply_ADmethod(sandy_full, test_data, method = "COF", flag = "flag_TC")
r6 <- apply_ADmethod(sandy_full, test_data, method = "LDOF", flag = "flag_TC")
r7 <- apply_ADmethod(sandy_full, test_data, method = "LOF", flag = "flag_TC")
r8 <- apply_ADmethod(sandy_full, test_data, method = "RKOF", flag = "flag_TC")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("Original series", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)


# order : Der_TCL, Der_TC, OS_TCL, OS_TC
load(file = "data/sandy_OS_der.rda")
sandy_der <- rbind(rt[9:16, ], rt[1:8, ])

load(file = "data/sandy_OS.rda")
sandy_OS <- rbind(rt[9:16, ], rt[1:8, ])

load(file = "data/sandy_Original.rda")
sandy_Original <- rbind(rt[9:16, ], rt[1:8, ])

sandy_runtime <- rbind(sandy_der, sandy_OS, sandy_Original)

runtime_min <- round(sandy_runtime$min_t, 2)
runtime_max <- round(sandy_runtime$max_t, 2)
runtime_mu <- round(sandy_runtime$mu_t, 2)
results_table <- cbind(results_table, runtime_min, runtime_mu, runtime_max)

results_table[, 9] <- round(results_table[, 9], 1)
results_table[, 13:15] <- round(results_table[, 13:15], 0)
colnames(results_table) <- c(
  "Variables", "Transformation", "Method", "TN", "FN", "FP",
  "TP", "Accuracy", "ER", "GM", "OP", "PPV", "NPV", "min_t",
  "mu_t", "max_t"
)
results_table <- cbind(i = 1:nrow(results_table), results_table %>% arrange(desc(OP)))
kable(results_table[, -10], "latex", caption = "Performance metrics of outlier
      detection algorithms performed on multivariate water-quality time series
      data (T, turbidity; C, conductivity; L, river level) from in situ sensors at Sandy Creek, arranged in descending 
      order of OP values. See Sections 2.7-8 for performance metric codes and details.", booktabs = T) %>%
  kable_styling(latex_options = "scale_down")




## ---- draw_bestoutput_sandy

sandy_full <- as_tibble(sandy_full)

sandy_full$label_all <- ifelse((sandy_full$out_label_Level == "typical" &
  sandy_full$out_label_Cond == "typical" &
  sandy_full$out_label_Tur == "typical"),
"typical", "outlier"
)

test_data <- sandy_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
)] %>%
  drop_na()


r2 <- apply_ADmethod(sandy_full, test_data, method = "KNN-SUM", flag = "flag_TCL", adjustment = TRUE)
ind <- r2$outlier_index

out_data <- sandy_full
variable <- "Turbidity"
tur <- ind[(ind$var == variable), 1]
index <- which(data_sandy_out$Timestamp %in% tur)
out_data$col <- "typical"
out_data$col[index] <- "outlier"

p1 <- plot_output(out_data,
  variable = "Turbidity", label_type = "out_label_Tur",
  title = variable
) + ggtitle("(a)")

variable <- "Conductivity"
con <- ind[(ind$var == variable), 1]
index <- which(data_sandy_out$Timestamp %in% con)
out_data$col <- "typical"
out_data$col[index] <- "outlier"

p2 <- plot_output(out_data,
  variable = "Conductivity", label_type = "out_label_Cond",
  title = variable
) + ggtitle("(b)")

variable <- "Level"
lev <- ind[(ind$var == variable), 1]
index <- which(data_sandy_out$Timestamp %in% lev)
out_data$col <- "typical"
out_data$col[index] <- "outlier"

p3 <- plot_output(out_data, variable = "Level", label_type = "out_label_Level", title = variable) + ggtitle("(c)")


p <- ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = TRUE, legend = "bottom")
ggsave("./fig/sandy_outlier_plot.png")










###########################################################
#### Analysis of data from Pioneer River
###########################################################

## ---- filter_outliers_from_anormalies_pioneer
data("data_pioneer_anom")
## Change the format and the class of the existing variable
data_pioneer_anom$Timestamp <- lubridate::dmy_hm(data_pioneer_anom$Timestamp)
# Turbidity
data_pioneer_anom$out_type_Tur <- ifelse(data_pioneer_anom$type_Tur %in%
  c("B", "C", "E", "H", "L"), "0",
as.character(data_pioneer_anom$type_Tur)
)
# Conductivity
data_pioneer_anom$out_type_Cond <- ifelse(data_pioneer_anom$type_Cond %in%
  c("B", "C", "E", "H", "L"), "0",
as.character(data_pioneer_anom$type_Cond)
)
# Level
data_pioneer_anom$out_type_Level <- ifelse(data_pioneer_anom$type_Level %in%
  c("B", "C", "E", "H", "L"), "0",
as.character(data_pioneer_anom$type_Level)
)
# creat a new dataset including only outliers
data_pioneer_out <- data_pioneer_anom[, c("Timestamp", "Level", "Cond", "Tur", "out_type_Tur", "out_type_Cond", "out_type_Level")]
data_pioneer_out$out_label_Tur <- ifelse((data_pioneer_out$out_type_Tur == "0"), 0, 1)
data_pioneer_out$out_label_Cond <- ifelse((data_pioneer_out$out_type_Cond == "0"), 0, 1)
data_pioneer_out$out_label_Level <- ifelse((data_pioneer_out$out_type_Level == "0"), 0, 1)

data_pioneer_out[, c(8:10)] <- lapply(
  data_pioneer_out[, c(8:10)],
  function(x) {
    x <- as.factor(x)
    levels(x) <- c("typical", "outlier")
    x
  }
)

out1 <- data_pioneer_out[which(data_pioneer_out$out_label_Tur == "outlier"), 1]
out2 <- data_pioneer_out[which(data_pioneer_out$out_label_Cond == "outlier"), 1]
out3 <- data_pioneer_out[which(data_pioneer_out$out_label_Level == "outlier"), 1]
out <- unique(c(out1, out2, out3))

index <- which(data_pioneer_out$Timestamp %in% out)


## neighbours

crit_out_t <- out1
neg_c <- data_pioneer_out[which(data_pioneer_out$Cond <= 0), 1]
crit_out_c <- out2[ which(!(out2 %in% neg_c))]
crit_out_l <- out3













# Figure 8 (Top panel)
## ---- Visualise_outlier_pairs_original_data_pioneer
colnames(data_pioneer_out)[2:4] <- c("Level", "Conductivity", "Turbidity")

data_pioneer_out$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
                                        data_pioneer_out$out_label_Cond == "typical" &
                                        data_pioneer_out$out_label_Level == "typical"),
                                     "typical", "outlier"
)

data <- data_pioneer_out

n1 <- which(data$Timestamp %in% crit_out_t) - 1
n1b <- which(data$Timestamp %in% crit_out_t) + 1
n2 <- which(data$Timestamp %in% crit_out_c) - 1
n2b <- which(data$Timestamp %in% crit_out_c) + 1
n3 <- which(data$Timestamp %in% crit_out_l) - 1
n3b <- which(data$Timestamp %in% crit_out_l) + 1

data$label_all[unique(c(n1, n1b, n2, n2b, n3, n3b))] <- "neighbour"

# (Top panel)
p1 <- plot_trans(x = "Turbidity", y = "Conductivity", data = data ) + ggtitle("(a)")
p2 <- plot_trans(x = "Turbidity", y = "Level", data = data ) + ggtitle("(b)")
p3 <- plot_trans(x = "Conductivity", y = "Level", data = data ) + ggtitle("(c)")



data_pioneer_out$flag_TC <- ifelse((data_pioneer_out$Cond <= 0 |
                                      data_pioneer_out$Tur <= 0),
                                   "outlier", "typical"
)

data_pioneer_out$flag_TCL <- ifelse((data_pioneer_out$Cond <= 0 |
                                       data_pioneer_out$Tur <= 0 |
                                       data_pioneer_out$Level < 0), "outlier", "typical")

data_pioneer_out <- drop_na(data_pioneer_out)
pioneer_trans <- oddwater::transform_data(data_pioneer_out[, 1:4])
pioneer_full <- left_join(pioneer_trans, data_pioneer_out[, -c(2:4)], by = "Timestamp")


# Bottom panel
data <- pioneer_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity",
  "neg_der_log_bound_Level", "label_all"
)]

n1 <- which(data$Timestamp %in% crit_out_t)-1
n1b <- which(data$Timestamp %in% crit_out_t)+1
n2 <- which(data$Timestamp %in% crit_out_c)-1
n2b <- which(data$Timestamp %in% crit_out_c)+1
n3 <- which(data$Timestamp %in% crit_out_l)-1
n3b <- which(data$Timestamp %in% crit_out_l)+1

data$label_all[unique(c(n1, n1b, n2, n2b, n3, n3b))] <- "neighbour"

colnames(data)[2:4] <- c("Trans_Tur", "Trans_Cond", "Trans_level")

data <- data[!is.infinite(rowSums(data[
  ,
  c("Trans_Tur", "Trans_Cond", "Trans_level"),
  ])), ]

p4 <- plot_trans(x = "Trans_Tur", y = "Trans_Cond", data = data ) +
  ggtitle("(d)")
p5 <- plot_trans(x = "Trans_Tur", y = "Trans_level", data = data ) + ggtitle("(e)")
p6 <- plot_trans(x = "Trans_Cond", y = "Trans_level", data = data ) + ggtitle("(f)")

p <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6,
                       nrow = 2, ncol = 3, legend = "bottom",
                       common.legend = TRUE
)

ggsave("./fig/Visualise_outlier_pairs_trans_data_pioneer.png", plot = p)









# Figure 7
## ---- Visualise_outlier_pioneer

p1 <- plot_original(
  data = data_pioneer_out, y = "Tur", colour = "out_label_Tur",
  y_label = "(a) Turbidity"
)
p2 <- plot_original(
  data = data_pioneer_out, y = "Cond", colour = "out_label_Cond",
  y_label = "(b) Conductivity"
)
p3 <- plot_original(
  data = data_pioneer_out, y = "Level", colour = "out_label_Level",
  y_label = "(c) Level"
)

p <- ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = TRUE, legend = "bottom")
ggsave("./fig/Visualise_outlier_pioneer.png")



# Figure 9
## ---- trans_demo_TCL_pioneer






# Figure 4
## ---- trans_demo_TCL
Tur <- data_sandy_out %>%
  select(Timestamp, Turbidity, out_label_Tur) %>%
  mutate(Var = "Turbidity")
Cond <- data_sandy_out %>%
  select(Timestamp, Conductivity, out_label_Cond) %>%
  mutate(Var = "Conductivity")
Level <- data_sandy_out %>%
  select(Timestamp, Level, out_label_Level) %>%
  mutate(Var = "Level")
colnames(Tur) <- colnames(Cond) <- colnames(Level) <- c("Time", "Value", "Type", "Var")


Trans_Tur <- sandy_full %>%
  select(Timestamp, neg_der_log_bound_Turbidity, out_label_Tur) %>%
  mutate(Var = "Trans_Tur")
Trans_Cond <- sandy_full %>%
  select(Timestamp, pos_der_log_bound_Conductivity, out_label_Tur) %>%
  mutate(Var = "Trans_cond")
Trans_Level <- sandy_full %>%
  select(Timestamp, neg_der_log_bound_Level, out_label_Tur) %>%
  mutate(Var = "Trans_Level")
colnames(Trans_Tur) <- colnames(Trans_Cond) <- colnames(Trans_Level) <- c("Time", "Value", "Type", "Var")


data <- bind_rows(Tur, Trans_Tur, Cond,Trans_Cond, Level, Trans_Level) %>%
  mutate(Type = as.factor(Type))
data$Var_o <- factor(data$Var, levels = c("Turbidity","Trans_Tur", "Conductivity", 
                                          "Trans_cond" , "Level", "Trans_Level" ))

out_data <- data %>% filter(Type %in% c("neighbour", "outlier"))

p <- ggplot(data, aes(
  x = Time, y = Value, color = Type,
  shape = Type, size = Type, alpha = Type
)) +
  geom_point() +
  scale_colour_manual(values = c(
    "typical" = "black", "outlier" = "#d95f02", "neighbour" = "#1b9e77"
  )) +
  scale_size_manual(values = c("outlier" = 4.5, "typical" = 1, "neighbour" = 3)) +
  scale_shape_manual(values = c("outlier" = 17, "typical" = 19, "neighbour" = 15)) +
  scale_alpha_manual(values = c("outlier" = 1, "typical" = 0.5, "neighbour" = 0.8)) +
  geom_point(data = out_data, aes(x = Time, y = Value)) +
  facet_grid(rows = vars(Var_o), scales = "free") +
  labs(x = "", y = "") +
  theme(legend.position = "bottom", legend.title = element_blank())

p <- ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = TRUE, legend = "bottom")
ggsave("./fig/trans_demo_TCL.png")









plot_data <- pioneer_full

# neighbours
n1 <- which(plot_data$Timestamp %in% crit_out_t) - 1
n1b <- which(plot_data$Timestamp %in% crit_out_t) + 1
n2 <- which(plot_data$Timestamp %in% crit_out_c) - 1
n2b <- which(plot_data$Timestamp %in% crit_out_c) + 1
n3 <- which(plot_data$Timestamp %in% crit_out_l) - 1
n3b <- which(plot_data$Timestamp %in% crit_out_l) + 1

plot_data$out_label_Tur <- as.character(plot_data$out_label_Tur)
plot_data$out_label_Tur[c(n1, n1b)] <- "neighbour"
plot_data$out_label_Tur <- as.factor(plot_data$out_label_Tur)


plot_data$out_label_Cond <- as.character(plot_data$out_label_Cond)
plot_data$out_label_Cond[c(n2, n2b)] <- "neighbour"
plot_data$out_label_Cond <- as.factor(plot_data$out_label_Cond)


plot_data$out_label_Level <- as.character(plot_data$out_label_Level)
plot_data$out_label_Level[c(n3, n3b)] <- "neighbour"
plot_data$out_label_Level <- as.factor(plot_data$out_label_Level)

p1 <- plot_original(
  data = plot_data, y = "neg_der_log_bound_Turbidity",
  colour = "out_label_Tur", y_label = "(a) Turbidity"
)
p2 <- plot_original(
  data = plot_data, y = "pos_der_log_bound_Conductivity",
  colour = "out_label_Cond", y_label = "(b) Conductivity"
)
p3 <- plot_original(
  data = plot_data, y = "neg_der_log_bound_Level",
  colour = "out_label_Level", y_label = "(c) Level"
)


p <- ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = TRUE, legend = "bottom")
ggsave("./fig/trans_demo_TCL_pioneer.png")














# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: First  Derivative
# Variables: Turbidity, conductivity
# Site: Pioneer River
## ---- derivative_TC_pioneer

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical"),
"typical", "outlier"
)

p1 <- plot_original(
  data = pioneer_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))
p2 <- plot_original(
  data = pioneer_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))

test_data <- pioneer_full[, c(
  "Timestamp", "der_log_bound_Turbidity",
  "der_log_bound_Conductivity"
)] %>% drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c("der_log_bound_Turbidity", "der_log_bound_Conductivity"),
])), ]

p3 <- performance_plot(pioneer_full, test_data,
  method = "NN-HD",
  flag = "flag_TC", title = "HDo", variable = "out_score"
)

p4 <- performance_plot(pioneer_full, test_data,
  method = "KNN-AGG", flag = "flag_TC",
  title = "KNN-AGG", variable = "out_score"
)
p5 <- performance_plot(pioneer_full, test_data,
  method = "KNN-SUM", flag = "flag_TC",
  title = "KNN-SUM", variable = "out_score"
)
p6 <- performance_plot(pioneer_full, test_data,
  method = "LOF", flag = "flag_TC",
  variable = "out_score"
)
p7 <- performance_plot(pioneer_full, test_data,
  method = "COF", flag = "flag_TC",
  variable = "out_score"
)
p8 <- performance_plot(pioneer_full, test_data,
  method = "INFLO", flag = "flag_TC",
  variable = "out_score"
)
p9 <- performance_plot(pioneer_full, test_data,
  method = "LDOF", flag = "flag_TC",
  variable = "out_score"
)
p10 <- performance_plot(pioneer_full, test_data,
  method = "RKOF", flag = "flag_TC",
  variable = "out_score"
) +
  theme(legend.position = "bottom", legend.title = element_blank())

pA <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 10)
ggsave("./fig/derivative_TC__pioneer.png", height = 13, width = 12)







# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: First  Derivative
# Variables: Turbidity, conductivity, Level
# Site: Pioneer River
## ---- derivative_TCL_pioneer

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical" &
  data_pioneer_out$out_label_Level == "typical"),
"typical", "outlier"
)

p1 <- plot_original(
  data = pioneer_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))
p2 <- plot_original(
  data = pioneer_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))
p11 <- plot_original(
  data = pioneer_full, y = "Level",
  colour = "out_label_Level", y_label = "Level"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))

test_data <- pioneer_full[, c(
  "Timestamp", "der_log_bound_Turbidity",
  "der_log_bound_Conductivity", "der_log_bound_Level"
)] %>% drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c("der_log_bound_Turbidity", "der_log_bound_Conductivity", "der_log_bound_Level"),
])), ]

p3 <- performance_plot(pioneer_full, test_data,
  method = "NN-HD",
  flag = "flag_TC", title = "HDo", variable = "out_score"
)
p4 <- performance_plot(pioneer_full, test_data,
  method = "KNN-AGG", flag = "flag_TCL",
  title = "KNN-AGG", variable = "out_score"
)
p5 <- performance_plot(pioneer_full, test_data,
  method = "KNN-SUM", flag = "flag_TCL",
  title = "KNN-SUM", variable = "out_score"
)
p6 <- performance_plot(pioneer_full, test_data,
  method = "LOF", flag = "flag_TCL",
  variable = "out_score"
)
p7 <- performance_plot(pioneer_full, test_data,
  method = "COF", flag = "flag_TCL",
  variable = "out_score"
)
p8 <- performance_plot(pioneer_full, test_data,
  method = "INFLO", flag = "flag_TCL",
  variable = "out_score"
)
p9 <- performance_plot(pioneer_full, test_data,
  method = "LDOF", flag = "flag_TCL",
  variable = "out_score"
)
p10 <- performance_plot(pioneer_full, test_data,
  method = "RKOF", flag = "flag_TCL",
  variable = "out_score"
) +
  theme(legend.position = "bottom", legend.title = element_blank())

pC <- ggpubr::ggarrange(p1, p2, p11, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 11)
ggsave("./fig/derivative_TCL_pioneer.png", height = 13, width = 12)






# Figure 10
# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: One-sided First  Derivative
# Variables: Turbidity, conductivity
# Site: Pioneer River
## ---- one_sided_derivative_TC_pioneer

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical"),
"typical", "outlier"
)


p1 <- plot_original(
  data = pioneer_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red")) +
  ylab("(i) Turbidity") + theme(text = element_text(size = 11))

p2 <- plot_original(
  data = pioneer_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red")) +
  ylab("(ii) Conductivity") + theme(text = element_text(size = 11))


test_data <- pioneer_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity"
)] %>% drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c("neg_der_log_bound_Turbidity", "pos_der_log_bound_Conductivity"),
])), ]

p3 <- performance_plot(pioneer_full, test_data,
  method = "NN-HD",
  flag = "flag_TC", title = "HDo", variable = "out_score"
) + ylab("(a) NN-HD") + theme(text = element_text(size = 11))
p4 <- performance_plot(pioneer_full, test_data,
  method = "KNN-AGG", flag = "flag_TC",
  title = "KNN-AGG", variable = "out_score"
) + ylab("(b) KNN-AGG") + theme(text = element_text(size = 11))
p5 <- performance_plot(pioneer_full, test_data,
  method = "KNN-SUM", flag = "flag_TC",
  title = "KNN-SUM", variable = "out_score"
) + ylab("(c) KNN-SUM") + theme(text = element_text(size = 11))
p6 <- performance_plot(pioneer_full, test_data,
  method = "LOF", flag = "flag_TC",
  variable = "out_score"
) + ylab("(d) LOF") + theme(text = element_text(size = 11))
p7 <- performance_plot(pioneer_full, test_data,
  method = "COF", flag = "flag_TC",
  variable = "out_score"
) + ylab("(e) COF") + theme(text = element_text(size = 11))
p8 <- performance_plot(pioneer_full, test_data,
  method = "INFLO", flag = "flag_TC",
  variable = "out_score"
) + ylab("(f) INFLO") + theme(text = element_text(size = 11))
p9 <- performance_plot(pioneer_full, test_data,
  method = "LDOF", flag = "flag_TC",
  variable = "out_score"
) + ylab("(g) LDOF") + theme(text = element_text(size = 11))
p10 <- performance_plot(pioneer_full, test_data,
  method = "RKOF", flag = "flag_TC",
  variable = "out_score"
) + ylab("(h) RKOF") + theme(text = element_text(size = 11))


p10 <- p10 +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(breaks = NULL)


pB <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 10)
ggsave("./fig/one_sided_derivative_TC_pioneer.png", height = 13, width = 12)







# Classification of outlier scores produced from different algorithms as TN, TP, FN, FP
# Tranformation: One-sided First  Derivative
# Variables: Turbidity, conductivity, Level
# Site: Pioneer River
## ---- one_sided_derivative_TCL_pioneer

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical" &
  data_pioneer_out$out_label_Level == "typical"),
"typical", "outlier"
)


p1 <- plot_original(
  data = pioneer_full, y = "Turbidity",
  colour = "out_label_Tur", y_label = "Turb"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))
p2 <- plot_original(
  data = pioneer_full, y = "Conductivity",
  colour = "out_label_Cond", y_label = "Cond"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))
p11 <- plot_original(
  data = pioneer_full, y = "Level",
  colour = "out_label_Level", y_label = "Level"
) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_colour_manual(values = c("typical" = "black", "outlier" = "red"))

test_data <- pioneer_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
)] %>% drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c(
    "neg_der_log_bound_Turbidity",
    "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
  ),
])), ]

p3 <- performance_plot(pioneer_full, test_data,
  method = "NN-HD",
  flag = "flag_TC", title = "HDo", variable = "out_score"
)
p4 <- performance_plot(pioneer_full, test_data,
  method = "KNN-AGG", flag = "flag_TCL",
  title = "KNN-AGG", variable = "out_score"
)
p5 <- performance_plot(pioneer_full, test_data,
  method = "KNN-SUM", flag = "flag_TCL",
  title = "KNN-SUM", variable = "out_score"
)
p6 <- performance_plot(pioneer_full, test_data,
  method = "LOF", flag = "flag_TCL",
  variable = "out_score"
)
p7 <- performance_plot(pioneer_full, test_data,
  method = "COF", flag = "flag_TCL",
  variable = "out_score"
)
p8 <- performance_plot(pioneer_full, test_data,
  method = "INFLO", flag = "flag_TCL",
  variable = "out_score"
)
p9 <- performance_plot(pioneer_full, test_data,
  method = "LDOF", flag = "flag_TCL",
  variable = "out_score"
)
p10 <- performance_plot(pioneer_full, test_data,
  method = "RKOF", flag = "flag_TCL",
  variable = "out_score"
) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(breaks = NULL)

pA <- ggpubr::ggarrange(p1, p2, p11, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 11)
ggsave("./fig/one_sided_derivative_TCL_pioneer.png", height = 13, width = 12)




# Table 3
## ---- pioneerTable

## Pioneer Derivative TCL
results_table <- NULL
pioneer_full <- as_tibble(pioneer_full)

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical" &
  data_pioneer_out$out_label_Level == "typical"),
"typical", "outlier"
)



test_data <- pioneer_full[, c(
  "Timestamp", "der_log_bound_Turbidity",
  "der_log_bound_Conductivity", "der_log_bound_Level"
)] %>%
  drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c(
    "der_log_bound_Turbidity",
    "der_log_bound_Conductivity", "der_log_bound_Level"
  ),
])), ]

r1 <- apply_ADmethod(pioneer_full, test_data, method = "NN-HD", flag = "flag_TCL")
r2 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-AGG", flag = "flag_TCL")
r3 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-SUM", flag = "flag_TCL")
r4 <- apply_ADmethod(pioneer_full, test_data, method = "INFLO", flag = "flag_TCL")
r5 <- apply_ADmethod(pioneer_full, test_data, method = "COF", flag = "flag_TCL")
r6 <- apply_ADmethod(pioneer_full, test_data, method = "LDOF", flag = "flag_TCL")
r7 <- apply_ADmethod(pioneer_full, test_data, method = "LOF", flag = "flag_TCL")
r8 <- apply_ADmethod(pioneer_full, test_data, method = "RKOF", flag = "flag_TCL")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C-L", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("First Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)


## Pioneer Derivative TC

pioneer_full <- as_tibble(pioneer_full)
pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical"),
"typical", "outlier"
)
test_data <- pioneer_full[, c(
  "Timestamp", "der_log_bound_Turbidity",
  "der_log_bound_Conductivity"
)] %>% drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c(
    "der_log_bound_Turbidity",
    "der_log_bound_Conductivity"
  ),
])), ]

r1 <- apply_ADmethod(pioneer_full, test_data, method = "NN-HD", flag = "flag_TC")
r2 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-AGG", flag = "flag_TC")
r3 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-SUM", flag = "flag_TC")
r4 <- apply_ADmethod(pioneer_full, test_data, method = "INFLO", flag = "flag_TC")
r5 <- apply_ADmethod(pioneer_full, test_data, method = "COF", flag = "flag_TC")
r6 <- apply_ADmethod(pioneer_full, test_data, method = "LDOF", flag = "flag_TC")
r7 <- apply_ADmethod(pioneer_full, test_data, method = "LOF", flag = "flag_TC")
r8 <- apply_ADmethod(pioneer_full, test_data, method = "RKOF", flag = "flag_TC")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("First Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)

## Pioneer One sided Derivative TCL

pioneer_full <- as_tibble(pioneer_full)

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical" &
  data_pioneer_out$out_label_Level == "typical"),
"typical", "outlier"
)
test_data <- pioneer_full[, c(
  "Timestamp", "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
)] %>%
  drop_na()
test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c(
    "neg_der_log_bound_Turbidity",
    "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
  ),
])), ]

r1 <- apply_ADmethod(pioneer_full, test_data, method = "NN-HD", flag = "flag_TCL")
r2 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-AGG", flag = "flag_TCL")
r3 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-SUM", flag = "flag_TCL")
r4 <- apply_ADmethod(pioneer_full, test_data, method = "INFLO", flag = "flag_TCL")
r5 <- apply_ADmethod(pioneer_full, test_data, method = "COF", flag = "flag_TCL")
r6 <- apply_ADmethod(pioneer_full, test_data, method = "LDOF", flag = "flag_TCL")
r7 <- apply_ADmethod(pioneer_full, test_data, method = "LOF", flag = "flag_TCL")
r8 <- apply_ADmethod(pioneer_full, test_data, method = "RKOF", flag = "flag_TCL")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C-L", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("One sided Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)

## Pioneer One sided Derivative TC

pioneer_full <- as_tibble(pioneer_full)

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical"),
"typical", "outlier"
)
test_data <- pioneer_full[, c("Timestamp", "neg_der_log_bound_Turbidity", "pos_der_log_bound_Conductivity")] %>% drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c(
    "neg_der_log_bound_Turbidity",
    "pos_der_log_bound_Conductivity"
  ),
])), ]
r1 <- apply_ADmethod(pioneer_full, test_data, method = "NN-HD", flag = "flag_TC")
r2 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-AGG", flag = "flag_TC")
r3 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-SUM", flag = "flag_TC")
r4 <- apply_ADmethod(pioneer_full, test_data, method = "INFLO", flag = "flag_TC")
r5 <- apply_ADmethod(pioneer_full, test_data, method = "COF", flag = "flag_TC")
r6 <- apply_ADmethod(pioneer_full, test_data, method = "LDOF", flag = "flag_TC")
r7 <- apply_ADmethod(pioneer_full, test_data, method = "LOF", flag = "flag_TC")
r8 <- apply_ADmethod(pioneer_full, test_data, method = "RKOF", flag = "flag_TC")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("One sided Derivative", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)


## Pioneer Original TCL

pioneer_full <- as_tibble(pioneer_full)

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical" &
  data_pioneer_out$out_label_Level == "typical"),
"typical", "outlier"
)
test_data <- pioneer_full[, c(
  "Timestamp", "Turbidity",
  "Conductivity", "Level"
)] %>%
  drop_na()
test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c(
    "Turbidity",
    "Conductivity", "Level"
  ),
])), ]

r1 <- apply_ADmethod(pioneer_full, test_data, method = "NN-HD", flag = "flag_TCL")
r2 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-AGG", flag = "flag_TCL")
r3 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-SUM", flag = "flag_TCL")
r4 <- apply_ADmethod(pioneer_full, test_data, method = "INFLO", flag = "flag_TCL")
r5 <- apply_ADmethod(pioneer_full, test_data, method = "COF", flag = "flag_TCL")
r6 <- apply_ADmethod(pioneer_full, test_data, method = "LDOF", flag = "flag_TCL")
r7 <- apply_ADmethod(pioneer_full, test_data, method = "LOF", flag = "flag_TCL")
r8 <- apply_ADmethod(pioneer_full, test_data, method = "RKOF", flag = "flag_TCL")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C-L", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("Original Series", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)

## Pioneer Original TC

pioneer_full <- as_tibble(pioneer_full)

pioneer_full$label_all <- ifelse((data_pioneer_out$out_label_Tur == "typical" &
  data_pioneer_out$out_label_Cond == "typical"),
"typical", "outlier"
)
test_data <- pioneer_full[, c("Timestamp", "Turbidity", "Conductivity")] %>% drop_na()

test_data <- test_data[!is.infinite(rowSums(test_data[
  ,
  c(
    "Turbidity", "Conductivity"
  ),
])), ]
r1 <- apply_ADmethod(pioneer_full, test_data, method = "NN-HD", flag = "flag_TC")
r2 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-AGG", flag = "flag_TC")
r3 <- apply_ADmethod(pioneer_full, test_data, method = "KNN-SUM", flag = "flag_TC")
r4 <- apply_ADmethod(pioneer_full, test_data, method = "INFLO", flag = "flag_TC")
r5 <- apply_ADmethod(pioneer_full, test_data, method = "COF", flag = "flag_TC")
r6 <- apply_ADmethod(pioneer_full, test_data, method = "LDOF", flag = "flag_TC")
r7 <- apply_ADmethod(pioneer_full, test_data, method = "LOF", flag = "flag_TC")
r8 <- apply_ADmethod(pioneer_full, test_data, method = "RKOF", flag = "flag_TC")

results <- rbind(
  r1$measures, r2$measures, r3$measures,
  r4$measures, r5$measures, r6$measures,
  r7$measures, r8$measures
)
TimeSeries <- rep("T-C", 8)
Method <- c(
  "NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF",
  "LDOF", "LOF", "RKOF"
)
Tranformation <- rep("Original Series", 8)

results <- data.frame(TimeSeries, Tranformation, Method, results[, c(1:6, 12:15)])
results_table <- rbind(results_table, results)



# order : Der_TCL, Der_TC, OS_TCL, OS_TC
load(file = "data/pioneer_OS_der.rda")
pioneer_der <- rbind(rt[9:16, ], rt[1:8, ])

load(file = "data/pioneer_OS.rda")
pioneer_OS <- rbind(rt[9:16, ], rt[1:8, ])

load(file = "data/pioneer_Original.rda")
pioneer_Original <- rbind(rt[9:16, ], rt[1:8, ])

pioneer_runtime <- rbind(pioneer_der, pioneer_OS, pioneer_Original)

runtime_min <- round(pioneer_runtime$min_t, 2)
runtime_max <- round(pioneer_runtime$max_t, 2)
runtime_mu <- round(pioneer_runtime$mu_t, 2)
results_table <- cbind(results_table, runtime_min, runtime_mu, runtime_max)


results_table[, 9] <- round(results_table[, 9], 1)
results_table[, 13:15] <- round(results_table[, 13:15], 0)

colnames(results_table) <- c(
  "Variables", "Transformation", "Method", "TN", "FN", "FP",
  "TP", "Accuracy", "ER", "GM", "OP", "PPV", "NPV", "min_t",
  "mu_t", "max_t"
)

results_table <- cbind(i = 1:nrow(results_table), results_table %>% arrange(desc(OP)))
kable(results_table[, -10], "latex", caption = "Performance metrics of outlier
      detection algorithms performed on multivariate water-quality time series
      data (T, turbidity; C, conductivity; L, river level) from in situ sensors at Pioneer River, arranged in descending 
      order of OP values. See Sections 2.7-8 for performance metric codes and details.", booktabs = T) %>%
  kable_styling(latex_options = "scale_down")





## ---- runTime

test_data_PTCL <- pioneer_full[, c(
  "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
)] %>% drop_na()

test_data_PTCL <- test_data_PTCL[!is.infinite(rowSums(test_data_PTCL[
  ,
  c(
    "neg_der_log_bound_Turbidity", "pos_der_log_bound_Conductivity",
    "neg_der_log_bound_Level"
  ),
])), ]

test_data_PTC <- pioneer_full[, c(
  "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity"
)] %>% drop_na()

test_data_PTC <- test_data_PTC[!is.infinite(rowSums(test_data_PTC[
  ,
  c(
    "neg_der_log_bound_Turbidity",
    "pos_der_log_bound_Conductivity"
  ),
])), ]

test_data_STCL <- sandy_full[, c(
  "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
)] %>% drop_na()

test_data_STC <- sandy_full[, c(
  "neg_der_log_bound_Turbidity",
  "pos_der_log_bound_Conductivity"
)] %>%
  drop_na()

data <- sandy_full[, 8:9] %>% drop_na()

t1 <- runtime(test_data_STC)
t2 <- runtime(test_data_STCL)
t3 <- runtime(test_data_PTC)
t4 <- runtime(test_data_PTC)
t <- cbind(t1, t2, t3, t4)
kable(t)
