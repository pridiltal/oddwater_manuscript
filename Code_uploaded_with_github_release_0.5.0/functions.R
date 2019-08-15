plot_original <- function(data, y, colour = out_label_Tur, y_label, title = y_label) {
  p <- ggplot(data, aes_string(x = "Timestamp", y = y)) +
    #geom_line() +
    geom_point(data = data, aes_string(colour = colour, size = colour, shape = colour)) +
    scale_colour_manual(values = c(
      "typical" = "black", "outlier" = "#d95f02",
      "neighbour" = "#1b9e77"
    )) +
    scale_size_manual(values = c("outlier" = 4, "typical" = 1, "neighbour" = 2)) +
    ylab(title) +
    scale_shape_manual(values = c("outlier" = 17, "typical" = 19, "neighbour" = 15)) +
    ylab(title) +
    theme(legend.position = "none", axis.title.x = element_blank(), legend.title = element_blank())
  return(p)
}



plot_trans <- function(x, y, data) {
  
  unitize <- function(z) {
    zrange <- range(z, na.rm = TRUE)
    if (!(dif <- diff(zrange))) return(rep(0,length(z)))
    (z - zrange[1])/dif
  }
  
  data[,x]<- unitize(data[,x])
  data[,y]<- unitize(data[,y])
  out_data_n <- data %>% filter(label_all %in% c("neighbour"))
  out_data_o <- data %>% filter(label_all %in% c("outlier"))
  p <- ggplot(data, aes_string(x, y, colour = "label_all",
                                     shape = "label_all", size = "label_all", alpha = "label_all")) +
    geom_point() +
    scale_colour_manual(
      name = "label_all", values =
        c(
          "outlier" = "#d95f02", "typical" = "black",
          "neighbour" = "#1b9e77"
        )
    ) +
    scale_shape_manual(
      name = "label_all", values =
        c("outlier" = 17, "typical" = 19, "neighbour" = 15)
    ) +
    scale_size_manual(values = c("outlier" = 2, "typical" = 2, "neighbour" = 2)) +
    scale_alpha_manual(values = c("outlier" = 1, "typical" = 1, "neighbour" = 0.8)) +
    geom_point(data = out_data_n, aes_string(x,y))+
    geom_point(data = out_data_o, aes_string(x,y))+
    theme(
      aspect.ratio = 1, legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),  legend.title = element_blank()
    )
  return(p)
}

## THis is my level order: Levels: FN FP TN TP

performance_plot <- function(data, test_data, method, flag = c("flag_TC", "flag_TCL"),
                             title = method, label_type = c("label_all"), variable) {
  results <- apply_ADmethod(data, test_data, method = method, flag)
  colour_data <- results$full_data

  code <- function(x, y) {
    if (x == "outlier" & y == "outlier") {
      return("TP")
    }
    if (x == "outlier" & y == "typical") {
      return("FN")
    }
    if (x == "typical" & y == "outlier") {
      return("FP")
    }
    if (x == "typical" & y == "typical") {
      return("TN")
    }
  }
  colour_data$per <- mapply(
    code, colour_data[[label_type]],
    colour_data[[flag]]
  )
  colour_data$per <- as.factor(colour_data$per)

  # reorder the levels
  colour_data$per <- ordered(colour_data$per, levels = c("TN", "FN", "FP", "TP"))
  
  colour_data[, variable]<- sqrt(colour_data[, variable])
  
  p <- ggplot(colour_data, aes_string(x = "Timestamp", y = variable)) +
    geom_point(data = colour_data, aes_string(
      fill = "per", alpha = "per",
      size = "per", shape = "per",
      colour = "per"
    )) +
    scale_fill_manual(values = c(
      "TN" = "black", "TP" = "#fdcc8a",
      "FN" = "#fc8d59", "FP" = "#d7301f"
    )) +
    scale_colour_manual(values = c(
      "TN" = "black", "TP" = "black",
      "FP" = "black", "FN" = "black"
    )) +
    scale_alpha_manual(guide = "none", values = list(
      "TN" = 0.3, "TP" = 1,
      "FP" = 1, "FN" = 1
    )) +
    scale_size_manual(values = c("TN" = 3, "FN" = 5, "FP"= 5,"TP"=  5)) +
    scale_shape_manual(values = c("TN"= 20, "FN" =23, "FP" =  22, "TP"= 24)) +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_blank(), 
      text = element_text(size = 14)
    ) +
    ylab(title)
  
  tp_data <- subset(colour_data, per == "TP")
  
  p <- p + 
    geom_point(data = tp_data, aes_string( x = "Timestamp", y = variable),
                                           fill = "#fdcc8a", alpha = 1,
                                           size = 5, shape = 24,
                                           colour = "black"
    )
  
 
  return(p)
}

apply_ADmethod <- function(full_data, test_data, method =
                             c(
                               "NN-HD", "KNN-AGG", "KNN-SUM",
                               "INFLO", "COF", "LDOF", "LOF", "RKOF"
                             ),
                           flag = c("flag_TC", "flag_TCL"),
                           adjustment = TRUE) {
  outlier_index <-  NULL
  if (method == "NN-HD") {
    out_score <- oddwater::NN_HD(test_data[, 2:ncol(test_data)])
    t <- find_theshold(out_score)
  }
  else if (method == "KNN-AGG") {
    out_score <- DDoutlier::KNN_AGG(test_data[, 2:ncol(test_data)],
      k_min = 1,
      k_max = 10
    )
    t <- find_theshold(out_score)
  }
  else if (method == "KNN-SUM") {
    out_score <- DDoutlier::KNN_SUM(test_data[, 2:ncol(test_data)], k = 10)
    t <- find_theshold(out_score)
  }
  else if (method == "INFLO") {
    out_score <- DDoutlier::INFLO(test_data[, 2:ncol(test_data)], k = 10)
    t <- find_theshold(out_score)
  }
  else if (method == "COF") {
    out_score <- DDoutlier::COF(test_data[, 2:ncol(test_data)], k = 10)
    t <- find_theshold(out_score)
  }
  else if (method == "LDOF") {
    out_score <- DDoutlier::LDOF(test_data[, 2:ncol(test_data)], k = 10)
    t <- find_theshold(out_score)
  }
  else if (method == "LOF") {
    out_score <- DDoutlier::LOF(test_data[, 2:ncol(test_data)], k = 10)
    t <- find_theshold(out_score)
  }
  else if (method == "RKOF") {
    out_score <- DDoutlier::RKOF(test_data[, 2:ncol(test_data)],
      k = 10,
      C = 1, alpha = 1, sigma2 = 1
    )
    t <- find_theshold(out_score)
  }
  out_time <- test_data[t, 1]
  index <- which(full_data$Timestamp %in% out_time$Timestamp)
  

  if (!is_empty(index) & adjustment) {
    data_var <- as.data.frame(full_data[, c(
      "neg_der_log_bound_Turbidity",
      "pos_der_log_bound_Conductivity", "neg_der_log_bound_Level"
    )])
    data1 <- drop_na(data_var)
    median_data <- apply(data1, 2, median)
    diff <- abs(data_var[index, ] - t(replicate(length(index), median_data)))
    anom_var <- apply(diff, 1, which.max)

    var <- ifelse(anom_var == 1, "Turbidity", ifelse(anom_var == 2, "Conductivity", "log_Level"))

    for (i in 1:length(index))
    {
      t <- index[i]
      index_pre <- t - 1
      index_next <- t + 1
      val <- as.data.frame(full_data[c(index_pre, t, index_next), var[i]])
      diff <- val - colMeans(val)
      if (var[i] %in% c("Turbidity", "Level")) {
        if (all(sign(diff) == c(1, -1, 1))) {
          index[i] <- t
        }
        else if (all(sign(diff) == c(1, -1, -1))) {
          index[i] <- index_pre
        }
      }
      if (var[i] == "Conductivity") {
        if (all(sign(diff) == c(-1, 1, -1))) {
          index[i] <- t
        }
        else if (all(sign(diff) == c(-1, 1, 1))) {
          index[i] <- index_pre
        }
      }
    }
    
    full_data <- as.data.frame(full_data)
    outlier_index <- cbind( full_data[index,1:4], var)
    
  }



  full_data[index, flag] <- "outlier"
  # print("full outlier set")
  # print(knitr::kable( full_data[which(full_data$flag_TCL == "outlier"),
  #                              c("Timestamp", "out_type_Level", "out_type_Cond",
  #                                "out_type_Tur", "label_all", "flag_TC" )]))
  out <- calc_performance_metrics(
    y_truth = full_data$label_all,
    y_output = full_data[[flag]],
    pos_label = "outlier",
    neg_label = "typical",
    print_out = FALSE
  )

  test_data$out_score <- out_score
  score <- test_data[, c("Timestamp", "out_score")]

  full_data <- left_join(full_data, score, by = "Timestamp")

  measures <- t(round(unlist(out), 4))
  return(list(
    "t" = t,
    "measures" = measures, "full_data" = full_data,
    "performance_M" = t(round(unlist(out), 4)), 
    "outlier_index" = outlier_index
  ))
}

find_theshold <- function(outlier_score, alpha = 0.05) {
  n <- length(outlier_score)
  ord <- order(outlier_score)
  gaps <- c(0, diff(outlier_score[ord]))
  n4 <- max(min(50, floor(n / 4)), 2)
  J <- 2:n4
  start <- max(floor(n / 2), 1) + 1
  ghat <- numeric(n)
  for (i in start:n) ghat[i] <- sum((J / (n4 - 1)) * gaps[i - J + 1 ]) # check i - j +1
  # J <- 1:n4
  #  start <- max(floor(n/2), 1) + 1
  #  ghat <- numeric(n)
  #  for (i in start:n) ghat[i] <- sum((J/(n4)) * gaps[i - J+1 ]) # check i - j +1
  logAlpha <- log(1 / alpha)
  bound <- Inf

  for (i in start:n) {
    if (gaps[i] > logAlpha * ghat[i]) {
      bound <- outlier_score[ord][i - 1]
      break
    }
  }
  ex <- which(outlier_score >= bound)

  return(ex)
}


runtime <- function(data) {
  mu_t <- max_t <- min_t <- NULL

  for (i in 1:8)
  {
    if (i == 1) {
      t <- microbenchmark::microbenchmark(dbscan::kNN(as.matrix(test_data[, 2:ncol(test_data)]), 1)$dist)
    }
    if (i == 2) {
      t <- microbenchmark::microbenchmark(KNN_AGG(data, k_min = 1, k_max = 10))
    }
    if (i == 3) {
      t <- microbenchmark::microbenchmark(KNN_SUM(data, k = 10))
    }
    if (i == 4) {
      t <- microbenchmark::microbenchmark(INFLO(data, k = 10))
    }
    if (i == 5) {
      t <- microbenchmark::microbenchmark(COF(out_data[, 2:v], k=10))
    }
    if (i == 6) {
      t <- microbenchmark::microbenchmark( LDOF(out_data[, 2:v], k= 10))
    }
    if (i == 7) {
      t <- microbenchmark::microbenchmark(LOF(data, k = 10))
    }
    if (i == 8) {
      t <- microbenchmark::microbenchmark(RKOF(data, k = 10, C = 1, alpha = 1, sigma2 = 1))
    }
    min_t <- c(min_t, min(t$time) / 10^6)
    mu_t <- c(mu_t, mean(t$time) / 10^6)
    max_t <- c(max_t, max(t$time) / 10^6)
  }
  nvar <- rep(ncol(data), 8)
  method <- c("NN-HD", "KNN-AGG", "KNN-SUM", "INFLO", "COF", "LDOF", "LOF", "RKOF")
  output <- data.frame(nvar, method, min_t, mu_t, max_t)
  print(output)
  return(output)
  
}

plot_output <- function(data, variable, label_type, title)
{
  
  code <- function(x, y) {
    if (x == "outlier" & y == "outlier") {
      return("TP")
    }
    if (x == "outlier" & y == "typical") {
      return("FN")
    }
    if (x == "typical" & y == "outlier") {
      return("FP")
    }
    if (x == "typical" & y == "typical") {
      return("TN")
    }
  }
  
  out_data$per <- mapply(
    code, out_data[[label_type]],
    out_data[["col"]]
  )
  
  out_data$per <- as.factor(out_data$per)
  out_data$per <- ordered(out_data$per, levels = c("TN", "FN", "FP", "TP"))
  
  p <- ggplot(out_data, aes_string(x = "Timestamp", y = variable)) +
    geom_point(data = out_data, aes_string(
      fill = "per", alpha = "per",
      size = "per", shape = "per",
      colour = "per"
    )) +
    scale_fill_manual(values = c(
      "TN" = "black", "TP" = "#fdcc8a",
      "FN" = "#fc8d59", "FP" = "#d7301f"
    )) +
    scale_colour_manual(values = c(
      "TN" = "black", "TP" = "black",
      "FP" = "black", "FN" = "black"
    )) +
    scale_alpha_manual(guide = "none", values = list(
      "TN" = 0.5, "TP" = 1,
      "FP" = 1, "FN" = 1
    )) +
    scale_size_manual(values = c("TN" = 3, "FN" = 5,  "FP" =5, "TP" =5)) +
    scale_shape_manual(values = c("TN" = 20, "FN" = 23, "FP" = 22,  "TP" =24)) +
    theme(
      legend.position = "bottom", axis.title.x = element_blank(),
      axis.text.x = element_blank(), axis.text.y = element_blank(),
      text = element_text(size = 14),  legend.title = element_blank()
    ) +
    ylab(title)
  return(p)
}
