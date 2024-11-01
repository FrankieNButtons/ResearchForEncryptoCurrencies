# 收益率序列的初步处理与准备
library("readxl");
prices <- read_excel("./data/daily_prices.xlsx");
date <- prices$date;

cols_to_remove <- c();
for (col_name in colnames(prices)) {
    if (startsWith(col_name, "date")) {
        cols_to_remove <- c(cols_to_remove, col_name);
    }
}
cleaned_prices <- prices[ !names(prices) %in% cols_to_remove];

final_prices <- as.data.frame(lapply(cleaned_prices, as.numeric));

returns <- as.data.frame(diff(log(as.matrix(final_prices))), stringsAsFactors = FALSE);

# 模块测试
pmiine_test <- mine(returns); # 可用
pmic_test <- pmine_test$MIC;
mic_trial <- mine(returns$BTC, returns$AMP)$MIC;
r_sliced <- getSlicedData(returns, 180, 7);
r_sliced[1];



# 用于按一定的窗口和步长切分同期数据并计算相应窗口期中某种参数值的某种网络
# 数据拆分函数getSlicedData
getSlicedData <- function(df_prices, winwidth, step) {
  
  # **description** 
  # 该函数根据指定的窗口宽度和步长对给定的数据框进行切分，生成多个较小的数据框。
  
  # **params** 
  # df_prices: 包含价格或其他时间序列数据的数据框。
  # winwidth: 用于切分数据框的窗口宽度（整数）。
  # step: 用于移动窗口的步长（整数）。
  
  # **returns** 
  # 返回一个切分后的数据框列表。
  
  total_rows <- nrow(df_prices);
  
  sliced_data <- list();
  
  # 使用for循环切分数据
  for (start in seq(1, total_rows - winwidth + 1, by = step)) {
    # 根据起始索引和窗口宽度切分数据
    sliced_df <- df_prices[start:(start + winwidth - 1), ];
    # 将切分的数据框添加到列表中
    sliced_data <- append(sliced_data, list(sliced_df));
  }
  
  # 返回切分后的数据框列表
  return(sliced_data);
}


# 核心函数getNet
getNet <- function(df_prices, winwidth, step, net, cfg, num, rpt) {
  
  # **description** 
  # 该函数根据给定的网络类型（如 MIC 或 TE）和窗口参数，计算网络矩阵及其 p 值矩阵。
  
  # **params** 
  # df_prices: 包含价格或其他时间序列数据的数据框。
  # winwidth: 用于切分数据框的窗口宽度（整数）。
  # step: 用于移动窗口的步长（整数）。
  # net: 网络的类型，可以为 "MIC" 或 "TE"。
  # cfg: 配置参数（未使用，保留以备扩展）。
  # num: 整数，随机置换检验中每次置换的点位数。
  # rpt: 整数，随机置换检验的重复次数。
  
  # **returns** 
  # 返回一个包含两个列表的列表：网络矩阵列表和 p 值矩阵列表。
  
  # 截取窗口数据
  sliced_dfs <- getSlicedData(df_prices, winwidth, step);
  
  # 初始化网络矩阵列表和 p 值矩阵列表
  net_matrices <- list();
  p_value_matrices <- list();
  total_net_num <- length(sliced_dfs);
  
  # 判断网络类型
  if (is.na(net) || net == "MIC") {
    library(minerva);
    # 对每个切片数据进行 MIC 计算
    for (i in seq_along(sliced_dfs)) {
      start_time <- Sys.time();
      
      sliced_df <- sliced_dfs[[i]];
      mic_matrix <- matrix(0, nrow = ncol(sliced_df), ncol = ncol(sliced_df));  # 初始化 MIC 矩阵
      p_matrix <- matrix(1, nrow = ncol(sliced_df), ncol = ncol(sliced_df));  # 初始化 p 值矩阵
      
      # 对于每一对不同的变量，计算 MIC 和 p 值
      for (j in 1:(ncol(sliced_df) - 1)) {
        for (k in (j + 1):ncol(sliced_df)) {
          # 计算 MIC
          mic_result <- mine(sliced_df[, j], sliced_df[, k]);
          mic_value <- mic_result$MIC;
          mic_matrix[j, k] <- mic_value;
          mic_matrix[k, j] <- mic_value;
          
          # 进行随机置换检验计算 p 值
          perm_count <- 0;
          for (s in 1:rpt) {
            # 对变量 k 进行部分随机置换
            permuted_y <- sliced_df[, k];
            indices_to_permute <- sample(1:length(permuted_y), num);  # 随机选择 num 个点位
            permuted_y[indices_to_permute] <- sample(permuted_y[indices_to_permute]);  # 对选中的点位进行置换
            
            # 计算置换后的 MIC 值
            perm_mic_result <- mine(sliced_df[, j], permuted_y);
            perm_mic_value <- perm_mic_result$MIC;
            
            if (perm_mic_value >= mic_value) {
              perm_count <- perm_count + 1;
            }
          }
          p_value <- (perm_count + 1) / (rpt + 1);  # 防止 p 值为 0
          p_matrix[j, k] <- p_value;
          p_matrix[k, j] <- p_value;
        }
      }
      
      # 将计算结果添加到列表
      net_matrices <- append(net_matrices, list(mic_matrix));
      p_value_matrices <- append(p_value_matrices, list(p_matrix));
      
      end_time <- Sys.time();
      time_taken <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2);
      cat(sprintf("Calculated MIC Net %d/%d and p_value cost %.2f seconds\n", i, total_net_num, time_taken));
    }
  } else if (net == "TE") {
    library(RTransferEntropy);
    # 对每个切片数据进行转移熵计算
    for (i in seq_along(sliced_dfs)) {
      start_time <- Sys.time();
      
      sliced_df <- sliced_dfs[[i]];
      te_matrix <- matrix(0, nrow = ncol(sliced_df), ncol = ncol(sliced_df));  # 初始化 TE 矩阵
      p_matrix <- matrix(1, nrow = ncol(sliced_df), ncol = ncol(sliced_df));  # 初始化 p 值矩阵
      
      for (j in 1:ncol(sliced_df)) {
        for (k in 1:ncol(sliced_df)) {
          # 计算转移熵和 p 值
          te_result <- transfer_entropy(sliced_df[, j], sliced_df[, k], quiet = TRUE);
          te_matrix[j, k] <- te_result$coef;
          p_matrix[j, k] <- te_result$p_value;
        }
      }
      
      # 将计算结果添加到列表
      net_matrices <- append(net_matrices, list(te_matrix));
      p_value_matrices <- append(p_value_matrices, list(p_matrix));
      
      end_time <- Sys.time();
      time_taken <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2);
      cat(sprintf("Calculated TE Net %d/%d and p_value cost %.2f seconds\n", i, total_net_num, time_taken));
    }
  }
  
  # 返回网络矩阵和 p 值矩阵
  return(list(net_matrices, p_value_matrices));
}

## 实操计算
rst_new <- getNet(returns, 180, 7, "MIC", null, 5, 10);
