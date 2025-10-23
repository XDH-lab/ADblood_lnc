##AIBL

p_values_AIBL <- data.frame(
  p_value = numeric(ncol(beta_AIBL)),
  direction = character(ncol(beta_AIBL)),  # 新增方向性列
  stringsAsFactors = FALSE
)
for (i in 1:ncol(beta_AIBL)) {
  if (i %% 100 == 0) {
    cat(i, "\n")
  }
  beta <- cbind(beta_AIBL[i], clin_AIBL)
  beta$clin <- factor(beta$clin, levels = c("CN","AD"), labels = c("0", "1"))
  beta$sex <- factor(beta$sex, levels = c("Female","Male"), labels = c("0", "1"))
  beta$age <- as.numeric(beta$age)
  colnames(beta)[1] <- "beta"
  
  model <- glm(
    clin ~ beta + age + sex + Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + Neu + NK + Treg,
    data = beta,
    family = binomial
  )
  
  # 提取P值
  p_values_AIBL$p_value[i] <- summary(model)$coefficients["beta", "Pr(>|z|)"]
  
  # 提取回归系数（beta的估计值）并判断方向性
  coef_beta <- summary(model)$coefficients["beta", "Estimate"]
  if (coef_beta > 0) {
    p_values_AIBL$direction[i] <- "++"  # AD中上调
  } else if (coef_beta < 0) {
    p_values_AIBL$direction[i] <- "--"  # AD中下调
  } else {
    p_values_AIBL$direction[i] <- "NA"  # 无变化（理论上极少见）
  }
  
  rownames(p_values_AIBL)[i] <- colnames(beta_AIBL)[i]
}


##ADNI
p_values_ADNI <- data.frame(
  p_value = numeric(ncol(beta_ADNI)),
  direction = character(ncol(beta_ADNI)),  # 新增方向性列
  stringsAsFactors = FALSE
)
for (i in 1:ncol(beta_ADNI)) {
  if (i %% 100 == 0) {
    cat(i, "\n")
  }
  beta <- cbind(beta_ADNI[i], clin_ADNI)
  beta$clin <- factor(beta$clin, levels = c("CN","AD"), labels = c("0", "1"))
  beta$sex <- factor(beta$sex, levels = c("Female","Male"), labels = c("0", "1"))
  beta$age <- as.numeric(beta$age)
  colnames(beta)[1] <- "beta"
  model <- glm(
    clin ~ beta + age + sex + Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + Neu + NK + Treg,
    data = beta,
    family = binomial
  )
  # 提取P值
  p_values_ADNI$p_value[i] <- summary(model)$coefficients["beta", "Pr(>|z|)"]
  # 提取回归系数（beta的估计值）并判断方向性
  coef_beta <- summary(model)$coefficients["beta", "Estimate"]
  if (coef_beta > 0) {
    p_values_ADNI$direction[i] <- "++"  # AD中上调
    } else if (coef_beta < 0) {
    p_values_ADNI$direction[i] <- "--"  # AD中下调
    } else {
    p_values_ADNI$direction[i] <- "NA"  # 无变化（理论上极少见）
    }

    rownames(p_values_ADNI)[i] <- colnames(beta_ADNI)[i]
}

##AddNeuroMed
p_values_AddNeuroMed <- data.frame(
  p_value = numeric(ncol(beta_AddNeuroMed)),
  direction = character(ncol(beta_AddNeuroMed)),  # 新增方向性列
  stringsAsFactors = FALSE
)
for (i in 1:ncol(beta_AddNeuroMed)) {
  if (i %% 100 == 0) {
    cat(i, "\n")
  }
  beta <- cbind(beta_AddNeuroMed[i], clin_AddNeuroMed)
  beta$clin <- factor(beta$clin, levels = c("CN","AD"), labels = c("0", "1"))
  beta$sex <- factor(beta$sex, levels = c("Female","Male"), labels = c("0", "1"))
  beta$age <- as.numeric(beta$age)
  colnames(beta)[1] <- "beta"
  model <- glm(
    clin ~ beta + age + sex + Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + Neu + NK + Treg,
    data = beta,
    family = binomial
  )
  # 提取P值
  p_values_AddNeuroMed$p_value[i] <- summary(model)$coefficients["beta", "Pr(>|z|)"]
  # 提取回归系数（beta的估计值）并判断方向性
  coef_beta <- summary(model)$coefficients["beta", "Estimate"]
  if (coef_beta > 0) {
    p_values_AddNeuroMed$direction[i] <- "++"  # AD中上调
    } else if (coef_beta < 0) {
    p_values_AddNeuroMed$direction[i] <- "--"  # AD中下调
    } else {
    p_values_AddNeuroMed$direction[i] <- "NA"  # 无变化（理论上极少见）
    }
    rownames(p_values_AddNeuroMed)[i] <- colnames(beta_AddNeuroMed)[i]
}







##AIBL显著差异lncRNA
sum(p_values_AIBL$p_value < 0.05)
significant_lnc <- rownames(p_values_AIBL)[which(p_values_AIBL$p_value<=0.05)]

##AIBL显著lncRNA在ADNI、AddNeuroMed中方向相同的lncRNA
direction_AIBL <- p_values_AIBL[significant_lnc, "direction"]
direction_ADNI <- p_values_ADNI[significant_lnc, "direction"]
direction_AddNeuroMed <- p_values_AddNeuroMed[significant_lnc, "direction"]
consistent_direction <- (direction_AIBL == direction_ADNI) & (direction_AIBL == direction_AddNeuroMed)
consistent_lnc <- significant_lnc[consistent_direction]


consistent_lnc_AIBL <- beta_AIBL[,which(colnames(beta_AIBL)%in%consistent_lnc)]
consistent_lnc_AIBL <- cbind(clin_AIBL$clin,consistent_lnc_AIBL)
colnames(consistent_lnc_AIBL)[1] <- "clin"

saveRDS(consistent_lnc_AIBL,"/ADblood_lnc/datasets/consistent_lnc_AIBL.RDS")