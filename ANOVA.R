setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyr)
library(agricolae)
library(dplyr)


outputBCTnBK <- data.frame(Analysis=character(),
                     PValue=double(),
                     Significant=character(),
                     stringsAsFactors=FALSE)

###################### BCT Pathway Scores Data ANOVA ########################################
BCT_file <- read.csv("BCT_ANOVA.csv")
BCT_data <- t(BCT_file)
BCT_data <- as.data.frame(BCT_data)

BCT_data$Bacteroidia <- BCT_data$V1
BCT_data$Negativicutes <- BCT_data$V2
BCT_data$Clostridia <- BCT_data$V3
BCT_data$Bacilli <- BCT_data$V4
BCT_data$Fusobacteriia <- BCT_data$V5

BCT_data[, c(1:5)] <- NULL
BCT_data <- BCT_data[-c(1), ]

columns_to_select <- c(1:5)
BCT_rundata <- BCT_data %>% pivot_longer(cols=all_of(columns_to_select), names_to = 'treatment', values_to = 'value', values_drop_na = TRUE)

BCT_rundata$value <- as.numeric(BCT_rundata$value)

BCT_aov <- aov(value ~ treatment, BCT_rundata)
BCT_anova <- anova(BCT_aov)

BCT_LSD <- LSD.test(BCT_aov, "treatment", group = FALSE)

(BCT_LSD)

save(BCT_LSD, file ="BCT_LSD.RData")

outputBCTnBK[nrow(outputBCTnBK) + 1,] <- c("BCT", BCT_anova$'Pr(>F)')

###################### BK Pathway Scores Data ANOVA ########################################
BK_file <- read.csv("BK_ANOVA.csv")
BK_data <- t(BK_file)
BK_data <- as.data.frame(BK_data)

BK_data$Bacteroidia <- BK_data$V1
BK_data$Negativicutes <- BK_data$V2
BK_data$Clostridia <- BK_data$V3
BK_data$Bacilli <- BK_data$V4
BK_data$Fusobacteriia <- BK_data$V5

BK_data[, c(1:5)] <- NULL
BK_data <- BK_data[-c(1), ]

columns_to_select <- c(1:5)
BK_rundata <- BK_data %>% pivot_longer(cols=all_of(columns_to_select), names_to = 'treatment', values_to = 'value', values_drop_na = TRUE)

BK_rundata$value <- as.numeric(BK_rundata$value)

BK_aov <- aov(value ~ treatment, BK_rundata)
BK_anova <- anova(BK_aov)



BK_LSD <- LSD.test(BK_aov, "treatment", group = FALSE)

(BK_LSD)

save(BK_LSD, file ="BK_LSD.RData")

outputBCTnBK[nrow(outputBCTnBK) + 1,] <- c("BK", BK_anova$'Pr(>F)')

print(outputBCTnBK)

write.csv(outputBCTnBK, "ANOVA_AnalysesBCTnBK.csv")

