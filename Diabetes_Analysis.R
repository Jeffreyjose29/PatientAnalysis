library(rJava)
library(RJDBC)
library(plyr)
library(dplyr)
library(PerformanceAnalytics)
library(ggplot2)
require(ggplot2)
library(hrbrthemes)
library(fmsb)
library(mice)    #Potentially use mice to impute missing data within the labs table
library(VIM)
library(leaps)
#library(MASS)   #stops the select() from dplyr to stop working
library(lindia)   #Diagnostic Plots Won't work without these packages
library(gridExtra)
library(grid)
library(caret)
library(boot)
library(Rtsne)
library(cluster)
library(Amelia)
library(psycho)
options(warn = -1)

#######################################



THIS SECTION COVERS THE CODE TO EXTRACT THE DATA. DUE TO CONFIDENTIALITY, THIS CANNOT BE SHOWN.




######################################

setwd("C:/Users/jeffr/OneDrive/Desktop/Honours Project")
############ I TRIED TO IMPUTE ACR AND EGFR USING THE MICE PACKAGE (THESE ARE THE ONLY VALUES WHERE IT IS POTENTIALLY USEFUL TO HAVE ACTUAL DATA)
# trying_imputation <- tbl_labs
# md.pattern(trying_imputation)
# aggr_plot <- aggr(trying_imputation, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(trying_imputation), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
# marginplot(trying_imputation[c(2,3)])
# 
# # Imputation
# tempData <- mice(trying_imputation, m =5, meth = 'pmm', seed = 500)
# summary(tempData)
# 
# tempData$imp$acr
# tempData$method
# completedData <- complete(tempData, 1) #Imputation for acr and efgr
# completedData$Imputed <- "Yes"
# 
# labs <- tbl_labs
# labs$Imputed <- "No"
# imputedLabs <- left_join(labs, completedData, by = "nhi")
# str(imputedLabs)
# 
# summary(imputedLabs$hba1c.x)
# summary(imputedLabs$hba1c.y)
# summary(imputedLabs$acr.x)
# summary(imputedLabs$acr.y)
# summary(imputedLabs$egfr.x)
# summary(imputedLabs$egfr.y)
# #The summary statistic between the imputed and actual data does not seem to be that different
# 
# #Trying to plot imputed data
# imputedLabs$indicator <- seq(from = 1, to = nrow(imputedLabs))
# completedData <- imputedLabs %>% select(nhi, hba1c.y, acr.y, egfr.y, Imputed.y, indicator)
# labs <- imputedLabs %>% select(-c(hba1c.y, acr.y, egfr.y, Imputed.y))
# colnames(completedData)[2] <- "hba1c.x"
# colnames(completedData)[3] <- "acr.x"
# colnames(completedData)[4] <- "egfr.x"
# colnames(completedData)[5] <- "Imputed.x"
# 
# graph_imputed <- rbind(labs, completedData)
# col_plot <- c("#9814F9", "#F9CF14")
# graph_imputed %>% filter(indicator < 200) %>% ggplot(aes(x = indicator, y = acr.x, colour = as.factor(Imputed.x), group = Imputed.x)) + geom_line(position = "dodge") +
#   facet_grid(~Imputed.x) + theme_ft_rc() + xlab("Patient Indicator Value") + ylab("Acr") + scale_color_manual(values = col_plot, name = "Imputed Values:") +
#   ggtitle("Comparison Of Imputed And Non-Imputed Acr Readings For The First 200 Patients Within Our Dataset")
# graph_imputed %>% filter(indicator < 200) %>% ggplot(aes(x = indicator, y = egfr.x, colour = as.factor(Imputed.x), group = Imputed.x)) + geom_line(position = "dodge") +
#   facet_grid(~Imputed.x) + theme_ft_rc() + xlab("Patient Indicator Value") + ylab("eGFR") + scale_color_manual(values = col_plot, name = "Imputed Values:") +
#   ggtitle("Comparison Of Imputed And Non-Imputed eGFR Readings For The First 200 Patients Within Our Dataset")

trying_imputation <- tbl_labs
trying_imputation <- left_join(trying_imputation, tbl_measurements, by = "nhi")
trying_imputation <- trying_imputation %>% filter(between(bmi, 15.0, 40.0))

joined <- merge(tbl_patients, trying_imputation, by = "nhi")
joined <- merge(joined, tbl_events, by = "nhi")

joined$aucode[is.na(joined$aucode)] <- -1 #(-1 means cant geocode)
joined$gp[is.na(joined$gp)] <- 0
joined$nurse[is.na(joined$nurse)] <- 0
joined$other[is.na(joined$other)] <- 0
joined$ed[is.na(joined$ed)] <- 0
joined$ash[is.na(joined$ash)] <- 0
joined$portal[is.na(joined$portal)] <- 0
joined$op[is.na(joined$op)] <- 0
joined$beddays[is.na(joined$beddays)] <- 0
joined$ip[is.na(joined$ip)] <- 0

joined$systolic_pressure <- as.numeric(substr(joined$bp, 1, 3)) #Usually from 90 to 250 
joined$diastolic_pressure <- as.numeric(sapply(strsplit(joined$bp, "/"), "[", 2))

joined$aucode <- NULL; joined$bp <- NULL

tempData <- mice(joined, m =5, meth = 'pmm')
summary(tempData)
imputedData <- complete(tempData, 1)
completedData <- complete(tempData, 1) #Imputation for acr and efgr
completedData$Imputed <- "Yes"

labs <- tbl_labs
labs <- left_join(labs, tbl_measurements, by = "nhi")
labs <- labs %>% filter(between(bmi, 15.0, 40.0))
labs$systolic_pressure <- as.numeric(substr(labs$bp, 1, 3)) #Usually from 90 to 250 
labs$diastolic_pressure <- as.numeric(sapply(strsplit(labs$bp, "/"), "[", 2))
labs$bp <- NULL
labs$Imputed <- "No"
imputedLabs <- merge(labs, completedData, by = "nhi")
# x is not imputed
str(imputedLabs)

summary(imputedLabs$hba1c.x)
summary(imputedLabs$hba1c.y)
summary(imputedLabs$acr.x)
summary(imputedLabs$acr.y)
summary(imputedLabs$egfr.x)
summary(imputedLabs$egfr.y)

ks.test(imputedLabs$hba1c.x, imputedLabs$hba1c.y)
ks.test(imputedLabs$acr.x, imputedLabs$acr.y)
ks.test(imputedLabs$egfr.x, imputedLabs$egfr.y)

imputedLabs <- imputedLabs %>% select(nhi, hba1c.x, hba1c.y, acr.x, acr.y, egfr.x, egfr.y, systolic_pressure.x, systolic_pressure.y, diastolic_pressure.x, diastolic_pressure.y, Imputed.x, Imputed.y)
imputedLabs$indicator <- seq(from = 1, to = nrow(imputedLabs))
completedData <- imputedLabs %>% select(nhi, hba1c.y, acr.y, egfr.y, Imputed.y, indicator, systolic_pressure.y, diastolic_pressure.y)
labs <- imputedLabs %>% select(-c(hba1c.y, acr.y, egfr.y, Imputed.y, systolic_pressure.y, diastolic_pressure.y))

colnames(completedData)[2] <- "hba1c.x"
colnames(completedData)[3] <- "acr.x"
colnames(completedData)[4] <- "egfr.x"
colnames(completedData)[5] <- "Imputed.x"
colnames(completedData)[7] <- "systolic_pressure.x"
colnames(completedData)[8] <- "diastolic_pressure.x"

graph_imputed <- rbind(labs, completedData)

col_plot <- c("#9814F9", "#F9CF14")
graph_imputed %>% filter(indicator < 200) %>% ggplot(aes(x = indicator, y = acr.x, colour = as.factor(Imputed.x), group = Imputed.x)) + geom_line(position = "dodge") +
  facet_grid(~Imputed.x) + theme_ft_rc() + xlab("Patient Indicator Value") + ylab("Acr") + scale_color_manual(values = col_plot, name = "Imputed Values:") +
  ggtitle("Comparison Of Imputed And Non-Imputed Acr Readings For The First 200 Patients Within Our Dataset")
graph_imputed %>% filter(indicator < 200) %>% ggplot(aes(x = indicator, y = egfr.x, colour = as.factor(Imputed.x), group = Imputed.x)) + geom_line(position = "dodge") +
  facet_grid(~Imputed.x) + theme_ft_rc() + xlab("Patient Indicator Value") + ylab("eGFR") + scale_color_manual(values = col_plot, name = "Imputed Values:") +
  ggtitle("Comparison Of Imputed And Non-Imputed eGFR Readings For The First 200 Patients Within Our Dataset")

#########################################################################################
################## I ALSO TRIED IMPUTATION USING HOT.DECK ##############################
amelia.data <- tbl_labs
amelia.data <- joined; amelia.data$log_acr <- NULL; amelia.data$log_hba1c <- NULL
missmap(amelia.data) #Shows the percentage missing from our lab readings
amelia.data <- amelia.data %>% select(nhi, hba1c, acr, egfr, systolic_pressure, diastolic_pressure)
amelia.data.imputed <- amelia(amelia.data, m = 5, idvars = "nhi")
summary(amelia.data.imputed)
head(amelia.data.imputed$imputations$imp1)
write.amelia(obj = amelia.data.imputed, file.stem = "C:/Users/jeffr/OneDrive/Desktop/Honours Project/Amelia Imputations")

#Read back in every amelia file
temp_ = list.files(pattern="*.csv")
for (i in 1:length(temp_)) assign(temp_[i], read.csv(temp_[i]))

#Remove the column X from the amelia files 
`Amelia Imputations1.csv`$X <- NULL
`Amelia Imputations2.csv`$X <- NULL
`Amelia Imputations3.csv`$X <- NULL
`Amelia Imputations4.csv`$X <- NULL
`Amelia Imputations5.csv`$X <- NULL

`Amelia Imputations1.csv`$Imput_Num <- "1"
`Amelia Imputations2.csv`$Imput_Num <- "2"
`Amelia Imputations3.csv`$Imput_Num <- "3"
`Amelia Imputations4.csv`$Imput_Num <- "4"
`Amelia Imputations5.csv`$Imput_Num <- "5"

`Amelia Imputations1.csv`$indicator <- seq(from = 1, to = nrow(`Amelia Imputations1.csv`))
`Amelia Imputations2.csv`$indicator <- seq(from = 1, to = nrow(`Amelia Imputations2.csv`))
`Amelia Imputations3.csv`$indicator <- seq(from = 1, to = nrow(`Amelia Imputations3.csv`))
`Amelia Imputations4.csv`$indicator <- seq(from = 1, to = nrow(`Amelia Imputations4.csv`))
`Amelia Imputations5.csv`$indicator <- seq(from = 1, to = nrow(`Amelia Imputations5.csv`))


amelia.imputations <- rbind(`Amelia Imputations1.csv`, `Amelia Imputations2.csv`, `Amelia Imputations3.csv`, `Amelia Imputations4.csv`, `Amelia Imputations5.csv`)

amelia.imputations %>% filter(indicator <= 200) %>%
  ggplot(aes(x = indicator, y = acr, group = Imput_Num, colour = Imput_Num)) + geom_line(position = "dodge") + theme_ft_rc() + scale_colour_discrete(name = "Imputation Number: ") +
  ggtitle("Amelia Imputation On ACR For The First 200 Patients") + xlab("Patient Indicator") + ylab("Acr Reading")
amelia.imputations %>% filter(indicator <= 200) %>%
  ggplot(aes(x = indicator, y = egfr, group = Imput_Num, colour = Imput_Num)) + geom_line(position = "dodge") + theme_ft_rc() + scale_colour_discrete(name = "Imputation Number: ") +
  ggtitle("Amelia Imputation On eGFR For The First 200 Patients") + xlab("Patient Indicator") + ylab("eGFR Reading")
amelia.imputations %>% filter(indicator <= 200) %>%
  ggplot(aes(x = indicator, y = hba1c, group = Imput_Num, colour = Imput_Num)) + geom_line(position = "dodge") + theme_ft_rc() + scale_colour_discrete(name = "Imputation Number: ") +
  ggtitle("Amelia Imputation On HbA1c For The First 200 Patients") + xlab("Patient Indicator") + ylab("HbA1c Reading")

# Comment: Although amelia imputation is a sound accurate for of imputation, I will be using MICE imputated data for my regression
# analysis as this is gave me more accurate results.
##########################################################################################################
#Mutation on patient population
tbl_patient_ethnicity <- regression_data.csv %>% group_by(ethgroup) %>% summarise(count = n())
tbl_patient_ethnicity$percentage <- round((tbl_patient_ethnicity$count / 16308) * 100, digit = 2)
tbl_patient_eth_age <- regression_data.csv %>% group_by(ethgroup, agegroup) %>% summarise(count = n())

#Exploratory Plot

tbl_patients_test <- regression_data.csv %>% group_by(type) %>% summarise(count = n())
tbl_patients_test %>% ggplot(aes(x = type, y = count, group = type, fill = type)) + geom_col() + theme_ft_rc() + 
  ggtitle("Number Of Patients By Diabetes Group")


#The following two plots show the percentage of patients within our dataset subsetted by ethnicity and the number of patients in 
#each ethnicity subsetted by their age group.
mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
tbl_patient_ethnicity %>%
  ggplot(aes(x = ethgroup, y = count, fill = ethgroup, group = ethgroup, label = percentage)) + geom_bar(stat = "identity", colour = "white") + theme_ft_rc() + xlab("Ethnic Groups") + 
  ylab("Number Of Patients") + ggtitle("Number Of Patients By Ethnicity") + scale_fill_discrete(name = "Ethnic Groups:") + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label = paste0(percentage, "%")), position=position_dodge(width=0.9), vjust=-0.25, colour = "white")

ageByEthnicity <- function(ethnicity){
  tbl_patient_eth_age %>% filter(ethgroup == ethnicity) %>%
    ggplot(aes(x = 2, y = count, fill = agegroup)) + geom_bar(stat = "identity", color = "black", lwd = 0.6) + 
    coord_polar(theta = "y", start = 0) + scale_fill_manual(values = mycols) + theme_void() + xlim(0.5, 2.5) + ggtitle(paste0("Count Of Patients In Each Age Group By Ethnicity: ", ethnicity)) + 
    theme(plot.title = element_text(hjust = 0.5)) + scale_fill_discrete(name = "Age Group:")
}
ageByEthnicity("Other")

rm(tbl_patient_eth_age, tbl_patient_ethnicity)

#Transforming the hba1c and 
tbl_labs$log_hba1c <- log10(tbl_labs$hba1c)
tbl_labs$log_acr <- log10(tbl_labs$acr)

#Joining all tables together by nhi
joined_table <- tbl_patients %>%
  merge(tbl_events, by = "nhi") %>%
  merge(tbl_labs, by = "nhi") %>%
  merge(tbl_measurements, by = "nhi")

#Change structure of some columns to factor
change_col <- c("type", "agegroup", "ethgroup", "gender", "currentsmoker", "engaged", "huhc", "csc", "cancer", "respiratory",
                "heart", "smoker", "mental")
joined_table[change_col] <- lapply(joined_table[change_col], factor)

############################################################################################################################
joined_table$systolic_pressure <- as.numeric(substr(joined_table$bp, 1, 3)) #Usually from 90 to 250 
joined_table$diastolic_pressure <- as.numeric(sapply(strsplit(joined_table$bp, "/"), "[", 2))

#Descriptive Statistics On BMI Without The Filtering
summary(joined_table$bmi)

uncommon_bmi <- joined_table %>% filter(bmi < 14.9 | bmi >= 40.1) #patients that have a very high bmi reading (very outlier)
joined_table <- joined_table %>% filter(between(bmi, 15.0, 40)) #patient table count went from 21349 -> 19216

#Descriptive Statistics On BMI With The Filtering
summary(joined_table$bmi)
joined_table %>% ggplot(aes(x = bmi)) + geom_boxplot(colour = "#DB16F2", fill = "#6716F2") + theme_ft_rc() + ggtitle("Overall Distribution Of BMI")

#Spread of BMI Across Each Age Group, Gender And Ethnicity
bar_cols = c("#F9257B","#9209FB","#09A5FB","#00D519","#E6E900")
joined_table %>% ggplot(aes(x = bmi, fill = agegroup)) + geom_histogram(colour = "white") + facet_wrap(~agegroup, scale = "free") + theme_ft_rc() + ggtitle("Spread Of BMI Across The Different Age Group") +
  scale_fill_manual(values = bar_cols, name = "Age Group") + ylab("Count Of Patients") + xlab("BMI Reading")
joined_table %>% ggplot(aes(x = bmi, fill = gender)) + geom_histogram(colour = "white") + facet_wrap(~gender, scale = "free") + theme_ft_rc() + ggtitle("Spread Of BMI Across The Two Genders") +
  scale_fill_manual(values = bar_cols, name = "Gender:") + ylab("Count Of Patients") + xlab("BMI Reading")
joined_table %>% ggplot(aes(x = bmi, fill = ethgroup)) + geom_histogram(colour = "white") + facet_wrap(~ethgroup, scale = "free") + theme_ft_rc() + ggtitle("Spread Of BMI Across The Different Ethnicities") +
  scale_fill_manual(values = bar_cols, name = "Ethnicity:") + ylab("Count Of Patients") + xlab("BMI Reading")

# Simple plots to see any evident relatioship
type_and_hba1c <- joined_table %>% select(type, hba1c, log_hba1c)
plot_hba1c_distribution <- type_and_hba1c %>%
  group_by(type, hba1c, log_hba1c) %>%
  tally()
# The distribution of hba1c reading grouped by diabetes type
plot_hba1c_distribution %>% ggplot(aes(x = hba1c, y = n, colour = as.factor(type), group = as.factor(type))) + geom_line(position = "dodge") + 
  ggtitle("Distribution Of HbA1c Grouped By Diabetes Type") + xlab("HbA1c Readings") + ylab("Count Of Patients") +
  scale_color_discrete(name = "Diabetes Group: ") + theme_ft_rc()
# The distribution of log hba1c reading grouped by diabetes type
col <- c("#B300AB", "#DF005B", "#09CFF4")
plot_hba1c_distribution %>% ggplot(aes(x = log_hba1c, y = n, colour = as.factor(type), group = as.factor(type))) + geom_line(position = "dodge") + 
  ggtitle("Distribution Of Log HbA1c Grouped By Diabetes Type") + xlab("HbA1c Readings") + ylab("Count Of Patients") +
  scale_colour_manual(values = col, name = "Diabetes Group: ") + theme_ft_rc()

rm(plot_hba1c_distribution)

hba1c_median <- type_and_hba1c %>% group_by(type) %>% summarise(average_measure = median(hba1c, na.rm = TRUE))
hba1c_mean <- type_and_hba1c %>% group_by(type) %>% summarise(average_measure = mean(hba1c, na.rm = TRUE))
hba1c_median$average <- "median"
hba1c_mean$average <- "mean"
hba1c_average <- rbind(hba1c_mean, hba1c_median)

#See if there is any major difference between the two averages (median and mean) for the haemoglobin levels of patients within
#the different diabetic groups
hba1c_average %>% ggplot(aes(x = type, y = average_measure, group = average, fill = average)) + geom_col(position = "dodge", stat = "identity") +
  theme_ft_rc() + ggtitle("HbA1c Averages Grouped By Diabetes Type") + scale_fill_manual(name = "Average Used:", values = c("mean" = "red", "median" = "orange")) + labs(y = "Average HbA1c", x = "Diabetes Type") + 
  labs(caption = "DM is a non-defined type")

rm(hba1c_median, hba1c_mean)

# HBA1C by secondary care usage
joined_table %>% ggplot(aes(x = hba1c, y = , fill = as.factor(huhc), group = as.factor(huhc))) + geom_bar(colour = "black") + theme_ft_rc() +
  ggtitle("HbA1c For Patients With High User Health Cards vs. Those Patients Without High User Health Card") + 
  scale_fill_manual(values = mycols, name = "High User Health Card?") + ylab("Number Of Patients") + xlab("HbA1c Readings")

joined_table %>% ggplot(aes(x = log_hba1c, y = egfr, color = as.factor(agegroup))) + geom_point() + theme_ft_rc() + ggtitle("Highlighting Any Existing Relationship Between HbA1c & eGFR") +
  scale_colour_discrete(name = "Age Group:") + ylab("eGFR") + xlab("log HbA1c")
numerical_cols <- joined_table[,c(28,29,30,31,32)]
numerical_cols <- na.omit(numerical_cols)
#Correlation between all the numerical values
round(cor(numerical_cols, use = "complete.obs"), 2)         #"Try correlation after removing records where hba1c < 30" - Dr Ryan
rm(numerical_cols)

#Quick analysis of the smoking habits of maori patients
maori_population <- joined_table %>% filter(ethgroup == "Maori")
maori_smoke <- maori_population %>% group_by(smoker, agegroup, gender) %>% summarise(count = n())
maori_smoke$porportion <- maori_smoke$count / sum(maori_smoke$count) #Averaged over the whole population (is it better to average over total in the age group)

#Percentage within the maori population that smoke subsetted by their age group
maori_smoke %>% ggplot(aes(x = agegroup, y = porportion, fill = smoker, group = smoker)) + geom_col(position = "dodge") + 
  facet_grid(.~gender, scales = "free") + theme_ft_rc() + scale_fill_manual(name = "Smoker?", values = c("No" = "Green", "Yes" = "Yellow")) + ggtitle("Porportion Of Smokers vs Non-Smokers Among Genders For Maori Patients") +
  labs(y = "Porportion Of Maori Patients", x = "Age Group") 

smoking_habit <- joined_table %>% group_by(smoker, ethgroup, agegroup) %>% summarise(count = n())
smoking_habit %>% ggplot(aes(x = agegroup, y = count, fill = agegroup, group = agegroup)) + geom_col() + facet_wrap(~ethgroup, scale = "free") + theme_ft_rc() +
  scale_fill_discrete(name = "Age Group:") + xlab("Age Group") + ylab("Number Of Patients") + ggtitle("Number Of Patients Who Smoke From Each Ethnic Group Grouped By Age Category")

#According to American Association of Kidney Patients (AAKP), studies have shown that smoking is harmful for the kidneys and can
#cause kidney disease to progress and increases the risk of proteinuria.
#Check for this using the smoking habits and acr reading
smoker_acr_egfr <- joined_table %>% select(smoker, acr)
smoker_acr_egfr_tally <- smoker_acr_egfr %>% group_by(smoker, acr) %>% tally()

#I've put a filter on this so far
cols = c("#E9FB12", "#12CDFB")
smoker_acr_egfr_tally %>% filter(acr < 100) %>% ggplot(aes(x = acr, y = n, colour = smoker, group = smoker)) + geom_line(position = "dodge") + theme_ft_rc() +
  ggtitle("ACR Levels For Patients Who Have Previously Smoked vs. Patients That Have Not Previously Smoked") + scale_color_manual(values = cols, name = "Previous Smoker:") +
  xlab("ACR Level") + ylab("Number Of Patients")

smokAcr_aov <- aov(acr ~ smoker, data = smoker_acr_egfr)
summary(smokAcr_aov) #p-value is a measure of the correlation (p < 0.05 means statistically significant and cannot reject null hypothesi)s
smokEgfr_aov <- aov(egfr ~ smoker, data = joined_table)
summary(smokEgfr_aov)
summary(aov(acr ~ egfr, data = joined_table))

#hba1c < 35 -might not have diabetes

#Correlation between high care needs and high user health cards
summary(aov(hcn ~ huhc, data = joined_table))
huhc.vs.hcn <- joined_table %>% select(hcn, huhc, agegroup)
huhc.vs.hcn_plot <- huhc.vs.hcn %>% group_by(hcn, huhc, agegroup) %>% tally()
huhc.vs.hcn_plot %>% ggplot(aes(x = hcn, y = n, fill = agegroup, group = agegroup)) + geom_col(position = "dodge", stat = "identity") + theme_ft_rc() +
  facet_wrap(~huhc, scale = "free") + xlab("HCN") + ylab("Number Of Patients") + ggtitle("Distribution Of High Care Needs For Patients Based On Whether They Are A High Case User Card Holder") +
  scale_fill_manual(values = bar_cols, name = "Patient Age:")

rm(huhc.vs.hcn, huhc.vs.hcn_plot)

#Risk factors for Waikato patients with diabetes not accessing secondary diabetes care
youth <- tbl_patients %>% filter(agegroup == "15-24") #dataset has 293 youths in total
youth_type_count <- youth %>% group_by(type) %>% tally()
youth_type_count %>% ggplot(aes(x = type, y = n, fill = type, group = type)) + geom_col() + theme_ft_rc() + 
  geom_text(aes(label = paste0(n, " Patients")), position=position_dodge(width=0.9), vjust=-0.25, colour = "white") +
  scale_fill_manual(values = bar_cols, name = "Diabetes Type") + xlab("Diabetes Type") + ylab("Number Of Patients") + 
  ggtitle("Youth Numbers By Diabetes Type") #High number of dms so I will try assign them to a diabetes type
rm(youth_type_count)
youth_labs <- left_join(youth, tbl_labs, by = "nhi")
not_recorded <- youth_labs[is.na(youth_labs$hba1c),]

recorded_cols <- c("#DAF7A6","#FFC300","#FF5733","#C70039")
t1_patients <- tbl_patients %>% filter(type == "T1") #dataset has 327 patients with T1 diabetes in total
t1_patients %>% ggplot(aes(x = agegroup, group = agegroup, fill = agegroup)) + geom_bar() + theme_ft_rc() + xlab("Age Group") + ylab("Number Of Patients") +
  scale_fill_discrete(name = "Age Group:") + ggtitle("Type 1 Patients By Age Group")
t1_patients %>% ggplot(aes(x = ethgroup, group = ethgroup, fill = ethgroup)) + geom_bar() + theme_ft_rc() + xlab("Ethnicity") + ylab("Number Of Patients") +
  scale_fill_discrete(name = "Ethnicity:") + ggtitle("Type 1 Patients By Ethnic Group")
t1_patients %>% ggplot(aes(x = agegroup, group = agegroup, fill = agegroup)) + geom_bar() + theme_ft_rc() + xlab("Age Group") + ylab("Number Of Patients") +
  scale_fill_manual(values = bar_cols, name = "Age Group:") + ggtitle("Age Of Type 1 Patients Among Each Ethnic Group") + facet_wrap(~ethgroup, scale = "free")
t1_patients %>% ggplot(aes(x = hcn, fill = agegroup, group = agegroup)) + geom_bar(position = "dodge") + theme_ft_rc() +
  scale_fill_manual(values = recorded_cols, name = "Age Group:") + xlab("HCN") + ylab("Number Of Patients") + ggtitle("HCN Levels For Patients With T1 Diabetes Subsetted By Age")

t1Patients_labs <- left_join(t1_patients, tbl_labs, by = "nhi")
not_recorded_t1Patients <- t1Patients_labs[is.na(t1Patients_labs$hba1c),]  
  
  
allPatients_labs_not_recorded <- joined_table[is.na(joined_table$hba1c),]  
  
#Renal (All patients with uacr > 70)
filtered_acr <- tbl_labs %>% filter(acr > 70)
filtered_patients <- tbl_patients %>% filter(type %in% c("T1", "T2"))
acr_patients_filtered <- left_join(filtered_acr, filtered_patients, by = "nhi")
# NOTE:
# Out Of all the patients that are diagnosed with T1 and T2 diabetes, only 483 have recorded acr values above 70 and 
# seen by WRDS (I cant find patients that are above 70 and not seen by wrds, might need a seperate data)

filtered_egfr <- tbl_labs %>% filter(between(egfr, 30, 60))
egfr_patients_filtered <- left_join(filtered_egfr, filtered_patients, by = "nhi")
# NOTE:
# Out of all the patients that are diagnosed with T1 and T2 diabetes, only 2189 patients have a recorded egfr value between
# 30 and 60 and are seen by WRDS (I can't find patients that are in between and not seen by wrds, might need a seperate data.)
# Try get general dataset






# WHICH PATIENTS WITH DIABETES SHOULD BE TARGETED TO REDUCE HOSPITALISATION
head(joined_table)
project.2 <- joined_table %>% filter(agegroup != "15-24") #Only want patients aged 30+

project.2 <- project.2 %>% rowwise() %>% mutate(num_admissions = sum(ash, ed, ip, op, na.rm = TRUE))
project.2$High_Users <- "No"

#Identifying the high users of secondary care via the method given on the information booklet
for(i in 1:nrow(project.2)){
  if(project.2[i, 37] > 3){
    project.2[i, 38] <- "Yes"
  }
}

high_users <- project.2 %>% filter(High_Users == "Yes")
low_users <- project.2 %>% filter(High_Users == "No")

summary_function <- function(){
  print("Statistic Summary For HBA1C Of Patients That Are Classified As High Users:")
  print(summary(high_users$hba1c))
  print("Statistic Summary For HBA1C Of Patients That Are Classified As Low Users:")
  print(summary(low_users$hba1c))
  print("Statistic Summary For ACR Of Patients That Are Classified As High Users:")
  print(summary(high_users$acr))
  print("Statistic Summary For ACR Of Patients That Are Classified As Low Users:")
  print(summary(low_users$acr))
  print("Statistic Summary For EGFR Of Patients That Are Classified As High Users:")
  print(summary(high_users$egfr))
  print("Statistic Summary For EGFR Of Patients That Are Classified As Low Users:")
  print(summary(low_users$egfr))
}
#Looking at HBA1C, ACR and eGFR for patients who are classified as high users vs low users
project.2 %>% ggplot(aes(x = hba1c, fill = High_Users, group = High_Users)) + geom_bar() + facet_wrap(~High_Users) + theme_ft_rc() +
  xlab("HbA1c") + ylab("Number Of Patients") + ggtitle("HbA1c For Patients Who Are High Users vs. Patients Who Are Not High Users") +
  scale_fill_discrete(name = "High Users?")

project.2 %>% filter(acr < 100) %>% ggplot(aes(x = acr, fill = High_Users, group = High_Users)) + geom_bar() + facet_wrap(~High_Users) + theme_ft_rc() +
  xlab("ACR") + ylab("Number Of Patients") + ggtitle("ACR For Patients Who Are High Users vs. Patients Who Are Not High Users") +
  scale_fill_discrete(name = "High Users?")

project.2 %>% ggplot(aes(x = egfr, fill = High_Users, group = High_Users)) + geom_bar() + facet_wrap(~High_Users) + theme_ft_rc() +
  xlab("eGFR") + ylab("Number Of Patients") + ggtitle("eGFR For Patients Who Are High Users vs. Patients Who Are Not High Users") +
  scale_fill_discrete(name = "High Users?")

summary_function()

#Note almost identical distribution for all three among the three variables

#I just wanted to look at the number of admissions for patients that are classified as high users
project.2 %>% ggplot(aes(x = num_admissions, group = High_Users, fill = High_Users)) + geom_bar() + facet_wrap(~High_Users, scales = "free") +
  scale_fill_manual(values = bar_cols, name = "High Users?") + theme_ft_rc() + xlab("Number Of Admissions") + ylab("Number Of Patients") +
  ggtitle("Number Of Admissions For High User Patients vs. Low User Patients")

project.2$Record_HU <- ifelse(project.2$High_Users == "Yes", 1, 0)
#Mean hba1c by practice
mean_hba1c_overall <- project.2
patPrac <- mean_hba1c_overall %>% group_by(practiceid) %>% tally()
patPrac$pop_prop <- patPrac$n / sum(patPrac$n)
hba1cPractice <- mean_hba1c_overall %>% group_by(practiceid) %>% summarise(mean_hba1c = mean(hba1c, na.rm = TRUE))
mean(hba1cPractice$mean_hba1c)
mean_hba1c_overall %>% ggplot(aes(x = hba1c)) + geom_bar(fill = "blue", colour = "white") + theme_ft_rc() + xlab("Log HBA1C") +
  ylab("Count Of Patients") + ggtitle("Mean HbA1c By Practice")

#Mean hba1c by practice for the maori population
patPrac_maori <- mean_hba1c_overall %>% filter(ethgroup == "Maori")
patPrac_maori <- patPrac_maori %>% group_by(practiceid) %>% summarise(mean_hba1c = mean(hba1c, na.rm = TRUE))
mean(patPrac_maori$mean_hba1c)
mean_hba1c_overall %>% filter(ethgroup == "Maori") %>% ggplot(aes(x = hba1c)) + geom_bar(fill = "red", colour = "white") + theme_ft_rc() + xlab("HBA1C") +
  ylab("Count Of Maori Patients") + ggtitle("Mean HbA1c By Practice For Maori Patients")

#Percentage of maori within each practice
maori_in_practice <- mean_hba1c_overall %>% filter(ethgroup == "Maori")
maori_in_practice <- maori_in_practice %>% group_by(practiceid) %>% tally()
colnames(maori_in_practice)[2] <- "Maori Total"
others_in_practice <- mean_hba1c_overall %>% filter(ethgroup != "Maori")
others_in_practice <- others_in_practice %>% group_by(practiceid) %>% tally()
colnames(others_in_practice)[2] <- "Others Total"

#Join the maori population and other population together
population_comparison <- left_join(maori_in_practice, others_in_practice, by = "practiceid")
population_comparison$maoriInPractice <- population_comparison$`Maori Total` / (population_comparison$`Others Total` + population_comparison$`Maori Total`)
population_comparison$percentageMaori <- population_comparison$maoriInPractice * 100
mean(population_comparison$percentageMaori)
orderedPopulation <- population_comparison[order(population_comparison$percentageMaori),]
orderedPopulation$indicator <- seq(from = 0, to = nrow(orderedPopulation) - 1)
orderedPopulation %>% ggplot(aes(x = reorder(practiceid, indicator), y = percentageMaori)) + geom_col(col = "#000000", fill = "#FFD000") + 
  theme_modern_rc() + ylab("Percentage Of Maori") + xlab("Practice") + ggtitle("Percentage Of Maori Patients Within Each Practice") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rm(orderedPopulation, population_comparison)

#Percentage of high secondary healthcare users
high_users_comparison <- project.2 %>% select(practiceid, hba1c, High_Users, Record_HU)
low_users_comparison <- high_users_comparison %>% filter(High_Users == "No")
high_users_comparison <- high_users_comparison %>% filter(High_Users == "Yes")
#number of high user patients within each practice
just_practiceid_high <- high_users_comparison %>% group_by(practiceid) %>% tally()
colnames(just_practiceid_high)[2] <- "High Users"
#number of low user patients within each practice
just_practiceid_low <- low_users_comparison %>% group_by(practiceid) %>% tally()
colnames(just_practiceid_low)[2] <- "Low Users"

high_users_comparison <- merge(just_practiceid_high, just_practiceid_low, by = "practiceid")
high_users_comparison$percentageOfHighUsers <- high_users_comparison$`High Users` / (high_users_comparison$`High Users` + high_users_comparison$`Low Users`)
high_users_comparison$percentageOfHighUsers <- high_users_comparison$percentageOfHighUsers * 100
high_users_comparison <- high_users_comparison[order(high_users_comparison$percentageOfHighUsers),]
high_users_comparison$indicator <- seq(from = 0, to = nrow(high_users_comparison)-1)
high_users_comparison %>% ggplot(aes(x = reorder(as.factor(practiceid), percentageOfHighUsers), y = percentageOfHighUsers)) + geom_col(colour = "#000000", fill = "#FFD000") + theme_ft_rc() + xlab("Practices") +
  ylab("Percentage Of High Users") + ggtitle("Percentage Of High Users Between Practices") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
mean(high_users_comparison$percentageOfHighUsers)

rm(just_practiceid_high, just_practiceid_low, low_users_comparison)

#mean hba1c readings for patients who are regarded high users within each practice
just_hba1c <- project.2 %>% select(practiceid, hba1c, High_Users, Record_HU)
just_hba1c <- just_hba1c %>% group_by(practiceid) %>% summarise(mean_hba1c = mean(hba1c, na.rm = TRUE))
just_hba1c <- just_hba1c[order(just_hba1c$practiceid),]
just_hba1c$indicator <- seq(0, nrow(just_hba1c)-1)
just_hba1c %>% ggplot(aes(x = reorder(as.factor(practiceid), mean_hba1c), y = mean_hba1c)) + geom_col(colour = "black", fill = "#FFD000") +
  theme_ft_rc() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Practices") + ylab("Mean HbA1c") +
  ggtitle("Mean HbA1c For High User Patients Across The Various Practices")
mean(just_hba1c$mean_hba1c)

rm(just_hba1c)

#################################################### REGRESSION ###########################################################

# NOTE: MIGHT NEED TO ALTER THE DATA A TINY BIT JUST SO THAT IT IS SUITABLE FOR REGRESSION

#### Correlation Plot (Takes a while to compute)
numeric_col <- select_if(project.2, is.numeric)
chart.Correlation(numeric_col, method = "spearman", histogram = TRUE, pch = 16)
#####################################

#Stepwise Regression
#Use a full dataset that doesn't has no na values for regression (use imputed dataset from labs)

#Create copies of original tables
regression_joined <- joined_table

#Set NA to proper values
regression_joined$aucode[is.na(regression_joined$aucode)] <- -1 #(-1 means cant geocode)
regression_joined$gp[is.na(regression_joined$gp)] <- 0
regression_joined$nurse[is.na(regression_joined$nurse)] <- 0
regression_joined$other[is.na(regression_joined$other)] <- 0
regression_joined$ed[is.na(regression_joined$ed)] <- 0
regression_joined$ash[is.na(regression_joined$ash)] <- 0
regression_joined$portal[is.na(regression_joined$portal)] <- 0
regression_joined$op[is.na(regression_joined$op)] <- 0
regression_joined$beddays[is.na(regression_joined$beddays)] <- 0
regression_joined$ip[is.na(regression_joined$ip)] <- 0
#The assumption I made is that NA means that the patients has not seen the variable name

regression_joined <- regression_joined[, -c(28:33)]
regression_joined <- imputedData

#Have to impute systolic pressure and diastolic pressure (I will be using the mice package for imputation again)
# rm(tempData)
# tempData <- mice(regression_joined, m = 10, meth = 'pmm', seed = 500)
# summary(tempData)
# regression_joined <- complete(tempData, 1)
# 
# #Join imputed labs table to the imputed regression table
# regression_joined <- left_join(regression_joined, completedData, by = "nhi")
# regression_joined <- regression_joined[, -c(34:35)]
# colnames(regression_joined)[31] <- 'hba1c'
# colnames(regression_joined)[32] <- 'acr'
# colnames(regression_joined)[33] <- 'egfr'

regression_joined <- regression_joined %>% filter(agegroup != "15-24") #Only want patients aged 30+
regression_joined <- regression_joined %>% select(-c(log_hba1c, log_acr, practiceid, nhi))

regression_joined <- regression_joined %>% rowwise() %>% mutate(num_admissions = sum(ash, ed, ip, op, na.rm = TRUE))
regression_joined$High_Users <- "No"

#Identifying the high users of secondary care via the method given on the information booklet
for(i in 1:nrow(regression_joined)){
  if(regression_joined[i, 31] > 3){
    regression_joined[i, 32] <- "Yes"
  }
}

regression_joined$Record_HU <- if_else(regression_joined$High_Users == "Yes", 1, 0)


################### NOTE: TABLE IS NOW IMPUTED AND READY FOR REGRESSION ANALYSIS ##############################

#Full model (Wont work as the data is too large and the computation power is excessive to the RAM size)
#full.model <- lm(Record_HU ~ ., data = regression_joined)
#Changing the reference level
levels(regression_joined$type)
regression_joined$type <- factor(regression_joined$type, c("T2", "T1", "DM"))
levels(regression_joined$agegroup)
regression_joined$agegroup <- factor(regression_joined$agegroup, c("65+", "25-44", "45-64", "15-24"))
levels(regression_joined$ethgroup)
regression_joined$ethgroup <- factor(regression_joined$ethgroup, c("European", "Asian", "Maori", "Other", "Pacific Island"))

#Start by looking at the null model
null_model <- glm(Record_HU ~ 1, data = regression_joined, family = "binomial")
summary(null_model)

#Notes on the null model:
# 1. A null model is highly related to a null hypothesis
# 2. It is used to jusitfy the use of GLM, as to reach certainty that there is no discrepancy in the second level
# 3. If there is no significant variation at level 2, you do not need to use GLM. If there is a variation, you then try to explain
# it by running the full model. 

#Look at the full model
#NOTE: I am going to do the full model analysis on only 10,000 patients as the size of the data is too big
full.model_table <- regression_joined #%>% slice(1:10000)
full.model_table <- full.model_table[, -c(1:2, 13)]

full_model <- glm(Record_HU ~ ., data = full.model_table, family = "binomial")
summary(full_model)

library(MASS)
step_model <- stepAIC(full_model, direction = "both", trace = FALSE)
summary(step_model)

models <- regsubsets(Record_HU ~ ., data = full.model_table, nvmax = 5, method = "seqrep")
summary(models)

#Setup repeate k-fold cross-validation
train.control <- trainControl(method = "cv", number = 10)
#Train the model
model <- train(Record_HU ~ ., data = full.model_table, method = "lmStepAIC", trace = FALSE,
               trControl = train.control)
model$results
model$bestTune
model$finalModel #Final Model Coefficients
summary(model$finalModel) #Summary of the model

linear.model <- lm(Record_HU ~ factor(type) + factor(agegroup) + factor(ethgroup) + factor(gender) + hba1c + quintile, data = regression_joined)
summary(linear.model)
plt <- gg_diagnose(linear.model, plot.all = FALSE)
names(plt)
plot_all(plt)  #Plotting the diagnostic plots for the linear model

rm(plt)    #Remove the information on the linear model

#Recording High Users Values (0 = No, 1 = Yes)
# regression_joined <- regression_joined[, -c(1:2)]
# regression_joined <- regression_joined[, -c(32:33)]
# regression_joined <- regression_joined[, -c(11)]

################################### JUST TESTING QUINTILE AS A FACTOR ######################################
factored_data <- regression_joined
factored_data$quintile <- as.factor(factored_data$quintile)
levels(factored_data$quintile)
factored_data$quintile <- factor(factored_data$quintile, c("5", "4", "3", "2", "1", "0", "-1"))
quintile_factor_mod <- glm(Record_HU ~ factor(type) + factor(agegroup) + factor(ethgroup) + factor(gender) + hba1c + factor(quintile), data = factored_data, 
            family = binomial(link = "logit"))
summary(quintile_factor_mod)
confint(quintile_factor_mod)

######################################################################################################
log_reg <- glm(HighUser ~ factor(type) + factor(agegroup) + factor(ethgroup) + factor(gender) + hba1c + quintile, data = regression_joined, 
               family = binomial(link = "logit"))
summary(log_reg)
confint(log_reg)
log_reg.diag <- glm.diag(log_reg)
glm.diag.plots(log_reg, log_reg.diag) #Check this model diagnostic plot with Lyn

############################################ USE CLUSTERING ON THE SIGNIFICANT VARIABLES #################################
simplified.glm <- regression_joined %>% select(type, agegroup, ethgroup, quintile, hba1c, HighUser)
simplified.glm.model <- glm(Record_HU ~ type + agegroup + ethgroup + quintile + log10(hba1c), data = simplified.glm, family = binomial)
summary(simplified.glm.model) #This is the simplified version of log_reg model

deviance(log_reg);deviance(simplified.glm.model) #log_reg is the best model in this situation

#Using the log hba1c
log_model.data <- simplified.glm
log_model.data$log_hba1c <- log10(log_model.data$hba1c)
log_model.data$hba1c <- NULL
simplified.log.model <- glm(Record_HU ~ ., data = log_model.data, family = binomial())
summary(simplified.log.model)

deviance(simplified.log.model); deviance(log_reg) #log_reg is the best model in terms of deviance but only by a small margin


#Plotting the hba1c, log_hba1c and standardized hba1c
hba1c_variants <- cbind(simplified.glm$hba1c, log_model.data$log_hba1c)
hba1c_variants <- as.data.frame(hba1c_variants)
colnames(hba1c_variants)[1] <- "hba1c"; colnames(hba1c_variants)[2] <- "log_hba1c"
hba1c_variants$standardized <- 0
for (i in 1:nrow(hba1c_variants)) {
  hba1c_variants[i, 3] <- (hba1c_variants[i,1] - mean(hba1c_variants$hba1c)) / sd(hba1c_variants$hba1c)
}
hba1c_variants$indicator <- seq.int(nrow(hba1c_variants))
hba1c_variants %>% ggplot(aes(x = hba1c)) + geom_bar(fill = "red", colour = "red") + theme_ft_rc() + ggtitle("Raw Imputed HbA1c Readings")
hba1c_variants %>% ggplot(aes(x = log_hba1c)) + geom_bar(fill = "yellow", colour = "yellow") + theme_ft_rc() + ggtitle("Log Of HbA1c")
hba1c_variants %>% ggplot(aes(x = standardized)) + geom_bar(fill = "orange", colour = "orange") + theme_ft_rc() + ggtitle("Standardized HbA1c")


#Clustering
gower_dist <- daisy(simplified.glm, metric = "gower")
gower_mat <- as.matrix(gower_dist)
#Printing the most similar patients
simplified.glm[which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]), arr.ind = TRUE)[1, ], ]
#Printing the most dissimilar patients
simplified.glm[which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]), arr.ind = TRUE)[1, ], ]

sil_width <- c(NA)
for(i in 2:8){
  pam_fit <- pam(gower_dist, diss = TRUE, k = i)
  sil_width[i] <- pam_fit$silinfo$avg.width
}

plot(1:8, sil_width, xlab = "Number Of Clusters", ylab = "Silhouettte Width")
lines(1:8, sil_width)

#Summary of each cluster (7 clusters has the highest silhouette width. 5 is simpler and almost as good. Lets pick k = 5)
k <- 8
pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- cbind(simplified.glm, pam_fit$clustering)
pam_results <- pam_results %>% group_by(`pam_fit$clustering`) %>% do(the_summary = summary(.))
pam_results$the_summary

#Visualisation in a lower dimensional space
tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
tsne_data <- tsne_obj$Y %>% data.frame() %>% setNames(c("X", "Y")) %>% mutate(cluster = factor(pam_fit$clustering))

ggplot(aes(x = X, y = Y), data = tsne_data) + geom_point(aes(color = cluster)) + theme_ft_rc() + 
  ggtitle(expression(atop("Patients Observed In A Lower Dimensional Space", atop(italic("(Clustering Using Gower Distancing)"), "")))) +
  scale_color_discrete(name = "Cluster Group:")


########################################################################################################################
# Cross-Validation For Model Accuracy

#Split the data into training and test set
set.seed(123)
training.samples <- simplified.glm$Record_HU %>% createDataPartition(p = 0.8, list = FALSE)
train.data <- simplified.glm[training.samples, ]
test.data <- simplified.glm[-training.samples, ]
#Make predictions and compute the R2, RMSE and MAE
predictions <- simplified.glm.model %>% predict(test.data)
data.frame(R2 = R2(predictions, test.data$Record_HU),
           RMSE = RMSE(predictions, test.data$Record_HU),
           MAE = MAE(predictions, test.data$Record_HU))
RMSE(predictions, test.data$Record_HU) / mean(test.data$Record_HU) #The prediction error rate

# K-fold cross-validation

#Defining training control
set.seed(123)
train.control <- trainControl(method = "LOOCV")
model_cv <- train(Record_HU ~ factor(type) + factor(agegroup) + factor(ethgroup) + factor(gender) + hba1c + quintile, data = regression_joined, 
               family = binomial(link = "logit"), trControl = train.control)

########### NOTES ON LOGISTIC REGRESSION ######################
# The referrent of a variable is the lowest value.
setwd("C:/Users/jeffr/OneDrive/Desktop/Honours Project")

regression_joined$HighUser <- "NA"
for(i in 1:nrow(regression_joined)){
  if(regression_joined[i, 33] == 1){
    regression_joined[i, 34] <- "Yes"
  } else{
    regression_joined[i, 34] <- "No"
  }
}

#Remove high user binary column
regression_joined$Record_HU <- NULL; regression_joined$High_Users <- NULL
write.csv(regression_joined, "regression_data_v1.csv", row.names = FALSE)

simplified.glm$HighUser <- "NA"
for(i in 1:nrow(simplified.glm)){
  if(simplified.glm[i, 6] == 1){
    simplified.glm[i, 7] <- "Yes"
  } else{
    simplified.glm[i, 7] <- "No"
  }
}

write.csv(simplified.glm, "simplified_glm_v1.csv", row.names = FALSE)


#Comparing attributes for patients that are high users as opposed to patients that are not
#non_users <- regression_joined %>% filter(HighUser == "No"); non_users$Record_HU <- NULL
# (high_users & low_users)
users <- rbind(high_users, low_users)
users$aucode[is.na(users$aucode)] <- -1 #(-1 means cant geocode)
users$gp[is.na(users$gp)] <- 0
users$nurse[is.na(users$nurse)] <- 0
users$other[is.na(users$other)] <- 0
users$ed[is.na(users$ed)] <- 0
users$ash[is.na(users$ash)] <- 0
users$portal[is.na(users$portal)] <- 0
users$op[is.na(users$op)] <- 0
users$beddays[is.na(users$beddays)] <- 0
users$ip[is.na(users$ip)] <- 0

#Updated table of high users and low users
high_uses_v2 <- users %>% filter(High_Users == "Yes")
low_users_v2 <- users %>% filter(High_Users == "No")

# Using test to compare betweent the two groups
users %>% ggplot(aes(x = hba1c, group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of HbA1c For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

users %>% ggplot(aes(x = log(hba1c), group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of Log HbA1c For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

summary(low_users_v2$log_hba1c)
summary(high_uses_v2$log_hba1c)

t.test(log(hba1c) ~ High_Users, data = users) #hba1c

users %>% ggplot(aes(x = num_admissions, group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of Number Of Admission For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

wilcox.test(num_admissions ~ High_Users, data = users)
wilcox.test(acr ~ High_Users, data = users)

users %>% ggplot(aes(x = log(acr), group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of Log AcR For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

users %>% ggplot(aes(x = egfr, group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of eGFR For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

wilcox.test(egfr ~ High_Users, data = users)

users %>% ggplot(aes(x = systolic_pressure, group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of Systolic Pressure For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

wilcox.test(systolic_pressure ~ High_Users, data = users)

wilcox.test(bmi ~ High_Users, data = users)

users %>% ggplot(aes(x = bmi, group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of Systolic Pressure For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

summary(low_users_v2$bmi)
summary(high_uses_v2$bmi)

wilcox.test(hcn ~ High_Users, data = users)

users %>% ggplot(aes(x = systolic_pressure, group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of Systolic Pressure For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

wilcox.test(users$systolic_pressure ~ users$High_Users)

users %>% ggplot(aes(x = diastolic_pressure, group = High_Users, fill = High_Users)) + geom_density(alpha = 0.3, colour = "white") + 
  theme_modern_rc() + scale_fill_discrete(name = "High User?") + ggtitle("Distribution of Diastolic Pressure For High & Low Users Of Secondary Care") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

wilcox.test(users$diastolic_pressure ~ users$High_Users)

# Porpotion test for high users

#Gender 
data.frame(table(high_uses_v2$gender)) #Male = 3587 and Female = 2869
prop.test(x = 3587, n = 6456, p = 0.5, correct = FALSE, alternative = "greater")

#H0: Proportion of T2 patients are the same for high users as opposed to low users
#Ha: Proportion of T2 patients are higher for high users as opposed to low users
data.frame(table(low_users_v2$type)) #Type 2 (High) = 5640 : Type 2 (Low) : 9122
prop.test(x = c(5640, 9122), n = c(6456, 9642), alternative = "greater")

#H0: Proportion of T1 patients are the same for high users as opposed to low users
#Ha: Proportion of T1 patients are higher for high users as opposed to low users
data.frame(table(high_uses_v2$type)) #Type 2 (High) = 148 : Type 2 (Low) : 79
prop.test(x = c(148, 79), n = c(6456, 9642), alternative = "greater")



#Notes:
# eGFR = The best test to measure your level of kidney function and determine your stage of kidney disease.
# acr = Urine albumin to creatinine ratio, helps identify kideny disease that can occur as a complication of diabetes (< 30 is normal)
# hba1c = Indicates a person's average blood glucose (sugar) levels over the previous 2 to 3 months


# Do cluster analysis to find commonalities between patients that are classified as high users of secondary care users



# Project 1: Risk factors for Waikato patients with diabetes not accessing secondary care
library(tidyr)
obesiety$type <- as.factor(obesiety$type)
obesiety$ethgroup <- as.factor(obesiety$ethgroup)
obesiety$agegroup <- as.factor(obesiety$agegroup)
patients <- merge(tbl_patients, tbl_events, by = "nhi")
patients$aucode[is.na(patients$aucode)] <- -1 #(-1 means cant geocode)
patients$gp[is.na(patients$gp)] <- 0
patients$nurse[is.na(patients$nurse)] <- 0
patients$other[is.na(patients$other)] <- 0
patients$ed[is.na(patients$ed)] <- 0
patients$ash[is.na(patients$ash)] <- 0
patients$portal[is.na(patients$portal)] <- 0
patients$op[is.na(patients$op)] <- 0
patients$beddays[is.na(patients$beddays)] <- 0
patients$ip[is.na(patients$ip)] <- 0

#Criteria for WRDS
generalDiabetes <- patients %>% filter(type == "T1")
youthGroup <- patients %>% filter(agegroup == "15-24")


patient_lab <- merge(tbl_labs, tbl_patients, by = "nhi")
removed <- patient_lab %>% drop_na(hba1c, acr, egfr)
new_DF <- patient_lab[is.na(patient_lab$hba1c),]

#tbl_patients = 21310
#tbl_events = 21310
#tbh_measurements = 21311
#tbl_labs = 21348


#Obesiety
imputedData$log_acr <- log10(imputedData$acr)
imputedData$log_hba1c <- log10(imputedData$hba1c)

obesiety <- imputedData


#BMI of patients by age
test <- obesiety %>% group_by(agegroup, specific, type) %>% summarise(counts = n())
test %>%
  ggplot(aes(x = specific, y = counts, group = type, fill = type)) + geom_col(position = "dodge", colour = "white") + facet_wrap(~ agegroup, scale = "free") +
  theme_modern_rc() + xlab("Number Of Patients") + ylab("BMI Category") + ggtitle("Effect Of The Type Of Diabetes On The Obesity Of Patients Within Age Groups") +
  theme(strip.background =element_rect(fill=NA)) + theme(strip.text = element_text(colour = 'white'))
age <- obesiety %>% group_by(agegroup, bmi) %>% summarise(counts = n())
ggplot(age, aes(bmi, colour = agegroup, group = agegroup)) + geom_density(fill = NA) + theme_modern_rc() +
  ggtitle("Distribution of BMI Across The Different Age Groups") + scale_color_discrete(name = "Age Group:")
ggplot(age, aes(bmi, colour = agegroup, group = agegroup)) + geom_boxplot(fill = NA) + theme_modern_rc() +
  ggtitle("Box-Plot of BMI Across The Different Age Groups") + scale_fill_discrete(name = "Age Group:") + xlab("BMI")
tapply(age$bmi, age$agegroup, summary)

pairwise.t.test(obesiety$bmi, obesiety$agegroup, p.adjust.method = "bonferroni", paired = FALSE, pool.sd = FALSE)

#BMI of patients by diabetes type
type <- obesiety %>% group_by(type, bmi) %>% summarise(counts = n())
ggplot(type, aes(bmi, colour = type, group = type)) + geom_density(fill = NA) + theme_modern_rc() +
  ggtitle("Distribution of BMI Across The Different Diabetes Type") 
ggplot(type, aes(bmi, colour = type, group = type)) + geom_boxplot(fill = NA) + theme_modern_rc() +
  ggtitle("Box-plot of BMI Across The Different Diabetes Type") + scale_fill_discrete(name = "Diabetes Group:") + xlab("BMI")
tapply(type$bmi, type$type, summary)

#Pairwise T-tests for multiple groups
pairwise.t.test(obesiety$bmi, obesiety$type, p.adjust.method = "bonferroni")

#BMI of patients by ethnicity
#test <- obesiety %>% group_by(ethgroup, agegroup, specific) %>% summarise(count = n())
x <- obesiety %>% select(ethgroup, gender, bmi) 
x$gender <- as.factor(x$gender)
x <- x %>% filter(ethgroup == "Pacific Island")
tapply(x$bmi, x$gender, summary)

x %>% ggplot(aes(x = bmi, colour = gender, group = gender)) + geom_boxplot(fill = NA) + facet_wrap(~ethgroup, scale = "free") + theme_modern_rc() +
  theme(strip.background =element_rect(fill=NA)) + theme(strip.text = element_text(colour = 'white')) + xlab("BMI") + ggtitle("Distribution Of BMI For Gender Based On Ethnicity") +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + scale_colour_discrete(name = "Gender:")

test <- obesiety %>% group_by(ethgroup, specific, gender) %>% summarise(count = n())
test <- test %>% filter(ethgroup == "Pacific Island")

eth_perv <- test %>% filter(specific %in% c("Overweight", "Obese"))
eth_perv$prevalence <- "NULL"

eth_perv$prevalence <- c(2.91, 0.73, 96.36, 3.30, 0.77, 95.92, 5.02, 1.03, 93.95,
                         9.20, 1.50, 89.30, 3.78, 0.83, 95.39, 6.94, 1.73, 91.32,
                         1.01, 98.99, 6.02, 2.41, 91.57, 4.25, 1.31, 94.44, 6.14,
                         1.75, 92.11)
eth_perv$count <- "NULL"
eth_perv$count <- c(24.69, 24.69, 24.69, 46.23, 46.23, 46.23, 53.76, 53.76, 53.76,
                    32.79, 32.79, 32.79, 69.10, 69.10, 69.10, 23.84, 23.84, 23.84,
                    46.92, 46.92, 39.34, 39.34, 39.34, 66.09, 66.09, 66.09,
                    24.62, 24.62, 24.62)

eth_perv %>% ggplot(aes(x = ethgroup, y = count, fill = type, group = type)) + geom_col(position = "dodge")

for(i in 1:nrow(eth_perv)){
  if(eth_perv[i, 1] == "Asian"){
    eth_perv[i, 4] = eth_perv[i, 3] / 1114
  }else if (eth_perv[i, 1] == "European"){
    eth_perv[i, 4] = eth_perv[i, 3] / 11374
  }else if (eth_perv[i, 1] == "Maori"){
    eth_perv[i, 4] = eth_perv[i, 3] / 3142
  }else if(eth_perv[i, 1] == "Other"){
    eth_perv[i, 4] = eth_perv[i, 3] / 211
  }else{
    eth_perv[i, 4] = eth_perv[i, 3] / 463
  }
}
eth_perv$prevalence <- as.double(eth_perv$prevalence)
eth_perv$prevalence <- eth_perv$prevalence * 100; eth_perv$prevalence <- format(round(eth_perv$prevalence, 2), nsmall = 2)
eth_perv$prevalence <- as.double(eth_perv$prevalence)
eth_perv %>%
  ggplot(aes(x = ethgroup, y = prevalence, fill = type, group = type)) + geom_col(colour = NA) + theme_modern_rc() + xlab("Ethnic Group") + ylab("Prevalence (%)") + 
  scale_fill_discrete(name = "BMI High End Category:") + ggtitle("Ethnicity-Obesity Prevalence For Obese And Overweight Patients") + facet_wrap(~specific) + 
  theme(strip.text = element_text(colour = 'white'))



ethnicity <- obesiety %>% group_by(ethgroup, bmi) %>% summarise(counts = n())
ggplot(ethnicity, aes(bmi, colour = ethgroup, group = ethgroup)) + geom_density(fill = NA) + theme_modern_rc() +
  ggtitle("Distribution of BMI Across The Different Ethnic Groups") 
ggplot(ethnicity, aes(bmi, colour = ethgroup, group = ethgroup)) + geom_boxplot(fill = NA) + theme_modern_rc() +
  ggtitle("Box-plot of BMI Across The Different Ethnicity Groups") + scale_fill_discrete(name = "Ethnicity Group:") + xlab("BMI")
tapply(ethnicity$bmi, ethnicity$ethgroup, summary)

pairwise.t.test(obesiety$bmi, obesiety$ethgroup, p.adjust.method = "bonferroni", paired = FALSE, pool.sd = FALSE)


#Prevalence
obesiety$specific <- "Test"
for(i in 1:nrow(obesiety)){
  if (obesiety[i, 23] > 29.9){
    obesiety[i, 35] <- "Obese"
  } else if(between(obesiety[i, 23], 25, 29.9)){
    obesiety[i, 35] <- "Overweight"
  } else if(between(obesiety[i, 23], 18.5, 24.9)){
    obesiety[i, 35] <- "Normal"
  } else{
    obesiety[i, 35] <- "Underweight"
  }
}
obesiety$specific <- as.factor(obesiety$specific)


overweightAgePervalence <- obesiety %>% select(agegroup, ethgroup, specific)
overweightAgePervalence <- overweightAgePervalence %>% filter(specific %in% c("Overweight", "Obese"))
#overweightAgePervalence <- overweightAgePervalence %>% group_by(agegroup, ethgroup) %>% summarise(count = n())
overweightAgePervalence <- overweightAgePervalence %>% group_by(agegroup, specific) %>% summarise(count = n())
overweightAgePervalence$prevalence <- "NULL"
for(i in 1:nrow(overweightAgePervalence)){
  if(overweightAgePervalence[i, 1] == "15-24"){
    overweightAgePervalence[i, 4] <- overweightAgePervalence[i, 3] / 108
  } else if(overweightAgePervalence[i, 1] == "25-44"){
    overweightAgePervalence[i, 4] <- overweightAgePervalence[i, 3] / 1008
  } else if(overweightAgePervalence[i, 1] == "45-64"){
    overweightAgePervalence[i, 4] <- overweightAgePervalence[i, 3] / 5084
  } else{
    overweightAgePervalence[i, 4] <- overweightAgePervalence[i, 3] / 7956
  }
}

for(i in 1:nrow(overweightAgePervalence)){
  overweightAgePervalence[i, 3] <- overweightAgePervalence[i, 2] / 19
}

overweightAgePervalence_ <- data.frame("agegroup" = c("15-24", "15-24", "25-44","25-44", "45-64", "45-64", "65+", "65+"), 
                                       "specific" = c("Overweight", "Obese", "Overweight", "Obese", "Overweight", "Obese", "Overweight", "Obese"), 
                                       "count" = c(49, 59, 340, 669, 1541, 3543, 3261, 4695), 
                                       "prevalence" = c(23.87, 28.64, 27.27, 53.64, 27.59, 63.44, 35.19, 50.67))
overweightAgePervalence$prevalence <- as.double(overweightAgePervalence$prevalence)
overweightAgePervalence$prevalence <- overweightAgePervalence$prevalence * 100
overweightAgePervalence %>%
  ggplot(aes(x = agegroup, y = prevalence, colour = ethgroup, group = ethgroup)) + geom_line() + theme_modern_rc() +
  xlab("Age Group") + ylab("Prevalence (%)") + scale_color_discrete(name = "Ethnicity Group:") + 
  ggtitle("Age-Specific Prevalence of Obesity For The Different Ethnicity Groups")

overweightAgePervalence_ %>%
  ggplot(aes(x = agegroup, y = prevalence, fill = specific, group = specific)) + geom_col(col = "white") + theme_modern_rc() + xlab("Age Group") + ylab("Prevalence (%)") + 
  scale_fill_discrete(name = "BMI High End Category:") + ggtitle("Age-Obesity Prevalence For Obese And Overweight Patients")


smoker <- obesiety %>% select(currentsmoker, smoker, specific)
smoker <- smoker %>% filter(smoker == "Yes")

smoker <- smoker %>% group_by(smoker, specific) %>% summarise(count = n())
smoker$prevalence <- "NULL"
for(i in 1:nrow(smoker)){
  if(smoker[i, 2] == "Normal"){
    smoker[i, 4] <- smoker[i, 3] / 2079
  } else if(smoker[i, 2] == "Obese"){
    smoker[i, 4] <- smoker[i, 3] / 8966
  } else if(smoker[i, 2] == "Overweight"){
    smoker[i, 4] <- smoker[i, 3] / 5191
  } else{
    smoker[i, 4] <- smoker[i, 3] / 68
  }
}
smoker <- na.omit(smoker)
smoker$prevalence <- as.double(smoker$prevalence)
smoker$prevalence <- smoker$prevalence * 100

smoker %>%
  ggplot(aes(x = specific, y = prevalence, fill = smoker, group = smoker)) + geom_col(stat = "identity", col = "white", position = "dodge") +
  theme_modern_rc() + xlab("BMI Group") + ylab("Prevalence (%)") + scale_fill_discrete(name = "Smoker:") + 
  ggtitle("Percentage of Smokers & Non-Smokers Within BMI Group") +
  labs(caption = "NOTE: Cross sectional study suggests that heavy smoking could be associated with higher risk of obesiety. (https://academic.oup.com/ajcn/article/87/4/801/4633357)")


#Looking at both the genders
test <- obesiety %>% group_by(gender, type) %>% summarise(count = n())
test %>% ggplot(aes(x = gender, y = count, fill = type, group = type)) + geom_col(colour = "white", position = "dodge") + theme_modern_rc() +
  xlab("Number Of Patients") + ylab("Gender") + ggtitle("Number Of Patients Within Each Diabetic Group By Gender") + scale_fill_discrete(name = "Type:")
genders_obesity <- obesiety %>% select(gender, specific)
genders_obesity <- genders_obesity %>% group_by(gender, specific) %>% summarise(count = n())
genders_obesity %>% ggplot(aes(x = specific, y = count, group = gender, fill = gender)) + geom_col(position = "dodge", color = "white") + theme_modern_rc() +
  xlab("BMI Category") + ylab("Number of Patients") + ggtitle("Number of patients within each BMI category") + scale_fill_discrete(name = "Gender:")
genders_obesity$prevalence <- "NULL"
for(i in 1:nrow(genders_obesity)){
  if(genders_obesity[i, 1] == "Female"){
    genders_obesity[i, 4] <- genders_obesity[i, 3] / 7049
  } else{
    genders_obesity[i, 4] <- genders_obesity[i, 3] / 9255
  }
}
genders_obesity$prevalence <- as.double(genders_obesity$prevalence)
genders_obesity$prevalence <- genders_obesity$prevalence * 100

genders_obesity %>% ggplot(aes(x = specific, y = prevalence, group = gender, colour= gender)) + geom_line() + theme_modern_rc() +
  xlab("BMI Category") + ylab("Prevalence (%)") + scale_color_discrete(name = "Gender:") + 
  ggtitle("Within Gender Sample Percenage Within BMI Category")


# Looking at obesity by high users and low users
obesiety$num_admissions <- NULL
obesiety$HighUser <- "NULL"
for(i in 1:nrow(obesiety)){
  if(obesiety[i, 27] + obesiety[i, 28] + obesiety[i, 30] + obesiety[i, 32] > 3){
    obesiety[i, 36] <- "Yes"
  } else{
    obesiety[i, 36] <- "No"
  }
}

highUser_obesity <- obesiety %>% select(HighUser, specific)
highUser_obesity <- highUser_obesity %>% group_by(HighUser, specific) %>% summarise(count = n())
highUser_obesity %>% ggplot(aes(x = specific, y = count, group = HighUser, fill = HighUser)) + geom_col(position = "dodge", color = "white") + theme_modern_rc() +
  xlab("BMI Category") + ylab("Number of Patients") + ggtitle("Number of patients within each BMI category (High User vs. Non-High Users)") + scale_fill_discrete(name = "Gender:")
highUser_obesity$prevalence <- "NULL"
for(i in 1:nrow(highUser_obesity)){
  if(highUser_obesity[i, 1] == "No"){
    highUser_obesity[i, 4] <- highUser_obesity[i, 3] / 9685
  } else{
    highUser_obesity[i, 4] <- highUser_obesity[i, 3] / 6619
  }
}
highUser_obesity$prevalence <- as.double(highUser_obesity$prevalence)
highUser_obesity$prevalence <- highUser_obesity$prevalence * 100

highUser_obesity %>% ggplot(aes(x = specific, y = prevalence, group = HighUser, colour= HighUser)) + geom_line() + theme_modern_rc() +
  xlab("BMI Category") + ylab("Prevalence (%)") + scale_color_discrete(name = "High User:") + 
  ggtitle("Within High User Sample Percenage Within BMI Category")

highUser_obesity$prevalence <- round(highUser_obesity$prevalence, 2)
highUser_obesity %>% 
  ggplot(aes(x = HighUser, y = prevalence, group = specific, fill = specific)) + geom_col(position = "dodge", colour = "white") + theme_modern_rc() +
  xlab("High User?") + ylab("Prevalence (%)") + scale_fill_discrete(name = "BMI Category:") + 
  geom_text(aes(label= paste0(prevalence, "%")), position=position_dodge(width=0.9), vjust=-0.25, colour = "white") + 
  ggtitle("Percentage of Patients Within Each Category By Secondary User")
       
#Disconnecting to the database
dbDisconnect(con)
