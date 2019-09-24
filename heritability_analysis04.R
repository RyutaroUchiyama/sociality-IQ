

require(reshape2)
require(lme4)
###require(car)
###library(ggplot2)
###library(stargazer) 

###qqp(dat.long$heritability, "norm")
# Q-Q plots show that heritability is normally distributed


# Import Uz (2015) tightness-looseness data 
uz <- read.csv("~/Desktop/Muthukrishna_Lab/sociality-IQ/DATA/uz_tightness_withISO.csv")
# Import Varieties of Democracy 'V-Dem Full+Others' data set
vdem <- readRDS("~/Desktop/Muthukrishna_Lab/sociality-IQ/DATA/Country_Year_V-Dem_Full+others_R_v9/V-Dem-CY-Full+others-v9.rds")
# Set working directory to one containing all MaTCH ICF/ICD10 subchapter (by country) tables 
setwd("~/Desktop/Muthukrishna_Lab/sociality-IQ/DATA/MaTCH all traits no constraints")

# Extract all country names from all MaTCH tables  
filenames <- list.files()
country <- as.character(matrix(NA, 1, length(filenames)))
for (i in 1:length(filenames)){
  file <- read.csv(filenames[i], sep=";")
  if (dim(file)[1]>0){
    country[i] <- as.character(read.csv(filenames[i], sep=";")$Code) }}
country <- unique(country)
country <- as.factor(country[!is.na(country)])


# Extract h2_all values from all MaTCH tables, collate into one table where rows are countries and columns are subchapter traits
dat.pre <- data.frame(country)
for (i in 1:length(filenames)){
  import <- read.csv(filenames[i], sep=";")[,c(1,21)]
  import$h2_all <- replace(import$h2_all, which(import$h2_all==0), NA) # where heritability is 0, replace with NA, assuming that all instances of 0 are in fact NA  
  colnames(import)[2] <- substr(filenames[i], 1, nchar(filenames[i])-4 ) # remove '.csv' from file name to turn into variable name, and make it the name of that column 
  dat.pre <- merge(dat.pre, import, by.x="country", by.y="Code", all.x=T)
  }

# Remove columns containing no data for h2_all, left with 73 traits (first column is just countries)
columns.nodata <- matrix(0,0,0)
for (i in 1:dim(dat.pre)[2]){
if (sum(!is.na(dat.pre[,i])) == 0){
  columns.nodata <- c(columns.nodata, i)
  }}
dat.pre <- dat.pre[,-columns.nodata]

# Merge Uz 2015 tightness-looseness data with MaTCH data
dat <- merge(uz[uz$CountryISO %in% country,], dat.pre, by.x="CountryISO", by.y="country")
# **** Uz contains all MaTCH countries except AU (Australia) & NO (Norway), so these get dropped from data.
# **** It might be possible to compute or approximate Uz values for these countries on the basis of WVS 
dat$CTL_DS <- as.numeric(as.character(dat$CTL_DS)) # CTL values imported as factors so converting into numeric 
dat$CTL_DG <- as.numeric(as.character(dat$CTL_DG))
dat$CTL_C <- as.numeric(as.character(dat$CTL_C))

# subset of V-Dem dataset, for 2006 (median year of all papers included in MaTCH meta-analysis) and for countries remaining from MaTCH data
vdem.sub <- rbind(vdem[which(vdem$country_name == "Belgium" & vdem$year == 2006),], 
                  vdem[which(vdem$country_name == "Canada" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "China" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Czech Republic" & vdem$year == 2006),],    
                  vdem[which(vdem$country_name == "Denmark" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Finland" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "France" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "United Kingdom" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "India" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Italy" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Japan" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "South Korea" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Netherlands" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Russia" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Sweden" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "United States of America" & vdem$year == 2006),]
) 
# further subsetting, taking only variables of interest from fill V-Dem dataset
col.select <- c(  
which(colnames(vdem.sub) == "country_name"),
which(colnames(vdem.sub) == "e_migdppc"), #GDP per capita
which(colnames(vdem.sub) == "e_peaveduc"), # years education for citizens older than 15 (missing for Croatia)
which(colnames(vdem.sub) == "e_wb_pop"), # population
which(colnames(vdem.sub) == "e_regionpol") # region (politico-geographic). Less sub-divisions than 'geographic region' (e_regiongeo)
)
dat <- cbind(dat, vdem.sub[,col.select])
dat <- dat[,-which(names(dat)== "CountryISO")] # remove columns for "CountryISO" and "Country", as we use "country_name" in the final data frame
dat <- dat[,-which(names(dat)== "Country")]

dat.long <- melt(dat, na.rm=T, id.vars=c("country_name", "CTL_DS", "CTL_DG", "CTL_C", "e_migdppc", "e_peaveduc", "e_wb_pop", "e_regionpol"))
colnames(dat.long)[c(9,10)] <- c("trait","heritability") # renaming "variable" and "value"

lmm <- lmer(heritability ~ CTL_DG+ e_migdppc+ e_peaveduc+ e_wb_pop+ trait+ (1+ trait+ CTL_DG| e_regionpol))



