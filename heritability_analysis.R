


assign("last.warning", NULL, envir = baseenv())
# reset errors (sigh)


# This script set for 145 traits, 


load(url("https://joachim-gassen.github.io/data/rdf_ests.RData"))


if (!require(tidyverse)) {
  install.packages('tidyverse')
}
if (!require(purrr)) {
  install.packages('purrr')
}
if (!require(broom)) {
  install.packages('broom')
}
if (!require(cowplot)) {
  install.packages('cowplot')
}


#----
library(devtools)
devtools::install_github("joachim-gassen/rdfanalysis")
library(rdfanalysis)

library(reshape2)
library(lme4)
library(readxl)
library(readstata13) 
library(dplyr)
###library(car)
library(ggplot2)
library(stargazer) 
library(stringr)
library(MuMIn)
library(lmerTest)
#library(formattable)

# Import WVS data used for Muthukrishna et al. Cultural Distance analysis
### wvs.alleles <- read.csv("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/wvs_allelesv02_merged.csv")
# Import raw (pre-'allelized') WVS data used for Muthukrishna et al. Cultural Distance analysis
wvs <- read.dta13("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/WVS_Longitudinal_1981_2014_stata_v2015_04_18.dta")
# Import EVS 2008 data
evs <- read.dta13("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural-evolution-genetic-heritability/DATA/EVS longitudinal dataset/ZA4804_v3-0-0.dta")
# Import WVS variable characteristics used for Muthukrishna et al. Cultural Distance analysis
dimensions <- read.csv("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/Additional_Supplementary_Materials/allele-dimensions-data.csv")
# Import list of WVS variabled used in Muthukrishna et al. 
wvsvariables <- read.csv("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/Additional_Supplementary_Materials/allele-dimensions-data.csv")
# Import Uz (2015) tightness-looseness data 
uz <- read.csv("~/Desktop/Muthukrishna_Lab/cultural-evolution-genetic-heritability/DATA/uz_tightness_withISO.csv")
# Import Varieties of Democracy 'V-Dem Full+Others' data set
vdem <- readRDS("~/Desktop/Muthukrishna_Lab/cultural-evolution-genetic-heritability/DATA/Country_Year_V-Dem_Full+others_R_v9/V-Dem-CY-Full+others-v9.rds")
# Import table 33 of Polderman et al. (2013), in supplementary material (.xlsx file) 
trait_hierarchy <- read_excel("~/Desktop/Muthukrishna_Lab/cultural-evolution-genetic-heritability/DATA/Polderman2013 - heritability metaanalysis/polderman2013-supp2.xlsx", sheet = "Supplementary Table 33", range = "A2:C328")
# save default directory
original_wd <- "/Users/ryutaro/Documents/R"
# Set to directory containing all MaTCH ICF/ICD10 subchapter (by country) tables 
setwd("~/Desktop/Muthukrishna_Lab/cultural-evolution-genetic-heritability/DATA/MaTCH all traits no constraints")
# Enter countries whose cultural variance scores you want to use in analysis (whether for pairing with heritability or for interpolation of other countries)
relevant.countries <- c("Australia", "Belgium","Brazil","Canada", "China","Croatia","Czech Republic","Denmark","Finland", "France","Gambia","Germany","Great Britain","Greece","Hungary","India","Israel","Italy","Jamaica","Japan","Jordan","Lebanon",
                        "Namibia","Netherlands","Norway","Poland","Russia","South Korea","Spain","Sri Lanka","Sweden","Taiwan","United States")
# Same as above; make sure the 2-letter ISOs are matching the above country vector
relavant.iso <- c("AU","BE","BR","CA","CN","HR","CZ","DK","FI","FR","GM","DE","GB","GR","HU","IN","IL","IT","JM","JP","JO","LB","NA","NL","NO","PL","RU","KR","ES","LK","SE","TW","US")


#°·°·°·°·°·°·°· Organise heritability data °·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·

# Extract all country codes included in all country-based tables available for download from MaTCH  
filenames <- list.files()
country <- as.character(list(NA))
for (i in 1:length(filenames)){
  file <- read.csv(filenames[i], sep=";")
#  if (dim(file)[1]>0){
    country.add <- as.character(read.csv(filenames[i], sep=";")$Code)
    country <- union(country, country.add)
}
country <- unique(country)
country <- as.factor(country[!is.na(country)])

#prior list: AU FI DK US CA GB NL SE IT CZ IN BE JP NO RU CN FR KR

# Extract h2 values (for each of the 4 h2 types) from all MaTCH tables, collate into one table where rows are countries and columns are subchapter traits
dat.pre.h2_all <- data.frame(country)
dat.pre.h2_ss <- data.frame(country)
dat.pre.h2_m <- data.frame(country)
dat.pre.h2_f <- data.frame(country)
for (i in 1:length(filenames)){
  import <- read.csv(filenames[i], sep=";")[,c(1,21,23,25,27)]
  if (all(is.na(import)) == FALSE){ # if there are numerical values in the matrix (needed because there are some matrices with none)
  import[import==0] <- NA } # where heritability is 0, replace with NA, because importing raw data causes blank cells to become 0. (We are assuming that all instances of 0 are in fact NA)  
  import.traitname <- substr(filenames[i], 9, nchar(filenames[i])-4 ) # delete 'country_' as well as '.csv' from file name to turn into variable name, and make it the name of that column 
  dat.pre.h2_all <- merge(dat.pre.h2_all, import[,c(1,2)], by.x="country", by.y="Code", all.x=T)
  dat.pre.h2_ss <- merge(dat.pre.h2_ss, import[,c(1,3)], by.x="country", by.y="Code", all.x=T)
  dat.pre.h2_m <- merge(dat.pre.h2_m, import[,c(1,4)], by.x="country", by.y="Code", all.x=T)
  dat.pre.h2_f <- merge(dat.pre.h2_f, import[,c(1,5)], by.x="country", by.y="Code", all.x=T)
  colnames(dat.pre.h2_all)[i+1] <- import.traitname # rename column as name of traits as extracted from filename 
  colnames(dat.pre.h2_ss)[i+1] <- import.traitname
  colnames(dat.pre.h2_m)[i+1] <- import.traitname
  colnames(dat.pre.h2_f)[i+1] <- import.traitname
  }
setwd(original_wd) # done reading files, so setting working directory back to what it was initially 

dat.pre <-  rbind(dat.pre.h2_all, dat.pre.h2_ss, dat.pre.h2_m, dat.pre.h2_f)
dat.pre <- cbind(rep(c("h2_all","h2_ss","h2_m","h2_f"), each= dim(dat.pre.h2_all)[1])  , dat.pre)
colnames(dat.pre)[1] <- "h2_type"
 
# Remove columns containing no data in any of the h2 types - left with 121 traits (first column is just countries)
columns.nodata <- matrix(0,0,0)
for (i in 1:dim(dat.pre)[2]){
  if (sum(!is.na(dat.pre[,i])) == 0){ columns.nodata <- c(columns.nodata, i)}}
dat.pre <- dat.pre[,-columns.nodata] #633 h2 estimates

#°·°·°·°·°·°·°· Organise WVS & EVS data °·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·

include_vars145 <- as.character(wvsvariables$WVS_Var[which(wvsvariables$"CAT.All"==1)]) # All variables that appear to be actually used in Muthukrishna et al. analysis
include_vars191 <-  c("A001","A002","A003","A004","A005","A006","A007","A025","A026","A029","A030","A032","A034","A035","A038","A039","A040","A041","A042","A165","A168","A057","A058","A059","A060","A061","A062","B001","B002","B003","B008","B009","A169",
                      "A064","A065","A066","A067","A068","A069","A070","A071","A072","A073","A074","A075","A076","A077","A079","A081","A082","A083","A084","A085","A086","A087","A088","A089","A090","A091","A092","A093","A094","C001","C002","C006","A170","A173",
                      "C008","C009","C010","C011","C012","C013","C014","C015","C016","C017","C018","C019","C020","C021","C036","C037","C038","C039","C040","C041","C059","C060","C061","X007","X011","D017","D018","D019","D022","D023","D054","D055","D056","D057","D058","D059","D060",
                      "E001","E002","E003","E004","E005","E006","E012","E014","E015","E016","E018","E019","E022","E023","E025","E026","E027","E028","E029","E033","E034","E035","E036","E037","E039","E046","E143","E112","E114","E115","E116","E117","E110",
                      "E120","E121","E122","E123","E124","E125","E128","E129","E135","E136","E137","E138","E139","F001","F022","F025","F028","F034","F035","F036","F037","F038","F050","F051","F052","F053","F054","F063","F064","F065","F066","F102","F103",
                      "F104","F105","F114","F115","F116","F117","F118","F119","F120","F121","F122","F123","G001","G002","G006","E150","G015","G016","E179","E180","E182")

#wvs <- wvs[,(names(wvs) %in% include_vars191 | names(wvs) %in% c("S001", "S002", "S003"))] # wvs/evs, year, country; but region not included ("X048WVS") because that column not present in EVS data
wvs <- wvs[,(names(wvs) %in% include_vars145 | names(wvs) %in% c("S001", "S002", "S003"))] # wvs/evs, year, country; but region not included ("X048WVS") because that column not present in EVS data
evs.uniqnames <- setdiff(names(evs), names(wvs))
evs <- evs[, -which(names(evs) %in% evs.uniqnames)] # remove variables that are unique to EVS 
wvs.uniqnames <- setdiff(names(wvs), names(evs))
###wvs <- cbind(wvs, matrix(NA,dim(wvs)[1], length(evs.uniqnames)) )
if (length(wvs.uniqnames)>0){
  evs <- cbind(evs, matrix(NA,dim(evs)[1], length(wvs.uniqnames))) # add empty columns that only exist in WSV, for rbind
  colnames(evs)[ (dim(evs)[2]-length(wvs.uniqnames)+1) : dim(evs)[2]] <- wvs.uniqnames # give names to new empty columns
}

colnames(wvs)[which(names(wvs) == "S001")] <- "survey"
colnames(wvs)[which(names(wvs) == "S002")] <- "wave"
colnames(wvs)[which(names(wvs) == "S003")] <- "pop"
colnames(evs)[which(names(evs) == "S001")] <- "survey"
colnames(evs)[which(names(evs) == "S002")] <- "wave"
colnames(evs)[which(names(evs) == "S003")] <- "pop"
# wvs$wave[which(wvs$survey=="EVS")] <- rep("2008", dim(evs)[1]) # insert wave for evs
wvs <- wvs[,order(names(wvs))]
col_idx <- c(which(names(wvs)=="survey"), which(names(wvs)=="wave"), which(names(wvs)=="pop")) # To move these variables to the front
N.extravars <- length(col_idx)
wvs <- wvs[,c( col_idx, setdiff(1:length(names(wvs)),col_idx) ) ]

for (i in (length(col_idx)+1):dim(wvs)[2]){ # convert all responses into integer, except for above 3 participant attributes
wvs[,i] <- as.integer(wvs[,i])
evs[,i] <- as.integer(evs[,i])
}
wevs <- rbind(wvs, evs) # merge wvs and evs
wevs <- wevs[str_detect(wevs$wave,"^2005") | str_detect(wevs$wave,"^2010") | str_detect(wevs$survey,"EVS"),] # use just these WVS waves (or any EVS observation, all from same wave)
wvsonly <- wvs
wvsonly <- wvsonly[str_detect(wvsonly$wave,"^2005") | str_detect(wvsonly$wave,"^2010"),] # use just these WVS waves (or any EVS observation, all from same wave)


# °·°·°·°·°·°·°· Cultural variance for WVS-EVS integrated data °·°·°·°·°·°·°·  

positiveresponses.wevs <- rep(NA,dim(wevs)[2])
for (i in (N.extravars+1):dim(wevs)[2]){ positiveresponses.wevs[i] <- sum(wevs[,i]>5, na.rm=T)} 
wevs.noresponses <- names(wevs)[which(positiveresponses.wevs %in% 0)] # variables with no positive responses (i.e., >5 in this data, or positive number in WVS codebook)
# ^ EMPTY COLUMNS--24 A057 no reponses; 54 C008, no responses

# NOMINAL VARIABLES in the set of 145 WVS variables——50 B008;  55 C009; 74 C061; 88 E001; 89 E003; 90 E005; 122 F022; 
# (A025 & A026 on the verge of nominal, depends on whether we treat "neither" as in-between the two other responses. Treating here as ordinal.)
wevs145.nominalvars <- c("A025","A026","B008","C009","C060","E001","E003","E005","F022")
#wevs191.nominalvars <- c("A025","A026","B008","C009","C010","C060","E001","E002","E003","E004","E005","E006","E135","E136","E137","E138","E139","F022", "X007")
wevs.removal <- union(wevs145.nominalvars, wevs.noresponses)
wevs <- wevs[,-which(names(wevs) %in% wevs.removal)] # remove columns of nominal and empty columns

# Borderline nominal: A026 (parent's responsibility), 

# # ===== Recoding of factors to make ordinal (necessary for rescaling of variables) =====
# # 12 A025 [6 (112882 "always respect") to21, 7(43399 "respect if earned") to23, 8(74 "neither") to22]
# # Putting "neither" in between "always respect" and "respect if earned"
# k <- which(names(wevs)=="A025")
# wevs[which(wevs[,k]==6), k] <- rep(21, length(which(wevs[,k]==6)))
# wevs[which(wevs[,k]==7), k] <- rep(23, length(which(wevs[,k]==7)))
# wevs[which(wevs[,k]==8), k] <- rep(22, length(which(wevs[,k]==8)))
# 
# # 13 A026 [8 (14552 "neither") to22, 10(114122 "do best") to21, 11(31492 "not sacrifice") to 23]
# # Putting "neither" in between "do best for children" and "should not sacrifice"
# k <- which(names(wevs)=="A026")
# wevs[which(wevs[,k]==8), k] <- rep(22, length(which(wevs[,k]==8)))
# wevs[which(wevs[,k]==10), k] <- rep(21, length(which(wevs[,k]==10)))
# wevs[which(wevs[,k]==11), k] <- rep(23, length(which(wevs[,k]==11)))

# 52 C001; [6(98446) to21, 7(164408) to23, 8(42914) to 22, ]
# Putting "neither" in between "agree" and "disagree"
k <- which(names(wevs)=="C001")
wevs[which(wevs[,k]==6), k] <- rep(21, length(which(wevs[,k]==6)))
wevs[which(wevs[,k]==7), k] <- rep(23, length(which(wevs[,k]==7)))
wevs[which(wevs[,k]==8), k] <- rep(22, length(which(wevs[,k]==8)))

# 53 C002; [6(203625) to21, 7(60178) to23, 8(32795) to22, ]
# Putting "neither" in between "agree" and "disagree"
k <- which(names(wevs)=="C002")
wevs[which(wevs[,k]==6), k] <- rep(21, length(which(wevs[,k]==6)))
wevs[which(wevs[,k]==7), k] <- rep(23, length(which(wevs[,k]==7)))
wevs[which(wevs[,k]==8), k] <- rep(22, length(which(wevs[,k]==8)))

# 75 C061; [6(53189) to21, 7(61789) to23, 8(42086)to22]
# Putting "depends" between "following instructions" and "must be convinced first"
k <- which(names(wevs)=="C061")
wevs[which(wevs[,k]==6), k] <- rep(21, length(which(wevs[,k]==6)))
wevs[which(wevs[,k]==7), k] <- rep(23, length(which(wevs[,k]==7)))
wevs[which(wevs[,k]==8), k] <- rep(22, length(which(wevs[,k]==8)))

# 80 D023; [6(99836) to21, 7(94672) to23, 8(43574)to22]
# Putting "depends" between "disapprove" and "approve"
k <- which(names(wevs)=="D023")
wevs[which(wevs[,k]==6), k] <- rep(21, length(which(wevs[,k]==6)))
wevs[which(wevs[,k]==7), k] <- rep(23, length(which(wevs[,k]==7)))
wevs[which(wevs[,k]==8), k] <- rep(22, length(which(wevs[,k]==8)))

# 97 E022; [6(73752) to21, 7(22759) to23, 8(41889) to22]
# Putting "some of each" in between "will help" and "will harm" 
k <- which(names(wevs)=="E022")
wevs[which(wevs[,k]==6), k] <- rep(21, length(which(wevs[,k]==6)))
wevs[which(wevs[,k]==7), k] <- rep(23, length(which(wevs[,k]==7)))
wevs[which(wevs[,k]==8), k] <- rep(22, length(which(wevs[,k]==8)))

# 123 F028;   
# 6:13 map onto c(34208,55291,33950,54152,5785,20508,32767,84949)
# need to convert 10s (5785 "other specific holidays") into 9s (54152 "only on special holy days")
# and then need to subtract 1 from each of 11:13, to make intervals even.
# Summary: Merging "other specific holidays" with "Only on special holy days/Christmas/Easter days"
k <- which(names(wevs)=="F028")
wevs[which(wevs[,k]==10), k] <- rep(9, length(which(wevs[,k]==10)))
wevs[which(wevs[,k]==11), k] <- rep(10, length(which(wevs[,k]==11)))
wevs[which(wevs[,k]==12), k] <- rep(11, length(which(wevs[,k]==12)))
wevs[which(wevs[,k]==13), k] <- rep(12, length(which(wevs[,k]==13)))

# 188 G001; [12(23572) to7]
# Putting "Region of country" into position of "region" which is empty for WVS
k <- which(names(wevs)=="G001")
wevs[which(wevs[,k]==12), k] <- rep(7, length(which(wevs[,k]==12)))

# 189 G001; [12(44635) to7]
# Putting "Region of country" into position of "region" which is empty for WVS
#k <- which(names(wevs)=="G002")
#wevs[which(wevs[,k]==12), k] <- rep(7, length(which(wevs[,k]==12)))

wevs.copy <- wevs

# ===== Min-max rescaling =====
for (i in (N.extravars+1):dim(wevs)[2] ){  
  wevs[,i] <- as.integer(wevs[,i]) # first convert variable from factor to integer
  wevs.norm <- wevs[,i]
  wevs.norm[which(wevs.norm<6)] <- NA # all negative-valued responses get turned into NA
  wevs[,i] <- (wevs.norm - min(wevs.norm, na.rm=T))/ (max(wevs.norm, na.rm=T) - min(wevs.norm, na.rm=T)) # perform min-max normalisation
}
pops.order <- order(as.character(unique(wevs$pop)))
pops <- unique(wevs$pop)[pops.order]

# Get modified ('non-allelised')  Muthukrishna cultural variance
cultvar.matrix.wevs <- matrix(NA, length(pops), dim(wevs)[2]-N.extravars) #rows are populations, columns are wvs/evs questions
for (col in 1:dim(cultvar.matrix.wevs)[2]){
  print(paste(col,"of",dim(cultvar.matrix.wevs)[2]))
  cultvar.matrix.wevs[,col] <-  sapply(1:length(pops), function(y){ var(wevs[which(wevs$pop==pops[y]), col+N.extravars ], na.rm=T)})
  }
cultvar.means.wevs <- rowMeans(cultvar.matrix.wevs, na.rm=T)
var.societies.wevs <- data.frame(pops, "CountryISO"= rep(NA, length(cultvar.means.wevs)) , cultvar.means.wevs) 
var.societies.wevs <- var.societies.wevs[order(var.societies.wevs$pops),]
levels(var.societies.wevs$pops)[which(levels(var.societies.wevs$pops)=="Great Britain_(826)")] <- "Great Britain"
# attaching the ISO codes only to countries for which we have heritability data 
var.countrymatch <-  match(relevant.countries, var.societies.wevs$pops)
iso.array <- relavant.iso[!is.na(var.countrymatch)] # NA in var.countrymatch means that the country is in MaTCH but not in WVS, so remove
var.countrymatch <- var.countrymatch[!is.na(var.countrymatch)] 
var.societies.wevs$CountryISO[var.countrymatch] <- iso.array
var.societies.wevs$CountryISO <- as.factor(var.societies.wevs$CountryISO)
#var.societies.reduced.wevs <- var.societies[var.societies$CountryISO %in% country,] # only inlcude countries for which we have heritability data


# °·°·°·°·°·°·°· Cultural variance for WVS only data [repeated from above for WVS only] °·°·°·°·°·°·°·  

positiveresponses.wvsonly <- rep(NA,dim(wvsonly)[2])
for (i in (N.extravars+1):dim(wvsonly)[2]){ positiveresponses.wvsonly[i] <- sum(wvsonly[,i]>5, na.rm=T) } 
wvsonly.noresponses <- names(wvsonly)[which(positiveresponses.wvsonly %in% 0)] # variables with no positive responses (i.e., >5 in this data, or positive number in WVS codebook)

# NOMINAL VARIABLES in the set of 145 WVS variables——50 B008;  55 C009; 74 C061; 88 E001; 89 E003; 90 E005; 122 F022; 
# (A025 & A026 on the verge of nominal, depends on whether we treat "neither" as in-between the two other responses. Treating here as ordinal.)
wevs145.nominalvars <- c("A025","A026","B008","C009","C060","E001","E003","E005","F022")
#wevs191.nominalvars <- c("A025","A026","B008","C009","C010","C060","E001","E002","E003","E004","E005","E006","E135","E136","E137","E138","E139","F022", "X007")
wvsonly.removal <- union(wevs145.nominalvars, wvsonly.noresponses)
wvsonly <- wvsonly[,-which(names(wvsonly) %in% wvsonly.removal)] # remove columns of nominal and empty columns

# ===== Recoding of factors to make ordinal (necessary for rescaling of variables) =====

# 52 C001; [6(98446) to21, 7(164408) to23, 8(42914) to 22, ]
# Putting "neither" in between "agree" and "disagree"
k <- which(names(wvsonly)=="C001")
wvsonly[which(wvsonly[,k]==6), k] <- rep(21, length(which(wvsonly[,k]==6)))
wvsonly[which(wvsonly[,k]==7), k] <- rep(23, length(which(wvsonly[,k]==7)))
wvsonly[which(wvsonly[,k]==8), k] <- rep(22, length(which(wvsonly[,k]==8)))

# 53 C002; [6(203625) to21, 7(60178) to23, 8(32795) to22, ]
# Putting "neither" in between "agree" and "disagree"
k <- which(names(wvsonly)=="C002")
wvsonly[which(wvsonly[,k]==6), k] <- rep(21, length(which(wvsonly[,k]==6)))
wvsonly[which(wvsonly[,k]==7), k] <- rep(23, length(which(wvsonly[,k]==7)))
wvsonly[which(wvsonly[,k]==8), k] <- rep(22, length(which(wvsonly[,k]==8)))

# 75 C061; [6(53189) to21, 7(61789) to23, 8(42086)to22]
# Putting "depends" between "following instructions" and "must be convinced first"
k <- which(names(wvsonly)=="C061")
wvsonly[which(wvsonly[,k]==6), k] <- rep(21, length(which(wvsonly[,k]==6)))
wvsonly[which(wvsonly[,k]==7), k] <- rep(23, length(which(wvsonly[,k]==7)))
wvsonly[which(wvsonly[,k]==8), k] <- rep(22, length(which(wvsonly[,k]==8)))

# 80 D023; [6(99836) to21, 7(94672) to23, 8(43574)to22]
# Putting "depends" between "disapprove" and "approve"
k <- which(names(wvsonly)=="D023")
wvsonly[which(wvsonly[,k]==6), k] <- rep(21, length(which(wvsonly[,k]==6)))
wvsonly[which(wvsonly[,k]==7), k] <- rep(23, length(which(wvsonly[,k]==7)))
wvsonly[which(wvsonly[,k]==8), k] <- rep(22, length(which(wvsonly[,k]==8)))

# 97 E022; [6(73752) to21, 7(22759) to23, 8(41889) to22]
# Putting "some of each" in between "will help" and "will harm" 
k <- which(names(wvsonly)=="E022")
wvsonly[which(wvsonly[,k]==6), k] <- rep(21, length(which(wvsonly[,k]==6)))
wvsonly[which(wvsonly[,k]==7), k] <- rep(23, length(which(wvsonly[,k]==7)))
wvsonly[which(wvsonly[,k]==8), k] <- rep(22, length(which(wvsonly[,k]==8)))

# 123 F028;   
# 6:13 map onto c(34208,55291,33950,54152,5785,20508,32767,84949)
# need to convert 10s (5785 "other specific holidays") into 9s (54152 "only on special holy days")
# and then need to subtract 1 from each of 11:13, to make intervals even.
# Summary: Merging "other specific holidays" with "Only on special holy days/Christmas/Easter days"
k <- which(names(wvsonly)=="F028")
wvsonly[which(wvsonly[,k]==10), k] <- rep(9, length(which(wvsonly[,k]==10)))
wvsonly[which(wvsonly[,k]==11), k] <- rep(10, length(which(wvsonly[,k]==11)))
wvsonly[which(wvsonly[,k]==12), k] <- rep(11, length(which(wvsonly[,k]==12)))
wvsonly[which(wvsonly[,k]==13), k] <- rep(12, length(which(wvsonly[,k]==13)))

# 188 G001; [12(23572) to7]
# Putting "Region of country" into position of "region" which is empty for WVS
k <- which(names(wvsonly)=="G001")
wvsonly[which(wvsonly[,k]==12), k] <- rep(7, length(which(wvsonly[,k]==12)))

# 189 G001; [12(44635) to7]
# Putting "Region of country" into position of "region" which is empty for WVS
#k <- which(names(wvsonly)=="G002")
#wvsonly[which(wvsonly[,k]==12), k] <- rep(7, length(which(wvsonly[,k]==12)))

wvsonly.copy <- wvsonly

# ===== Min-max rescaling =====
for (i in (N.extravars+1):dim(wvsonly)[2] ){  
  wvsonly[,i] <- as.integer(wvsonly[,i]) # first convert variable from factor to integer
  wvsonly.norm <- wvsonly[,i]
  wvsonly.norm[which(wvsonly.norm<6)] <- NA # all negative-valued responses get turned into NA
  wvsonly[,i] <- (wvsonly.norm - min(wvsonly.norm, na.rm=T))/ (max(wvsonly.norm, na.rm=T) - min(wvsonly.norm, na.rm=T)) # perform min-max normalisation
}
pops.order <- order(as.character(unique(wvsonly$pop)))
pops <- unique(wvsonly$pop)[pops.order]

# Get modified ('non-allelised')  Muthukrishna cultural variance
cultvar.matrix.wvsonly <- matrix(NA, length(pops), dim(wvsonly)[2]-N.extravars) #rows are populations, columns are wvs/evs questions
for (col in 1:dim(cultvar.matrix.wvsonly)[2]){
  print(paste(col,"of",dim(cultvar.matrix.wvsonly)[2]))
  cultvar.matrix.wvsonly[,col] <-  sapply(1:length(pops), function(y){ var(wvsonly[which(wvsonly$pop==pops[y]), col+N.extravars ], na.rm=T)})
}
cultvar.means.wvsonly <- rowMeans(cultvar.matrix.wvsonly, na.rm=T)
var.societies.wvsonly <- data.frame(pops, "CountryISO"= rep(NA, length(cultvar.means.wvsonly)) , cultvar.means.wvsonly) 
var.societies.wvsonly <- var.societies.wvsonly[order(var.societies.wvsonly$pops),]
levels(var.societies.wvsonly$pops)[which(levels(var.societies.wvsonly$pops)=="Great Britain_(826)")] <- "Great Britain"
# attaching the ISO codes only to countries for which we have heritability data 
var.countrymatch <-  match(relevant.countries, var.societies.wvsonly$pops)
iso.array <- relavant.iso[!is.na(var.countrymatch)] # NA in var.countrymatch means that the country is in MaTCH but not in WVS, so remove
var.countrymatch <- var.countrymatch[!is.na(var.countrymatch)] 
var.societies.wvsonly$CountryISO[var.countrymatch] <- iso.array
var.societies.wvsonly$CountryISO <- as.factor(var.societies.wvsonly$CountryISO)
#var.societies.reduced.wvsonly <- var.societies.wvsonly[var.societies.wvsonly$CountryISO %in% country,] # only inlcude countries for which we have heritability data


##°·°·°·°·°·°·°· Combine cultural variance with heritability °·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·

# Merge Muthukrishna and Uz 2015 indices with MaTCH data
dat <- merge(uz[uz$CountryISO %in% country,], dat.pre, by.x="CountryISO", by.y="country", all.y = T)
dat <- data.frame(dat[,1:2], matrix(NA,dim(dat)[1],4), dat$CTL_DS,dat$CTL_DG, dat$CTL_C, dat$CTL_DS, dat$CTL_DG, dat[,-c(1:5)])
colnames(dat)[3:11] <- c("var_muth_wvs","var_muth_wvs.ascribed","var_muth_wevs","var_muth_wevs.ascribed","var_uzDS","var_uzDG","var_uzC", "var_uzDG.ascribed","var_uzDS.ascribed")
dat$var_muth_wvs <- var.societies.wvsonly$cultvar.means[match(dat$CountryISO, var.societies.wvsonly$CountryISO)]
dat$var_muth_wevs <- var.societies.wevs$cultvar.means[match(dat$CountryISO, var.societies.wevs$CountryISO)]

# more wrangling
dat$var_uzDS[which(dat$var_uzDS=="n/a")] <- NA # Convert Uz "n/a"s into proper NAs
dat$var_uzDG[which(dat$var_uzDG=="n/a")] <- NA
dat$var_uzC[which(dat$var_uzC=="n/a")] <- NA
dat$var_uzDS <- as.numeric(as.character(dat$var_uzDS))# CTL values imported as factors so converting into numeric 
dat$var_uzDG <- as.numeric(as.character(dat$var_uzDG))# CTL values imported as factors so converting into numeric 
dat$var_uzC <- as.numeric(as.character(dat$var_uzC))# CTL values imported as factors so converting into numeric 
dat <- dat[order(as.character(dat$CountryISO)),] # reorder by country ISO code

# Which countries are missing from Muthukrishna cultural variance? But many of these are not in WVS for the waves being used anyway. See explanation at 'iso.array' above
dat$var_muth_wvs.ascribed <- dat$var_muth_wvs
dat$var_muth_wevs.ascribed <- dat$var_muth_wevs
dat$var_uzDG.ascribed <- dat$var_uzDG
dat$var_uzDS.ascribed <- dat$var_uzDS

# Manually entering countries to use for ascription of Muthukrishna cultural variance
dat$var_muth_wvs.ascribed[which(dat$CountryISO=="BE")] <- mean(var.societies.wvsonly$cultvar.means[which(var.societies.wvsonly$pops %in% c("France","Germany","Netherlands"))]) # for Belgium
dat$var_muth_wvs.ascribed[which(dat$CountryISO=="DK")] <- mean(var.societies.wvsonly$cultvar.means[which(var.societies.wvsonly$pops %in% c("Sweden","Germany"))]) # for Denmark
dat$var_muth_wvs.ascribed[which(dat$CountryISO=="IL")] <- mean(var.societies.wvsonly$cultvar.means[which(var.societies.wvsonly$pops %in% c("Jordan","Lebanon"))]) # for Israel
#dat$var_muth_wvs.ascribed[which(dat$CountryISO=="LK")] <- mean(var.societies.wvsonly$cultvar.means[which(var.societies.wvsonly$pops %in% c("India"))]) # for Sri Lanka

dat$var_uzDG.ascribed[which(dat$CountryISO=="NO")] <- mean(as.numeric(uz$CTL_DG[which(uz$Country %in% c("Sweden","Finland"))])) # for Norway
dat$var_uzDG.ascribed[which(dat$CountryISO=="IL")] <- mean(as.numeric(uz$CTL_DG[which(uz$Country %in% c("Jordan"))])) # for Israel
#dat$var_uzDG.ascribed[which(dat$CountryISO=="LK")] <- mean(as.numeric(uz$CTL_DG[which(uz$Country %in% c("India"))])) # for Sri Lanka
#dat$var_uzDG.ascribed[which(dat$CountryISO=="AU")] <- mean(as.numeric(uz$CTL_DG[which(uz$Country %in% c("Indonesia"))])) # for Australia
dat$var_uzDG.ascribed[which(dat$CountryISO=="CN")] <- mean(as.numeric(uz$CTL_DG[which(uz$Country %in% c("India","Kyrgyzstan","Russian Federation","VietNam","Pakistan"))])) # for China, only needed for DG
dat$var_uzDG.ascribed[which(dat$CountryISO=="BR")] <- mean(as.numeric(uz$CTL_DG[which(uz$Country %in% c("Argentina", "Peru"))])) # for Brazil, Venezuela missing in UzDG


dat$var_uzDS.ascribed[which(dat$CountryISO=="NO")] <- mean(as.numeric(uz$CTL_DS[which(uz$Country %in% c("Sweden","Finland"))])) # for Norway
dat$var_uzDS.ascribed[which(dat$CountryISO=="IL")] <- mean(as.numeric(uz$CTL_DS[which(uz$Country %in% c("Jordan"))])) # for Israel
dat$var_uzDS.ascribed[which(dat$CountryISO=="BR")] <- mean(as.numeric(uz$CTL_DS[which(uz$Country %in% c("Argentina", "Peru", "Venezuela"))])) # for Brazil, Venezuela missing in UzDG
#dat$var_uzDS.ascribed[which(dat$CountryISO=="LK")] <- mean(as.numeric(uz$CTL_DS[which(uz$Country %in% c("India"))])) # for Sri Lanka
#dat$var_uzDS.ascribed[which(dat$CountryISO=="AU")] <- mean(as.numeric(uz$CTL_DS[which(uz$Country %in% c("Indonesia"))])) # for Australia



#°·°·°·°·°·°·°· Organise V-Dem data °·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·

# subset of V-Dem dataset, for 2006 (median year of all papers included in MaTCH meta-analysis) and for countries remaining from MaTCH data
# Match the countries specified here to the ones in "dat$Country" (i.e., the countries automatically extracted from the MaTCH heritability data) !
vdem.year <- 2006
vdem.sub <- rbind(vdem[which(vdem$country_name == "Australia" & vdem$year == vdem.year),], 
                  vdem[which(vdem$country_name == "Belgium" & vdem$year == vdem.year),], 
                  vdem[which(vdem$country_name == "Brazil" & vdem$year == vdem.year),],                   
                  vdem[which(vdem$country_name == "Canada" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "China" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Czech Republic" & vdem$year == vdem.year),],  
                  vdem[which(vdem$country_name == "Germany" & vdem$year == vdem.year),],                    
                  vdem[which(vdem$country_name == "Denmark" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Spain" & vdem$year == vdem.year),],                  
                  vdem[which(vdem$country_name == "Finland" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "France" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "United Kingdom" & vdem$year == vdem.year),],
                  #                  vdem[which(vdem$country_name == "Gambia" & vdem$year == 2006),],    
                  # manually entering empty vector because Gambia is not in V-Dem dataset
                  c("Gambia", rep(NA, dim(vdem)[2]-1)),
                  vdem[which(vdem$country_name == "Greece" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Croatia" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Hungary" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Israel" & vdem$year == vdem.year),],                   
                  vdem[which(vdem$country_name == "India" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Italy" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Jamaica" & vdem$year == vdem.year),],                   
                  vdem[which(vdem$country_name == "Japan" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "South Korea" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Sri Lanka" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Namibia" & vdem$year == vdem.year),],                   
                  vdem[which(vdem$country_name == "Netherlands" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Norway" & vdem$year == vdem.year),],                   
                  vdem[which(vdem$country_name == "Poland" & vdem$year == vdem.year),],                  
                  vdem[which(vdem$country_name == "Russia" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Sweden" & vdem$year == vdem.year),],
                  vdem[which(vdem$country_name == "Taiwan" & vdem$year == vdem.year),], 
                  vdem[which(vdem$country_name == "United States of America" & vdem$year == vdem.year),]
) 
# further subsetting, taking only variables of interest from full V-Dem dataset
col.select <- c(  
which(colnames(vdem.sub) == "country_name"),
which(colnames(vdem.sub) == "e_migdppc"), #GDP per capita
which(colnames(vdem.sub) == "e_peaveduc"), # years education for citizens older than 15 (missing for Croatia? but Croatia data not in MaTCH anyway)
which(colnames(vdem.sub) == "e_wb_pop"), # population
which(colnames(vdem.sub) == "e_regionpol") # region (politico-geographic). Less sub-divisions than 'geographic region' (e_regiongeo)
)

vdem.selected <- vdem.sub[,col.select]
vdem.selected <- vdem.selected[rep(1:nrow(vdem.selected), each=4),] # multiplying V-Dem rows for each of the four h2 types.
vdem.selected$country_name <- as.factor(vdem.selected$country_name)
vdem.selected$e_migdppc <- as.numeric(vdem.selected$e_migdppc)
vdem.selected$e_peaveduc <- as.numeric(vdem.selected$e_peaveduc)
vdem.selected$e_wb_pop <- as.numeric(vdem.selected$e_wb_pop)
vdem.selected$e_regionpol <- as.numeric(vdem.selected$e_regionpol)


##°·°·°·°·°·°·°· Combine everything °·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·

dat <- cbind(dat, vdem.selected)
dat <- dat[,-which(names(dat)== "CountryISO")] # remove columns for "CountryISO" and "Country", as we use "country_name" in the final data frame
dat <- dat[,-which(names(dat)== "Country")]
colnames(dat)[which(colnames(dat)=="country_name")] <- "country" # rename variable
###dat.long <- melt(dat, na.rm=T, id.vars=c("h2_type","country_name", "var_uzDS", "var_uzDG", "var_uzC", "e_migdppc", "e_peaveduc", "e_wb_pop", "e_regionpol")) # Use melt() to reorganise matrix
# ^ 535 obs, 14 countries, 112 unique subchapters, 
dat.long <- melt(dat, na.rm=T, id.vars=c("h2_type","country", "var_muth_wvs","var_muth_wvs.ascribed","var_muth_wevs","var_muth_wevs.ascribed","var_uzDS","var_uzDG","var_uzC", "var_uzDS.ascribed","var_uzDG.ascribed","e_migdppc", "e_peaveduc", "e_wb_pop", "e_regionpol")) # Use melt() to reorganise matrix
dat.long[,c(dim(dat.long)[2],dim(dat.long)[2]-1)] <- dat.long[,c(dim(dat.long)[2]-1,dim(dat.long)[2])] # simply putting trait name after heritability
colnames(dat.long)[c(dim(dat.long)[2]-1,dim(dat.long)[2])] <- c("heritability","trait_string") # renaming "value" and "variable"

# trait labels in full verbal form (not just strings extracted from file names), from Polderman et al 2015 supplementary table 33. Using these to fetch higher-order categories as well the same table.  
# manually matching: unique(dat.long$trait_string)
supp32traitlabels <-  c("Acquired Absence of Organs, Not Elsewhere Classified","All-Cause Mortality","Allergy, Unspecified","Asthma","Atopic Dermatitis","Attention Functions","Basic Interpersonal Interactions","Bipolar Affective Disorder","Blood Pressure Functions","Blood Vessel Functions","Calculation Functions","Complex Interpersonal Interactions","Conduct Disorders","Coxarthrosis [Arthrosis of Hip]",
"Cystic Fibrosis","Dementia In Alzheimer Disease","Depressive Episode","Diseases of Pulp and Periapical Tissues","Disorders of Social Functioning with Onset Specific to Childhood and Adolescence","Dissocial Personality Disorder","Dorsalgia","Eating Disorders","Education","Emotional Disorders with Onset Specific to Childhood",
"Emotionally Unstable Personality Disorder","Endocrine Gland Functions","Energy and Drive Functions","Exercise Tolerance Functions","Experience of Self and Time Functions","Expression","Food","Function of Brain","Gene Expression","General Metabolic Functions","Gestational [Pregnancy-Induced] Hypertension Without Significant Proteinuria","Global Psychosocial Functions",
"Gout","Habit and Impulse Disorders","Haematological System Functions","Hearing Functions","Heart Functions","Height","Higher-Level Cognitive Functions","Hyperkinetic Disorders","Hyperplasia of Prostate","Immunological System Functions","Individual Attitudes of Strangers","Informal Social Relationships","Intellectual Functions",
"Intimate Relationships","Irritable Bowel Syndrome","Labour and Delivery Complicated By Umbilical Cord Complications","Leiomyoma of Uterus","Looking After One's Health","Memory Functions","Menstruation Functions","Mental and Behavioural Disorders Due to Multiple Drug Use and Use of Other Psychoactive Substances",
"Mental and Behavioural Disorders Due to Use of Alcohol","Mental and Behavioural Disorders Due to Use of Cannabinoids","Mental and Behavioural Disorders Due to Use of Cocaine","Mental and Behavioural Disorders Due to Use of Hallucinogens","Mental and Behavioural Disorders Due to Use of Opioids","Mental and Behavioural Disorders Due to Use of Sedatives Or Hypnotics",
"Mental and Behavioural Disorders Due to Use of Tobacco","Mental Functions, Unspecified","Migraine","Mild Mental Retardation","Mood [Affective] Disorders","Mortality From Heart Disease","Muscle Power Functions","Nonsuppurative Otitis Media","Obsessive-Compulsive Disorder","Osteoporosis In Diseases Classified Elsewhere","Other Anxiety Disorders",
"Other Functions of the Skin","Other Nontoxic Goitre","Other Symptoms and Signs Involving the Urinary System","Pain and Other Conditions Associated with Female Genital Organs and Menstrual Cycle","Parkinson Disease","Pervasive Developmental Disorders","Phobic Anxiety Disorders","Potential Health Hazards Related to Socioeconomic and Psychosocial Circumstances",
"Problems Related to Upbringing","Procreation Functions","Protective Functions of the Skin","Psychological and Behavioural Disorders Associated with Sexual Development and Orientation","Psychomotor Functions","Reaction to Severe Stress, and Adjustment Disorders","Recreation and Leisure","Recurrent Depressive Disorder","Religion and Spirituality","Schizophrenia","Schizotypal Disorder","Seeing Functions","Senile Cataract","Sensation of Pain","Sexual Functions","Sleep Disorders",
"Sleep Functions","Slow Fetal Growth and Fetal Malnutrition","Societal Attitudes","Specific Developmental Disorder of Motor Function","Specific Personality Disorders","Structure of Areas of Skin","Structure of Brain","Structure of Cardiovascular System","Structure of Eyeball","Structure of Head and Neck Region","Structure of Lower Extremity","Structure of Mouth","Structure of Pelvic Region","Structure of Trunk",
"Superficial Injuries Involving Multiple Body Regions","Temperament and Personality Functions","Vasomotor and Allergic Rhinitis","Vestibular Functions","Voice Functions","Walking and Moving","Water, Mineral and Electrolyte Balance Functions",
"Weight Maintenance Functions")

# Creating a table that lists traits across hierarchical levels  
colnames(trait_hierarchy) <- c("domain","chapter","subchapter")
trait_hierarchy <- unique(trait_hierarchy)
trait_hierarchy <- trait_hierarchy[order(trait_hierarchy$domain) , ]
trait_hierarchy <- trait_hierarchy[order(trait_hierarchy$chapter) , ]
trait_hierarchy <- trait_hierarchy[order(trait_hierarchy$subchapter) , ]
# ^ Duplicates: Basic Interpersonal Interactions; Exercise Tolerance Functions; Immunological System Functions
# these traits appear more than once because they fall under more than one super-category, but in the code below, we select only one instance of each 

# Aligning subchapter level traits in above table with strings 
trait_alignment <-  data.frame("string"= unique(dat.long$trait_string), "verbal"=supp32traitlabels )
dat.long$trait_subchapter <- trait_alignment$verbal[match(dat.long$trait_string, trait_alignment$string)]
dat.long$trait_chapter <- trait_hierarchy$chapter[match(dat.long$trait_subchapter, trait_hierarchy$subchapter)]
dat.long$trait_domain <- trait_hierarchy$domain[match(dat.long$trait_subchapter, trait_hierarchy$subchapter)]

# Individually specifying list of traits within dat.long that we consider as being culturally transmissable
culturaltraits_nonpsychiatric <-  which(dat.long$trait_string %in% c(
  "attention_functions", "basic_interpersonal_interactions","calculation_functions","complex_interpersonal_interactions","education","global_psychosocial_functions","higher.level_cognitive_functions","individual_attitudes_of_strangers","informal_social_relationships","intellectual_functions",
  "looking_after_ones_health","memory_functions","mild_mental_retardation","potential_health_hazards_related_to_socioeconomic_and_psychosocial_circumstances","problems_related_to_upbringing","psychomotor_functions","religion_and_spirituality","societal_attitudes")
  )
dat.long.unscaled <- dat.long # save an unscaled version

# Standardize and transform the variables, although not heritability
dat.long$var_muth_wvs <- scale(dat.long$var_muth_wvs)
dat.long$var_muth_wevs <- scale(dat.long$var_muth_wevs)
dat.long$var_muth_wvs.ascribed <- scale(dat.long$var_muth_wvs.ascribed)
dat.long$var_muth_wevs.ascribed <- scale(dat.long$var_muth_wevs.ascribed)
dat.long$var_uzDS <- scale(dat.long$var_uzDS)
dat.long$var_uzDG <- scale(dat.long$var_uzDG)
dat.long$var_uzC <- scale(dat.long$var_uzC)
dat.long$var_uzDS.ascribed <- scale(dat.long$var_uzDS.ascribed)
dat.long$var_uzDG.ascribed <- scale(dat.long$var_uzDG.ascribed)
dat.long$e_migdppc <- scale(dat.long$e_migdppc) 
dat.long$e_peaveduc <- scale(dat.long$e_peaveduc) 
dat.long$e_wb_pop <- scale(log(dat.long$e_wb_pop)) # log transform as well since highly skewed



# Adding to the above list all traits within Psychiatric domain except for 'sleep functions'
culturaltraits_withpsychiatric <- union(culturaltraits_nonpsychiatric, setdiff(which(dat.long$trait_domain=="Psychiatric"), which(dat.long$trait_subchapter=="Sleep Functions")))
dat.core <- dat.long[culturaltraits_nonpsychiatric,]
dat.expanded <- dat.long[culturaltraits_withpsychiatric,]
dat.nonprereg <- dat.long[ -culturaltraits_withpsychiatric, ]
dat.psychiatric <- dat.long[ setdiff(which(dat.long$trait_domain=="Psychiatric"), which(dat.long$trait_subchapter=="Sleep Functions")),]

dat.core.unscaled <- dat.long.unscaled[culturaltraits_nonpsychiatric,]
dat.expanded.unscaled <- dat.long.unscaled[culturaltraits_withpsychiatric,]
dat.nonprereg.unscaled <- dat.long.unscaled[ -culturaltraits_withpsychiatric, ]
dat.psychiatric.unscaled <- dat.long.unscaled[ setdiff(which(dat.long.unscaled$trait_domain=="Psychiatric"), which(dat.long.unscaled$trait_subchapter=="Sleep Functions")),]


# comparing sample sizes between imputed and unimputed data sets
sum(!is.na(dat.core$var_muth_wvs))
sum(!is.na(dat.core$var_muth_wvs.ascribed))
sum(!is.na(dat.core$var_uzDG))
sum(!is.na(dat.core$var_uzDG.ascribed))
sum(!is.na(dat.core$var_uzDS))
sum(!is.na(dat.core$var_uzDS.ascribed))

sum(!is.na(dat.expanded$var_muth_wvs))
sum(!is.na(dat.expanded$var_muth_wvs.ascribed))
sum(!is.na(dat.expanded$var_uzDG))
sum(!is.na(dat.expanded$var_uzDG.ascribed))
sum(!is.na(dat.expanded$var_uzDS))
sum(!is.na(dat.expanded$var_uzDS.ascribed))

sum(!is.na(dat.nonprereg$var_muth_wvs))
sum(!is.na(dat.nonprereg$var_muth_wvs.ascribed))
sum(!is.na(dat.nonprereg$var_uzDG))
sum(!is.na(dat.nonprereg$var_uzDG.ascribed))
sum(!is.na(dat.nonprereg$var_uzDS))
sum(!is.na(dat.nonprereg$var_uzDS.ascribed))


#°·°·°·°·°·°·°· Output results °·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·°·

# >>> Various tests (not including regression models) are ready to run in "heritability_analysis_results_0x.Rmd"

# Output for use in "heritability_analysis_results_0x.Rmd"
write.csv(dat.long, "heritability-dat-long.csv")
write.csv(dat.core, "heritability-dat-core.csv")
write.csv(dat.expanded, "heritability-dat-expanded.csv")
write.csv(dat.nonprereg, "heritability-dat-nonprereg.csv")
write.csv(dat.psychiatric, "heritability-dat-psychiatric.csv")
write.csv(var.societies.wvsonly, "heritability-var-societies-wvs.csv")
write.csv(var.societies.wevs, "heritability-var-societies-wevs.csv")

# convergence issues. See lme4::?convergence
###strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))

model.1 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.2 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.3 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.4 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.5 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.6 <- lmer(heritability ~ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.7 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.8 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.9 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.10 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.11 <- lmer(heritability ~ var_muth_wvs+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.12 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.13 <- lmer(heritability ~ var_muth_wvs+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.14 <- lmer(heritability ~ var_muth_wvs+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.15 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.16 <- lmer(heritability ~ var_muth_wvs+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.17 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.18 <- lmer(heritability ~ var_muth_wvs+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.19 <- lmer(heritability ~ var_muth_wvs+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.20 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.21 <- lmer(heritability ~ var_muth_wvs.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.22 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.23 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.24 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.25 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.26 <- lmer(heritability ~ var_muth_wvs.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.27 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.28 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.29 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.30 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.31 <- lmer(heritability ~ var_uzDG.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.32 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.33 <- lmer(heritability ~ var_uzDG.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.34 <- lmer(heritability ~ var_uzDG.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.35 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.36 <- lmer(heritability ~ var_uzDG.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.37 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.38 <- lmer(heritability ~ var_uzDG.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.39 <- lmer(heritability ~ var_uzDG.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.40 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.41 <- lmer(heritability ~ var_uzDS.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.42 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.43 <- lmer(heritability ~ var_uzDS.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.44 <- lmer(heritability ~ var_uzDS.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.45 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.46 <- lmer(heritability ~ var_uzDS.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.47 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.48 <- lmer(heritability ~ var_uzDS.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.49 <- lmer(heritability ~ var_uzDS.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.50 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.51 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.52 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.53 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.54 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.55 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.56 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.57 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.58 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.59 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.60 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.61 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.62 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.63 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.64 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.65 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.66 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.67 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.68 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.69 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.70 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.71 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.72 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.73 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.74 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.75 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.76 <- lmer(heritability ~ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.77 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.78 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.79 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.80 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.81 <- lmer(heritability ~ var_muth_wvs+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.82 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.83 <- lmer(heritability ~ var_muth_wvs+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.84 <- lmer(heritability ~ var_muth_wvs+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.85 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.86 <- lmer(heritability ~ var_muth_wvs+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.87 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.88 <- lmer(heritability ~ var_muth_wvs+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.89 <- lmer(heritability ~ var_muth_wvs+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.90 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.91 <- lmer(heritability ~ var_muth_wvs.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.92 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.93 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.94 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.95 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.96 <- lmer(heritability ~ var_muth_wvs.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.97 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.98 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.99 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.100 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.101 <- lmer(heritability ~ var_uzDG.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.102 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.103 <- lmer(heritability ~ var_uzDG.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.104 <- lmer(heritability ~ var_uzDG.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.105 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.106 <- lmer(heritability ~ var_uzDG.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.107 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.108 <- lmer(heritability ~ var_uzDG.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.109 <- lmer(heritability ~ var_uzDG.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.110 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.111 <- lmer(heritability ~ var_uzDS.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.112 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.113 <- lmer(heritability ~ var_uzDS.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.114 <- lmer(heritability ~ var_uzDS.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.115 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.116 <- lmer(heritability ~ var_uzDS.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.117 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.118 <- lmer(heritability ~ var_uzDS.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.119 <- lmer(heritability ~ var_uzDS.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.120 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.121 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.122 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.123 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.124 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.125 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.126 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.127 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.128 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.129 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.130 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.131 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.132 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.133 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.134 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.135 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.136 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.137 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.138 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.139 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.140 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.141 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.142 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.143 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.144 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.145 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.146 <- lmer(heritability ~ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.147 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.148 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.149 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.150 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)

model.151 <- lmer(heritability ~ var_muth_wvs+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.152 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.153 <- lmer(heritability ~ var_muth_wvs+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.154 <- lmer(heritability ~ var_muth_wvs+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.155 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.156 <- lmer(heritability ~ var_muth_wvs+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.157 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.158 <- lmer(heritability ~ var_muth_wvs+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.159 <- lmer(heritability ~ var_muth_wvs+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.160 <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)

model.161 <- lmer(heritability ~ var_muth_wvs.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.162 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.163 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.164 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.165 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.166 <- lmer(heritability ~ var_muth_wvs.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.167 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.168 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.169 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.170 <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)

model.171 <- lmer(heritability ~ var_uzDG.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.172 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.173 <- lmer(heritability ~ var_uzDG.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.174 <- lmer(heritability ~ var_uzDG.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.175 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.176 <- lmer(heritability ~ var_uzDG.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.177 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.178 <- lmer(heritability ~ var_uzDG.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.179 <- lmer(heritability ~ var_uzDG.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.180 <- lmer(heritability ~ var_uzDG.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)

model.181 <- lmer(heritability ~ var_uzDS.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.182 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.183 <- lmer(heritability ~ var_uzDS.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.184 <- lmer(heritability ~ var_uzDS.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.185 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.186 <- lmer(heritability ~ var_uzDS.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.187 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.188 <- lmer(heritability ~ var_uzDS.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.189 <- lmer(heritability ~ var_uzDS.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.190 <- lmer(heritability ~ var_uzDS.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)

model.191 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.192 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.193 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.194 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.195 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.196 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.197 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.198 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.199 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.200 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)

model.201 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.202 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.203 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.204 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.205 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.nonprereg)
model.206 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.207 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.208 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.209 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)
model.210 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.nonprereg)




#:::::::::::::::::::::::::::
#getting mean standard errors across models, as a representation of how much noise in each aggregate analysis
SE.core <- matrix(NA,60,2)
SE.expanded <- matrix(NA,60,2)
SE.nonprereg <- matrix(NA,60,2)
for (i in 11:70){
  ci <- confint(get(paste0("model.",as.character(i))))
  SE.core[i-10,1] <- ci[str_detect(names(ci[,1]),"var"),1]
  SE.core[i-10,2] <- ci[str_detect(names(ci[,1]),"var"),2]
  }
for (i in 81:140){
  ci <- confint(get(paste0("model.",as.character(i))))
  SE.expanded[i-80,1] <- ci[str_detect(names(ci[,1]),"var"),1]
  SE.expanded[i-80,2] <- ci[str_detect(names(ci[,1]),"var"),2]
}
for (i in 151:210){
  ci <- confint(get(paste0("model.",as.character(i))))
  SE.nonprereg[i-150,1] <- ci[str_detect(names(ci[,1]),"var"),1]
  SE.nonprereg[i-150,2] <- ci[str_detect(names(ci[,1]),"var"),2]
}

mean(abs(SE.core[,1]-SE.core[,2]))
mean(abs(SE.expanded[,1]-SE.expanded[,2]))
mean(abs(SE.nonprereg[,1]-SE.nonprereg[,2]))

model.1 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)

model.1 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)

# ---- >>>correlations<<<-----

cor1 <- cor.test(dat.core$heritability, dat.core$var_muth_wvs)
cor2 <- cor.test(dat.core$heritability, dat.core$var_uzDG)
cor3 <- cor.test(dat.core$heritability, dat.core$var_uzDS)
cor4 <- cor.test(dat.expanded$heritability, dat.expanded$var_muth_wvs)
cor5 <- cor.test(dat.expanded$heritability, dat.expanded$var_uzDG)
cor6 <- cor.test(dat.expanded$heritability, dat.expanded$var_uzDS)
cor7 <- cor.test(dat.nonprereg$heritability, dat.nonprereg$var_muth_wvs)
cor8 <- cor.test(dat.nonprereg$heritability, dat.nonprereg$var_uzDG)
cor9 <- cor.test(dat.nonprereg$heritability, dat.nonprereg$var_uzDS)


corN1 <- cor.test(dat.core[dat.core$var_muth_wvs<2 & dat.core$var_muth_wvs>-2,]$heritability,
                  dat.core[dat.core$var_muth_wvs<2 & dat.core$var_muth_wvs>-2,]$var_muth_wvs)
corN2 <- cor.test(dat.core[dat.core$var_uzDG<2 & dat.core$var_uzDG>-2,]$heritability,
                  dat.core[dat.core$var_uzDG<2 & dat.core$var_uzDG>-2,]$var_uzDG)
corN3 <- cor.test(dat.core[dat.core$var_uzDS<2 & dat.core$var_uzDS>-2,]$heritability,
                  dat.core[dat.core$var_uzDS<2 & dat.core$var_uzDS>-2,]$var_uzDS)

corN4 <- cor.test(dat.expanded[dat.expanded$var_muth_wvs<2 & dat.expanded$var_muth_wvs>-2,]$heritability,
                  dat.expanded[dat.expanded$var_muth_wvs<2 & dat.expanded$var_muth_wvs>-2,]$var_muth_wvs)
corN5 <- cor.test(dat.expanded[dat.expanded$var_uzDG<2 & dat.expanded$var_uzDG>-2,]$heritability,
                  dat.expanded[dat.expanded$var_uzDG<2 & dat.expanded$var_uzDG>-2,]$var_uzDG)
corN6 <- cor.test(dat.expanded[dat.expanded$var_uzDS<2 & dat.expanded$var_uzDS>-2,]$heritability,
                  dat.expanded[dat.expanded$var_uzDS<2 & dat.expanded$var_uzDS>-2,]$var_uzDS)

corN7 <- cor.test(dat.nonprereg[dat.nonprereg$var_muth_wvs<2 & dat.nonprereg$var_muth_wvs>-2,]$heritability,
                  dat.nonprereg[dat.nonprereg$var_muth_wvs<2 & dat.nonprereg$var_muth_wvs>-2,]$var_muth_wvs)
corN8 <- cor.test(dat.nonprereg[dat.nonprereg$var_uzDG<2 & dat.nonprereg$var_uzDG>-2,]$heritability,
                  dat.nonprereg[dat.nonprereg$var_uzDG<2 & dat.nonprereg$var_uzDG>-2,]$var_uzDG)
corN9 <- cor.test(dat.nonprereg[dat.nonprereg$var_uzDS<2 & dat.nonprereg$var_uzDS>-2,]$heritability,
                  dat.nonprereg[dat.nonprereg$var_uzDS<2 & dat.nonprereg$var_uzDS>-2,]$var_uzDS)

# correlation plots

cor_expanded_muth <- ggplot(dat.expanded, aes(var_muth_wvs,heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="orangered2") + geom_smooth(method='lm',formula= y~x, se=F)+
  ggtitle("Muthukrishna index")+ xlab("") + ylab("Inclusive set of traits") + xlim(-3.3,2)+ ylim(0,1)+
  annotate("text", x = -3, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor4$estimate,3))))+
  theme(plot.title = element_text(hjust = 0.5, size=13, face = "bold"), axis.title.y = element_text(size = 13, face = "bold"))
cor_expanded_uzDG <- ggplot(dat.expanded, aes(var_uzDG, heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="orangered2")+ geom_smooth(method='lm',formula= y~x, se=F)+
  ggtitle("Uz DG index")+ xlab("") + ylab("") + xlim(-2,3.5)+ ylim(0,1) +
  annotate("text", x = -1.7, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor5$estimate,3))))+
  theme(plot.title = element_text(hjust = 0.5, size=13, face = "bold"))
cor_expanded_uzDS <- ggplot(dat.expanded, aes(var_uzDS, heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="orangered2")+  geom_smooth(method='lm',formula= y~x, se=F)+
  ggtitle("Uz DS index")+ xlab("") + xlab("") + ylab("") + xlim(-3.5,2.5)+ ylim(0,1)+
  annotate("text", x = -3.2, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor6$estimate,3))))+
  theme(plot.title = element_text(hjust = 0.5, size=13, face = "bold"))

cor_core_muth <- ggplot(dat.core, aes(var_muth_wvs,heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="purple2") + geom_smooth(method='lm',formula= y~x, se=F)+
  xlab("") + ylab("Restricted set of traits") + xlim(-3.3,2)+ ylim(0,1)+
  annotate("text", x = -3, y = 0.94, hjust=0,  color="blue", label = paste0("r = ",as.character(round(cor1$estimate,3))))+
  theme(axis.title.y = element_text(size = 13, face = "bold"))
cor_core_uzDG <- ggplot(dat.core, aes(var_uzDG, heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="purple2")+ geom_smooth(method='lm',formula= y~x, se=F)+
  xlab("") + ylab("") + xlim(-2,3.5)+ ylim(0,1)+
  annotate("text", x = -1.7, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor2$estimate,3))))
cor_core_uzDS <- ggplot(dat.core, aes(var_uzDS, heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="purple2")+  geom_smooth(method='lm',formula= y~x, se=F)+
  ylab("") + xlim(-3.5,2.5)+ ylim(0,1)+
  annotate("text", x = -3.2, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor3$estimate,3))))

cor_nonprereg_muth <- ggplot(dat.nonprereg, aes(var_muth_wvs,heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="springgreen4") + geom_smooth(method='lm', formula= y~x, se=F)+
  xlab("") + ylab("Acultural set of traits") + xlim(-3.3,2)+ ylim(0,1)+
  annotate("text", x = -3, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor7$estimate,3))))+
  theme(axis.title.y = element_text(size = 13, face = "bold"))
cor_nonprereg_uzDG <- ggplot(dat.nonprereg, aes(var_uzDG, heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="springgreen4")+ geom_smooth(method='lm', formula= y~x, se=F)+
  xlab("") + ylab("") + xlim(-2,3.5)+ ylim(0,1)+
  annotate("text", x = -1.7, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor8$estimate,3))))
cor_nonprereg_uzDS <- ggplot(dat.nonprereg, aes(var_uzDS, heritability))+
  theme_bw()+ geom_point(alpha=0.5, color="springgreen4")+  geom_smooth(method='lm', formula= y~x, se=F)+
  xlab("") + ylab("") + xlim(-3.5,2.5)+ ylim(0,1)+
  annotate("text", x = -3.2, y = 0.94, hjust=0, color="blue", label = paste0("r = ",as.character(round(cor9$estimate,3))))
  


correlations_agg <-  cowplot::plot_grid(cor_expanded_muth, cor_expanded_uzDG, cor_expanded_uzDS,
                                        cor_core_muth, cor_core_uzDG, cor_core_uzDS,
                                    cor_nonprereg_muth, cor_nonprereg_uzDG, cor_nonprereg_uzDS,
                                    ncol = 3, nrow=3, align = "hv")
y_title <- ggdraw()+ draw_label("Reported heritability", fontface = 'bold', vjust = 0.5, angle=90, size=16) + theme(plot.margin = margin(0, 0, 0, 0))
#x_title <- ggdraw()+ draw_label("Cultural Variance", fontface = 'bold', hjust=1, size=16)  + theme(plot.margin = margin(0, 0, 0, 0))

correlations_agg <-  cowplot::plot_grid(y_title, correlations_agg, nrow=1, rel_widths = c(0.1, 1))
correlations_agg <- add_sub(correlations_agg, "Cultural variance",x=0.57,fontface = 'bold',size=16)
ggsave("correlations_agg_color8.png", correlations_agg, width=200, height=200, units="mm")

# ===============





#---------
#---------
# main results
# ~~~~~~~~~~~~~

results.core <- data.frame(matrix(NA, 30, 14))
colnames(results.core) <- c("cultural.variance", "controls", "survey", "countries.ascribed", "trait.level", "specification","estimate", "lowerCI", "upperCI", "p", "slopesign", "intercept", "marginal.R2", "conditional.R2")
results.core[,1] <- c(rep("Muthukrishna",10), rep("Uz DG",10), rep("Uz DS",10))
results.core[,2] <- rep(c("none", "GDPpc", "education", "population", "all three"), 6)
results.core[,3] <- c(rep("WVS",10), rep("WVS & EVS", 20))
results.core[,4] <- rep("not ascribed",30)
results.core[,5] <- rep(c(rep("subchapter",5), rep("domain",5)),3)
options(warn=0) 
for (i in 1:30){
  print(i)
  results.core[i,6] <- i
  results.core[i,7] <- summary(get(paste0("model.",as.character(c(11:20,51:70)[i]))))$coefficients[2,1]
  ci <- confint(get(paste0("model.",as.character(c(11:20,51:70)[i]))))
  results.core[i,8] <- ci[str_detect(names(ci[,1]),"var"),1]
  results.core[i,9] <- ci[str_detect(names(ci[,1]),"var"),2]
  results.core[i,10] <- summary(get(paste0("model.",as.character(c(11:20,51:70)[i]))))$coefficients[2,5]
  results.core[i,11] <- ifelse(results.core[i,7] >0, "+", "-")
  results.core[i,12] <- summary(get(paste0("model.",as.character(c(11:20,51:70)[i]))))$coefficients[1,1]
  results.core[i,13] <- r.squaredGLMM(get(paste0("model.",as.character(c(11:20,51:70)[i]))))[1]
  results.core[i,14] <- r.squaredGLMM(get(paste0("model.",as.character(c(11:20,51:70)[i]))))[2]
}  
results.core$specification <- 1:30
results.core$specification <- factor(results.core$specification, levels=unique(as.character(results.core$specification)))

results.expanded <- data.frame(matrix(NA, 30, 14))
colnames(results.expanded) <- c("cultural.variance", "controls", "survey", "countries.ascribed", "trait.level", "specification","estimate", "lowerCI", "upperCI", "p", "slopesign", "intercept", "marginal.R2", "conditional.R2")
results.expanded[,1] <- c(rep("Muthukrishna",10), rep("Uz DG",10), rep("Uz DS",10))
results.expanded[,2] <- rep(c("none", "GDPpc", "education", "population", "all three"), 6)
results.expanded[,3] <- c(rep("WVS",10), rep("WVS & EVS", 20))
results.expanded[,4] <- rep("not ascribed",30)
results.expanded[,5] <- rep(c(rep("subchapter",5), rep("domain",5)),3)
for (i in 1:30){
  print(i)
  results.expanded[i,6] <- i
  results.expanded[i,7] <- summary(get(paste0("model.",as.character(c(81:90,121:140)[i]))))$coefficients[2,1]
  ci <- confint(get(paste0("model.",as.character(c(81:90,121:140)[i]))))
  results.expanded[i,8] <- ci[str_detect(names(ci[,1]),"var"),1]
  results.expanded[i,9] <- ci[str_detect(names(ci[,1]),"var"),2]
  results.expanded[i,10] <- summary(get(paste0("model.",as.character(c(81:90,121:140)[i]))))$coefficients[2,5]
  results.expanded[i,11] <- ifelse(results.expanded[i,7] > 0, "+", "-")
  results.expanded[i,12] <- summary(get(paste0("model.",as.character(c(81:90,121:140)[i]))))$coefficients[1,1]
  results.expanded[i,13] <- r.squaredGLMM(get(paste0("model.",as.character(c(81:90,121:140)[i]))))[1]
  results.expanded[i,14] <- r.squaredGLMM(get(paste0("model.",as.character(c(81:90,121:140)[i]))))[2]
  }  
results.expanded$specification <- 1:30
results.expanded$specification <- factor(results.expanded$specification, levels=unique(as.character(results.expanded$specification)))


results.nonprereg <- data.frame(matrix(NA, 30, 14))
colnames(results.nonprereg) <- c("cultural.variance", "controls", "survey", "countries.ascribed", "trait.level", "specification","estimate", "lowerCI", "upperCI", "p", "slopesign", "intercept", "marginal.R2", "conditional.R2")
results.nonprereg[,1] <- c(rep("Muthukrishna",10), rep("Uz DG",10), rep("Uz DS",10))
results.nonprereg[,2] <- rep(c("none", "GDPpc", "education", "population", "all three"), 6)
results.nonprereg[,3] <- c(rep("WVS",10), rep("WVS & EVS", 20))
results.nonprereg[,4] <- rep("not ascribed",30)
results.nonprereg[,5] <- rep(c(rep("subchapter",5), rep("domain",5)),3)
for (i in 1:30){
  print(i)
  results.nonprereg[i,6] <- i
  results.nonprereg[i,7] <- summary(get(paste0("model.",as.character(c(151:160,191:210)[i]))))$coefficients[2,1]
  ci <- confint(get(paste0("model.",as.character(c(151:160,191:210)[i]))))
  results.nonprereg[i,8] <- ci[str_detect(names(ci[,1]),"var"),1]
  results.nonprereg[i,9] <- ci[str_detect(names(ci[,1]),"var"),2]
  results.nonprereg[i,10] <- summary(get(paste0("model.",as.character(c(151:160,191:210)[i]))))$coefficients[2,5]
  results.nonprereg[i,11] <- ifelse(results.nonprereg[i,7] > 0, "+", "-")
  results.nonprereg[i,12] <- summary(get(paste0("model.",as.character(c(151:160,191:210)[i]))))$coefficients[1,1]
  results.nonprereg[i,13] <- r.squaredGLMM(get(paste0("model.",as.character(c(151:160,191:210)[i]))))[1]
  results.nonprereg[i,14] <- r.squaredGLMM(get(paste0("model.",as.character(c(151:160,191:210)[i]))))[2]
  }  
results.nonprereg$specification <- 1:30
results.nonprereg$specification <- factor(results.nonprereg$specification, levels=unique(as.character(results.nonprereg$specification)))


table.core <- results.core[c("specification","estimate","lowerCI","upperCI","p","marginal.R2","conditional.R2")]
table.core[,2:7] <- round(table.core[,2:7], 3) 
table.expanded <- results.expanded[c("specification","estimate","lowerCI","upperCI","p","marginal.R2","conditional.R2")]
table.expanded[,2:7] <- round(table.expanded[,2:7], 3) 
table.nonprereg <- results.nonprereg[c("specification","estimate","lowerCI","upperCI","p","marginal.R2","conditional.R2")]
table.nonprereg[,2:7] <- round(table.nonprereg[,2:7], 3) 

write.csv(table.core, "table_core.csv")
write.csv(table.expanded, "table_expanded.csv")
write.csv(table.nonprereg, "table_nonprereg.csv")

mean(table.core$estimate)
mean(table.expanded$estimate)
mean(table.nonprereg$estimate)

# ++++++++//////////
# results, with ascription of missing countries
# ++++++++//////////

results.core.ascribed <- data.frame(matrix(NA, 30, 14))
colnames(results.core.ascribed) <- c("cultural.variance", "controls", "survey", "countries.ascribed", "trait.level", "specification","estimate", "lowerCI", "upperCI", "p", "slopesign", "intercept", "marginal.R2", "conditional.R2")
results.core.ascribed[,1] <- c(rep("Muthukrishna",10), rep("Uz DG",10), rep("Uz DS",10))
results.core.ascribed[,2] <- rep(c("none", "GDPpc", "education", "population", "all three"), 6)
results.core.ascribed[,3] <- c(rep("WVS",10), rep("WVS & EVS", 20))
results.core.ascribed[,4] <- rep("ascribed",30)
results.core.ascribed[,5] <- rep(c(rep("subchapter",5), rep("domain",5)),3)
for (i in 1:30){
  print(i)
  results.core.ascribed[i,6] <- i
  results.core.ascribed[i,7] <- summary(get(paste0("model.",as.character(c(21:50)[i]))))$coefficients[2,1]
  ci <- confint(get(paste0("model.",as.character(c(21:50)[i]))))
  results.core.ascribed[i,8] <- ci[str_detect(names(ci[,1]),"var"),1]
  results.core.ascribed[i,9] <- ci[str_detect(names(ci[,1]),"var"),2]
  results.core.ascribed[i,10] <- summary(get(paste0("model.",as.character(c(21:50)[i]))))$coefficients[2,5]
  results.core.ascribed[i,11] <- ifelse(results.core.ascribed[i,7] > 0, "+", "-")
  results.core.ascribed[i,12] <- summary(get(paste0("model.",as.character(c(21:50)[i]))))[1]
  results.core.ascribed[i,13] <- r.squaredGLMM(get(paste0("model.",as.character(c(21:50)[i]))))[1]
  results.core.ascribed[i,14] <- r.squaredGLMM(get(paste0("model.",as.character(c(21:50)[i]))))[2]
  }  
results.core.ascribed$specification <- 1:30
results.core.ascribed$specification <- factor(results.core.ascribed$specification, levels=unique(as.character(results.core.ascribed$specification)))


results.expanded.ascribed <- data.frame(matrix(NA, 30, 14))
colnames(results.expanded.ascribed) <- c("cultural.variance", "controls", "survey", "countries.ascribed", "trait.level", "specification","estimate", "lowerCI", "upperCI", "p", "slopesign", "intercept", "marginal.R2", "conditional.R2")
results.expanded.ascribed[,1] <- c(rep("Muthukrishna",10), rep("Uz DG",10), rep("Uz DS",10))
results.expanded.ascribed[,2] <- rep(c("none", "GDPpc", "education", "population", "all three"), 6)
results.expanded.ascribed[,3] <- c(rep("WVS",10), rep("WVS & EVS", 20))
results.expanded.ascribed[,4] <- rep("ascribed",30)
results.expanded.ascribed[,5] <- rep(c(rep("subchapter",5), rep("domain",5)),3)
for (i in 1:30){
  print(i)
  results.expanded.ascribed[i,6] <- i
  results.expanded.ascribed[i,7] <- summary(get(paste0("model.",as.character(c(91:120)[i]))))$coefficients[2,1]
  ci <- confint(get(paste0("model.",as.character(c(91:120)[i]))))
  results.expanded.ascribed[i,8] <- ci[str_detect(names(ci[,1]),"var"),1]
  results.expanded.ascribed[i,9] <- ci[str_detect(names(ci[,1]),"var"),2]
  results.expanded.ascribed[i,10] <- summary(get(paste0("model.",as.character(c(91:120)[i]))))$coefficients[2,5]
  results.expanded.ascribed[i,11] <- ifelse(results.expanded.ascribed[i,7] > 0, "+", "-")
  results.expanded.ascribed[i,12] <- summary(get(paste0("model.",as.character(c(91:120)[i]))))$coefficients[1,1]
  results.expanded.ascribed[i,13] <- r.squaredGLMM(get(paste0("model.",as.character(c(91:120)[i]))))[1]
  results.expanded.ascribed[i,14] <- r.squaredGLMM(get(paste0("model.",as.character(c(91:120)[i]))))[2]
}  
results.expanded.ascribed$specification <- 1:30
results.expanded.ascribed$specification <- factor(results.expanded.ascribed$specification, levels=unique(as.character(results.expanded.ascribed$specification)))


results.nonprereg.ascribed <- data.frame(matrix(NA, 30, 14))
colnames(results.nonprereg.ascribed) <- c("cultural.variance", "controls", "survey", "countries.ascribed", "trait.level", "specification","estimate", "lowerCI", "upperCI", "p", "slopesign", "intercept", "marginal.R2", "conditional.R2")
results.nonprereg.ascribed[,1] <- c(rep("Muthukrishna",10), rep("Uz DG",10), rep("Uz DS",10))
results.nonprereg.ascribed[,2] <- rep(c("none", "GDPpc", "education", "population", "all three"), 6)
results.nonprereg.ascribed[,3] <- c(rep("WVS",10), rep("WVS & EVS", 20))
results.nonprereg.ascribed[,4] <- rep("ascribed",30)
results.nonprereg.ascribed[,5] <- rep(c(rep("subchapter",5), rep("domain",5)),3)
for (i in 1:30){
  print(i)
  results.nonprereg.ascribed[i,6] <- i
  results.nonprereg.ascribed[i,7] <- summary(get(paste0("model.",as.character(c(161:190)[i]))))$coefficients[2,1]
  ci <- confint(get(paste0("model.",as.character(c(161:190)[i]))))
  results.nonprereg.ascribed[i,8] <- ci[str_detect(names(ci[,1]),"var"),1]
  results.nonprereg.ascribed[i,9] <- ci[str_detect(names(ci[,1]),"var"),2]
  results.nonprereg.ascribed[i,10] <- summary(get(paste0("model.",as.character(c(161:190)[i]))))$coefficients[2,5]
  results.nonprereg.ascribed[i,11] <- ifelse(results.nonprereg.ascribed[i,7] > 0, "+", "-")
  results.nonprereg.ascribed[i,12] <- summary(get(paste0("model.",as.character(c(161:190)[i]))))$coefficients[1,1]
  results.nonprereg.ascribed[i,13] <- r.squaredGLMM(get(paste0("model.",as.character(c(161:190)[i]))))[1]
  results.nonprereg.ascribed[i,14] <- r.squaredGLMM(get(paste0("model.",as.character(c(161:190)[i]))))[2]
}  
results.nonprereg.ascribed$specification <- 1:30
results.nonprereg.ascribed$specification <- factor(results.nonprereg.ascribed$specification, levels=unique(as.character(results.nonprereg.ascribed$specification)))


table.core.ascribed <- results.core.ascribed[c("specification","estimate","lowerCI","upperCI","p","marginal.R2","conditional.R2")]
table.core.ascribed[,2:7] <- round(table.core.ascribed[,2:7], 3) 
table.expanded.ascribed <- results.expanded.ascribed[c("specification","estimate","lowerCI","upperCI","p","marginal.R2","conditional.R2")]
table.expanded.ascribed[,2:7] <- round(table.expanded.ascribed[,2:7], 3) 
table.nonprereg.ascribed <- results.nonprereg.ascribed[c("specification","estimate","lowerCI","upperCI","p","marginal.R2","conditional.R2")]
table.nonprereg.ascribed[,2:7] <- round(table.nonprereg.ascribed[,2:7], 3) 

write.csv(table.core.ascribed, "table_core_ascribed_2.csv")
write.csv(table.expanded.ascribed, "table_expanded_ascribed_2.csv")
write.csv(table.nonprereg.ascribed, "table_nonprereg_ascribed_2.csv")

mean(table.core.ascribed$estimate)
mean(table.expanded.ascribed$estimate)
mean(table.nonprereg.ascribed$estimate)


# =================
# === PLOTTING === 
# ================

specif.panel <- data.frame(matrix(NA, 30*8, 3))
colnames(specif.panel) <- c("specification","features","present")
specif.panel[,1] <- rep(1:30,8)
specif.panel[,2] <- c(rep("control: population",30),rep("control: education",30),rep("control: GDP per capita",30),
                          rep("domain-level traits",30),rep("subchapter-level traits",30),
                          rep("Uz DS index",30),rep("Uz DG index",30),rep("Muthukrishna index",30))
specif.panel[,3] <- c(
  rep(c(NA,NA,NA,"■","■"), 6), # population
  rep(c(NA,NA,"■",NA,"■"), 6), # education
  rep(c(NA,"■",NA,NA,"■"), 6), # GDPpc 
  rep(c(rep(NA,5), rep("■",5)),3), #domain-level
  rep(c(rep("■",5), rep(NA,5)),3),  #subchapter-level
  c(rep(NA,20), rep("■",10)), # Uz DS index
  c(rep(NA,10), rep("■",10), rep(NA,10)), # Uz DG index
    c(rep("■",10), rep(NA,20)) # Muth index
)
specif.panel$specification <- rep(1:30,8)
specif.panel$specification <- factor(specif.panel$specification, levels=unique(as.character(specif.panel$specification)))
specif.panel$features <- factor(specif.panel$features, levels=unique(as.character(specif.panel$features)))


top.expanded = ggplot(results.expanded , aes(specification, estimate, color = slopesign)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI),  shape = "", alpha = .5) +
  geom_point(size = 1) +
  scale_color_manual(values = c("dodgerblue2", "darkorange1")) +
  labs(x = "", y = "") +
  ggtitle("Inclusive set of cultural traits") +
  scale_y_continuous(limits = c(-0.12, 0.16), breaks = c(-0.1,-0.05,0,0.05,0.1,0.15) ) +
  geom_hline(yintercept=0, linetype=5, color='blue', alpha=0.5) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

top.core = ggplot(results.core , aes(specification, estimate, color = slopesign)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI),  shape = "", alpha = .5) +
  geom_point(size = 1) +
  #scale_x_continuous(limits = c(1, 60), breaks = 1:60) + 
  scale_color_manual(values = c("dodgerblue2", "darkorange1")) +
  labs(x = "", y = "Standardized beta coefficient of cultural variance") +
  ggtitle("Restricted set of cultural traits") +
  geom_hline(yintercept=0, linetype=2, color='blue', alpha=0.5) +
  scale_y_continuous(limits = c(-0.12, 0.16), breaks = c(-0.1,-0.05,0,0.05,0.1,0.15) ) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

top.nonprereg = ggplot(results.nonprereg , aes(specification, estimate, color = slopesign)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI),  shape = "", alpha = .5) +
  geom_point(size = 1) +
  #scale_x_continuous(limits = c(1, 60), breaks = 1:60) + 
  scale_color_manual(values = c("dodgerblue2", "darkorange1")) +
  #labs(x = "", y = "Standardized beta") +
  labs(x = "", y = "") +
  ggtitle("Acultural set of traits") +
  scale_y_continuous(limits = c(-0.12, 0.16), breaks = c(-0.1,-0.05,0,0.05,0.1,0.15) ) +
  geom_hline(yintercept=0, linetype=5, color='blue', alpha=0.5) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# plot bottom panel
bottom.panel = ggplot(specif.panel, aes(specification, features)) +
  geom_text(aes(label = present), color="firebrick", size=5) +
  #scale_x_discrete(limits = c(1, 60), breaks = 1:60) +
  labs(x = "\nSpecification number", y = "Included variables") + 
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# join panels
specificationcurve <-  cowplot::plot_grid(top.expanded, top.core, top.nonprereg, bottom.panel, ncol = 1, align = "v")

ggsave("specificationcurve_04.png", specificationcurve, width=180, height=200, units="mm")

# specificationcurve_expanded <-  cowplot::plot_grid(top.expanded, bottom.panel, ncol = 1, align = "v")
# 
# ggsave("specificationcurve_expanded.png", specificationcurve_expanded, width=180, height=100, units="mm")
# 
# specificationcurve_nonprereg <-  cowplot::plot_grid(top.nonprereg, bottom.panel, ncol = 1, align = "v")
# 
# ggsave("specificationcurve_nonprereg.png", specificationcurve_nonprereg, width=180, height=100, units="mm")


top.expanded.ascribed = ggplot(results.expanded.ascribed , aes(specification, estimate, color = slopesign)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI),  shape = "", alpha = .5) +
  geom_point(size = 1) +
  scale_color_manual(values = c("dodgerblue2", "darkorange1")) +
  labs(x = "", y = "") +
  ggtitle("Inclusive set of cultural traits (cultural variance imputed)") +
  geom_hline(yintercept=0, linetype=2, color='blue', alpha=0.5) +
  scale_y_continuous(limits = c(-0.12, 0.16), breaks = c(-0.1,-0.05,0,0.05,0.1,0.15) ) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

top.core.ascribed = ggplot(results.core.ascribed , aes(specification, estimate, color = slopesign)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI),  shape = "", alpha = .5) +
  geom_point(size = 1) +
  scale_color_manual(values = c("dodgerblue2", "darkorange1")) +
  labs(x = "", y = "Standardized beta coefficient of cultural variance") +
  ggtitle("Restricted set of cultural traits (cultural variance imputed)") +
  scale_y_continuous(limits = c(-0.12, 0.16), breaks = c(-0.1,-0.05,0,0.05,0.1,0.15) ) +
  geom_hline(yintercept=0, linetype=5, color='blue', alpha=0.5) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

top.nonprereg.ascribed = ggplot(results.nonprereg.ascribed , aes(specification, estimate, color = slopesign)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI),  shape = "", alpha = .5) +
  geom_point(size = 1) +
  scale_color_manual(values = c("dodgerblue2", "darkorange1")) +
  labs(x = "", y = "") +
  ggtitle("Acultural set of traits (cultural variance imputed)") +
  scale_y_continuous(limits = c(-0.12, 0.16), breaks = c(-0.1,-0.05,0,0.05,0.1,0.15) ) +
  geom_hline(yintercept=0, linetype=5, color='blue', alpha=0.5) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# plot bottom panel
bottom.panel.ascribed = ggplot(specif.panel, aes(specification, features)) +
  geom_text(aes(label = present), color="firebrick", size=5) +
  #scale_x_discrete(limits = c(1, 60), breaks = 1:60) +
  labs(x = "\nSpecification number", y = "Included variables") + 
  theme_minimal(base_size = 11) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# join panels
specificationcurve.ascribed <-  cowplot::plot_grid(top.expanded.ascribed, top.core.ascribed, top.nonprereg.ascribed, bottom.panel.ascribed, ncol = 1, align = "v")

ggsave("specificationcurve_ascribed_05.png", specificationcurve.ascribed, width=180, height=200, units="mm")




# ====== ===== ==== ==== ====
# ==== ==== ==== ==== =====
#-----------



models <- paste0("model.", 1:280) # enter model numbers here
var.types <- c("Muth index WVS", "Muth index WVS (ascr.)", "Muth index WEVS", "Muth index WEVS (ascr.)", "Uz DG index", "Uz DS index") 
rsquared.list <- matrix(NA,280,3)

for (i in c(seq(11,70,by=10),seq(81,140,by=10),seq(151,210,by=10)) ){
  if (i<71){datascope <- "Core"}
  if (i>=71 & i<141){datascope <- "Expanded"}
  if (i>=141 & i<211){datascope <- "Non-preregistered"}
  if (i>=211){datascope <- "Psychiatric"}
  vartype.id <- (((i+9)/10)-1) %% (length(var.types)+1) 
  if (i %in% c(1,71,141,211)){ tablerows <- c("GDP per cap.", "education", "population")
    }else{ tablerows <- c(var.types[vartype.id], "GDP per cap.", "education", "population") }
  for (j in seq(0,9)){rsquared.list[i+j,1:2] <-  round(r.squaredGLMM(get(models[i+j])), 4)} # get R-squared
  if (i %in% 11:61){ rsquared.list[i:(i+9),3] <- round(rsquared.list[i:(i+9),1] - rsquared.list[1:10,1],4) } # get r-squared gains compared to models without cultural variance
  if (i %in% 81:131){ rsquared.list[i:(i+9),3] <- round(rsquared.list[i:(i+9),1] - rsquared.list[71:80,1],4)}
  if (i %in% 151:201){ rsquared.list[i:(i+9),3] <- round(rsquared.list[i:(i+9),1] - rsquared.list[141:150,1],4)}
  if (i %in% 221:271){ rsquared.list[i:(i+9),3] <- round(rsquared.list[i:(i+9),1] - rsquared.list[211:220,1],4)}
  
  # for some reason, the models that omit the cultural var predictors give an error at stargazer(), so skipping those at the for loop
  # need to detatch lmerTest
   if (!(i %in% c(1,71,141,211))){
  stargazer(get(models[i]),get(models[i+1]),get(models[i+2]),get(models[i+3]),get(models[i+4]),
            get(models[i+5]),get(models[i+6]),get(models[i+7]),get(models[i+8]),get(models[i+9]),
            column.labels = c(paste0(datascope," traits at Subchapter-level"), paste0(datascope," traits at Domain-level")), 
            column.separate = c(5, 5), title= paste0("Models ",i," to ", i+9, ": ", datascope, " trait list, ", var.types[vartype.id]), 
            covariate.labels = tablerows, align=TRUE, ci = TRUE,
            add.lines = list(c("R2 marginal", rsquared.list[i:(i+9),1] ),
                             c("R2 conditional", rsquared.list[i:(i+9),2] ),
                             c("R2 marginal gain", rsquared.list[i:(i+9),3] )),
            out= paste0("modelresults_",i,"_to_",i+9, ".html" ))
  }
}


# for manuscript figure
model.85m <- lmer(heritability ~ var_muth_wvs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.95m <- lmer(heritability ~ var_muth_wvs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.105m <- lmer(heritability ~ var_muth_wevs+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.115m <- lmer(heritability ~ var_muth_wevs.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.125m <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.135m <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
var.types <- c("Muth index WVS", "Muth index WVS (ascr.)", "Muth index WEVS", "Muth index WEVS (ascr.)", "Uz DG index", "Uz DS index") 
rsquared.list.m <- matrix(NA,6,2)
rsquared.list.m[1,1:2] <-  round(r.squaredGLMM(get("model.85m")), 4)
rsquared.list.m[2,1:2] <-  round(r.squaredGLMM(get("model.95m")), 4)
rsquared.list.m[3,1:2] <-  round(r.squaredGLMM(get("model.105m")), 4)
rsquared.list.m[4,1:2] <-  round(r.squaredGLMM(get("model.115m")), 4)
rsquared.list.m[5,1:2] <-  round(r.squaredGLMM(get("model.125m")), 4)
rsquared.list.m[6,1:2] <-  round(r.squaredGLMM(get("model.135m")), 4)
tablerows <- c("cultural variance", "GDP per cap.", "education", "population")

stargazer(get("model.85m"),get("model.95m"),get("model.105m"),get("model.115m"),get("model.125m"),get("model.135m"),          
          title= "Models results under Expanded trait list and subchapter-level traits", 
          covariate.labels = c("Muth index raw (WVS)", "Muth index ascribed (WVS)", "Muth index raw (W+EVS)", "Muth index ascribed (W+EVS)", 
                               "Uz DG index", "Uz DS index", "GDP per cap.", "education", "population"), 
          align=TRUE, ci = TRUE,
          add.lines = list(c("R2 marginal", rsquared.list.m[1:6,1] ),
                           c("R2 conditional", rsquared.list.m[1:6,2] )) ,
          out= "modelresults_manuscript.html" )

# stargazer(model.85m, title= "Models results under", align=TRUE, ci = TRUE,
#           add.lines = list(c("R2 marginal", rsquared.list[1:6,1] ),
#                            c("R2 conditional", rsquared.list[1:6,2] )),
#                            out= "modelresults_manuscript.html" )

# 85,95,105,115,125,135

# = )) = )) = )) = )) = )) = )) = )) = )) = )) = )) = )) = )) = )) = )) = ))


coef(model.41)









ggsave("hist.cultvar.1.png", hist.cultvar.1, width=5, height=2)


summary(model.11)$ngrps



$ ngrps       : Named num [1:3] 17 10 4
..- attr(*, "names")= chr [1:3] "trait_subchapter" "country" "h2_type"


# stack two tables?
#
#-------
results.list$m.R2gain <- results.list$m.R2 - rep( results.list$m.R2[c(1:10,51:60)], 4)

summary(model.99)

add.lines = list(c("Model number", seq(i,i+9))),


for (j in vec){class(get(model.list[99])) <- "lmerMod"}



class(get(model.list[99])) <- "lmerMod"

summary(model.100)$coefs


data(cake)

insertrow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

summary(M1 <- lmer(angle ~ temp + (1 | replicate) + (1|recipe:replicate), cake, REML= FALSE))
summary(M2 <- lmer(angle ~ factor(temperature) + (1 | replicate) + (1|recipe:replicate), cake, REML= FALSE))

Tables <- stargazer(M1, M2, style="ajps", title="An Illustrative Model Using Cake Data", dep.var.labels.include = FALSE, 
                    covariate.labels=c( "Temperature (Continuous)", "Temperature (Factor $<$ 185)", "Temperature (Factor $<$ 195)", "Temperature (Factor $<$ 205)", "Temperature (Factor $<$ 215)", "Temperature (Factor $<$ 225)")
)

Tables <- as.data.frame(Tables)
Tables$Tables <- as.character(Tables$Tables)
Tables

r <- 25

randomeffect <- "{\\bf Random Effect} & & \\\\"
hline <- "\\hline"
newline <- "\\\\"

Tables <- insertrow(Tables, hline, r)
Tables <- insertrow(Tables,randomeffect,r+1)
Tables <- insertrow(Tables,hline,r+2)

num.recipe.replicate <- sapply(ranef(M1),nrow)[1]
num.replicate <- sapply(ranef(M1),nrow)[2]

stddev.M1.recipe.replicate <- attributes(VarCorr(M1)$"recipe:replicate")$stddev
stddev.M1.replicate <- attributes(VarCorr(M1)$replicate)$stddev
stddev.M2.recipe.replicate <- attributes(VarCorr(M2)$"recipe:replicate")$stddev
stddev.M2.replicate <- attributes(VarCorr(M2)$replicate)$stddev

number.of.recipe.replicate <- paste("\\# of Recipe:Replicate & ", num.recipe.replicate, "&", num.recipe.replicate, "\\\\")
stddev.recipe.replicate <- paste("Recipe:Replicate Standard Deviation & ", round(stddev.M1.recipe.replicate, 3), "&", round(stddev.M2.recipe.replicate, 3), "\\\\")
number.of.replicate <-paste("\\# of Replicate & ", num.replicate, "&", num.replicate, "\\\\")
stddev.replicate <- paste("Replicate Standard Deviation & ", round(stddev.M1.replicate, 3), "&", round(stddev.M2.replicate, 3), "\\\\")

Tables <- insertrow(Tables,number.of.recipe.replicate,r+3)
Tables <- insertrow(Tables,stddev.recipe.replicate,r+4)
Tables <- insertrow(Tables,newline,r+5)
Tables <- insertrow(Tables,number.of.replicate,r+6)
Tables <- insertrow(Tables,stddev.replicate,r+7)

write.table(Tables,file="tables.tex",sep="",row.names= FALSE,na="", quote = FALSE, col.names = FALSE)



#--------



stargazer(model.51, out="heri_test.tex")

class(model.51)<- "lmerMod"


model.list <- paste0("model.", 1:100)
results.list <- data.frame(matrix(NA, 100, 8))
colnames(results.list) <- c("beta","SE","df","t","p","m.R2","c.R2","m.R2gain")
for (i in 1:100){
  print(i)
  results.list[i,1] <- i  
  if (i %in% c(11:50,61:100)  ){
    results.list[i,1:5] <- round(summary(get(model.list[i]))$coefficients[2,], 4) 
  } else {
    results.list[i,1:5] <- rep(NA,5)
  }
  results.list[i,6:7] <- round(r.squaredGLMM(get(model.list[i])), 4)
}
results.list$m.R2gain <- results.list$m.R2 - rep( results.list$m.R2[c(1:10,51:60)], 4)

full.list <- data.frame(matrix(NA,dim(results.list)[1],7) , results.list)
colnames(full.list)[1:7] <- c("model","trait.breadth","cult.var","trait.level","gdp.pc","education","population")
full.list$model <- 1:100
full.list$trait.breadth <- c(rep("core",50),rep("expanded",50))
full.list$cult.var <- rep(c(rep("none",10),rep("muth",10),rep("muth.ascr",10),rep("uz.DG",10),rep("uz.DS",10)),2)
full.list$trait.level <- rep(c(rep("subchap",5),rep("domain",5)),10)
full.list$gdp.pc <- rep(c(0,1,0,0,1),20)
full.list$education <- rep(c(0,0,1,0,1),20)
full.list$population <- rep(c(0,0,0,1,1),20)

formattable(full.list)

formattable(full.list, list(trait.breadth = color_tile("white", "pink")))

formattable(full.list, align =c(rep("c",15)), list(
  `beta`= color_tile("white", "green"),
  `SE`= color_tile("white", "green"),
  `df`= color_tile("white", "green"),
  `t`= color_tile("white", "green"),
  `p`= color_tile("white", "green"),
  `m.R2`= color_tile("white", "orange"),
  `c.R2`= color_tile("white", "orange"),
  `m.R2gain`= color_tile("white", "orange"),
  `gdp.pc`= color_tile("white", "darkmagenta"),
  `education`= color_tile("white", "darkmagenta"),
  `population`= color_tile("white", "darkmagenta")
))



summary(model.17)


#-------------







stargazer(model.33, out="heri_test.htm")

summary(model.34)$coef[1:12]




anova(model.39, model.99)

summary(model.11)


library(xtable)
# Get the table first.
summary(M1 <- lme4::lmer(angle ~ temp + (1 | replicate) + (1|recipe:replicate), cake, REML= FALSE))
summary(M2 <- lme4::lmer(angle ~ factor(temperature) + (1 | replicate) + (1|recipe:replicate), cake, REML= FALSE))

stargazer(M1, M2, style="ajps", title="An Illustrative Model Using Cake Data", dep.var.labels.include = FALSE, 
          covariate.labels=c( "Temperature (Continuous)", "Temperature (Factor $<$ 185)", "Temperature (Factor $<$ 195)", "Temperature (Factor $<$ 205)", "Temperature (Factor $<$ 215)", "Temperature (Factor $<$ 225)")
)



# now for lmerTest
summary(M1a <- lmer(angle ~ temp + (1 | replicate) + (1|recipe:replicate), cake, REML= FALSE))
summary(M2a <- lmer(angle ~ factor(temperature) + (1 | replicate) + (1|recipe:replicate), cake, REML= FALSE))
anovadf <- data.frame(anova(M1a,M2a))
xtable(anovadf)
print(anovadf, type = "html")






#-----------
model.list <- paste0("model.", 1:90)
results.list <- data.frame(matrix(NA, 90, 9))
colnames(results.list) <- c( "model", "cultvar.estimate","cultvar.SE","cultvar.df","cultvar.t-value","cultvar.p-value","marginal.R2","conditional.R2","marginal.R2gain")
for (i in 1:90){
  results.list[i,1] <- model.list[i]  
  if (i<81){
  results.list[i,2:6] <- round(summary(get(model.list[i]))$coefficients[2,], 4) 
  } else {
    results.list[i,2:6] <- rep(NA,5)
  }
    results.list[i,7:8] <- round(r.squaredGLMM(get(model.list[i])), 4)
}
results.list$marginal.R2gain <- results.list$marginal.R2 - rep( results.list$marginal.R2[81:90], 9)


results.list <- t(sapply(model.list, function(x)r.squaredGLMM(get(x))))
results.list <- data.frame(round(results.list, 3))
colnames(results.list) <- c("marginal.R2","conditional.R2")
results.list$marginal.R2gain <- rep(NA, dim(results.list)[1])
results.list$marginal.R2gain <- results.list$marginal.R2 - rep( results.list$marginal.R2[81:90], 9)

results.list



# var muth - expanded
# model 23 * (control variables: education) ; beta -7.0 !?
# model 24 . (control variables: population) ; -4.8

# var muth ascribed - expanded
# model 32 . ; -3.9
# model 33 * ; -7.3
# model 34 . ; -4.9
# model 35 . ; -6.3

# var muth not significant in any of the "domain"-level phenotype models
# var uz DG (models 40-49) doesn't work with dat.core

# var uz DG - expanded
# model 51 *
# model 52 *
# model 53 *
# model 54 *
# model 55 .
# model 56 *
# model 57 .
# model 58 .
# model 59 .

# var uz DS non-significant in every model













#results.list$conditional.R2gain <- rep(NA, dim(results.list)[1])
##results.list$marginal.R2gain <- results.list$marginal.R2

# results.list$conditional.R2[1] - results.list$conditional.R2[61]
# results.list$conditional.R2[1] / results.list$conditional.R2[61]
# results.list$conditional.R2[5] - results.list$conditional.R2[65]
# results.list$conditional.R2[5] / results.list$conditional.R2[65]
# 
# 
# 
# 
# baseline <- t(sapply(paste0("model.", 61:65), function(x)r.squaredGLMM(get(x))))
# rbind(rep(baseline,12), nrow=nrow(baseline), bycol=T)
# 
# 
# matrix( rep( t( baseline ) , 10 ) , ncol = ncol(m) , byrow = TRUE )
# 
# rsquared.list$test <- 
# # subchapter higher r2 than domain. Core list higher than expanded list. uzDG generally higher than uzDS, but uzDS often higher|Domain
# 
# 
# 
# #******
# colnames(dat.long)
# table(dat.cultural_withpsychiatric$country_name)
# table(dat.cultural$trait_chapter)
# 
# 
# table(dat.cultural_psychiatric$country_name)
# 
# hist(unique(dat.cultural$var_uzDG), breaks=20)
# 
# 
# table(dat.long$trait_domain)
# 
# 
# 
# 
# 
# dat.long.1 <- dat.long[culturaltraits_nonpsychiatric,]
# dat.long.2 <-dat.long[culturaltraits_withpsychiatric,]
# 
# tabulate(dat.long.1$e_regionpol)
# tabulate(dat.long.2$e_regionpol)
# 
# levels(as.factor(dat.long.2$country_name))
# tabulate(as.factor(dat.long.2$country_name))
# 
# 
# 
#   length(culturaltraits_withpsychiatric)
#   
#   x <- dat.long[culturaltraits_nonpsychiatric,]
#   length(unique(x$country_name))
#   
#     
#     
# )
# 
# 
# 
# trait_hierarchy2 <- unique( dat.long[,14:12] )
# trait_hierarchy2 <- trait_hierarchy2[order(trait_hierarchy2$trait_subchapter) , ]
# trait_hierarchy2 <- trait_hierarchy2[order(trait_hierarchy2$trait_chapter) , ]
# trait_hierarchy2 <- trait_hierarchy2[order(trait_hierarchy2$trait_domain) , ]
# write.csv(trait_hierarchy2, "actual_traits_foranalysis.csv")
# 
# trait_hierarchy3 <- unique( dat.long[,12:11] )
# trait_hierarchy3 <- trait_hierarchy3[order(trait_hierarchy3$trait_subchapter) , ]
# write.csv(trait_hierarchy3, "subchapter_and_trait.csv")
# 
# 
# 
# 
# lmm <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ trait+ (1+ trait+ var_uzDG| e_regionpol))
# 
# 
# 
# lmm <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ trait+ )
# 
# 
# 



