


assign("last.warning", NULL, envir = baseenv())
# reset errors (sigh)



#----
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
library(formattable)

original_wd <- getwd()
# Import WVS data used for Muthukrishna et al. Cultural Distance analysis
### wvs.alleles <- read.csv("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/wvs_allelesv02_merged.csv")
# Import raw (pre-'allelized') WVS data used for Muthukrishna et al. Cultural Distance analysis
wvs <- read.dta13("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/WVS_Longitudinal_1981_2014_stata_v2015_04_18.dta")
# Import WVS variable characteristics used for Muthukrishna et al. Cultural Distance analysis
dimensions <- read.csv("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/Additional_Supplementary_Materials/allele-dimensions-data.csv")
# Import Uz (2015) tightness-looseness data 
uz <- read.csv("~/Desktop/Muthukrishna_Lab/sociality-IQ/DATA/uz_tightness_withISO.csv")
# Import Varieties of Democracy 'V-Dem Full+Others' data set
vdem <- readRDS("~/Desktop/Muthukrishna_Lab/sociality-IQ/DATA/Country_Year_V-Dem_Full+others_R_v9/V-Dem-CY-Full+others-v9.rds")
# Import table 33 of Polderman et al. (2013), in supplementary material (.xlsx file) 
trait_hierarchy <- read_excel("~/Desktop/Muthukrishna_Lab/sociality-IQ/DATA/Polderman2013 - heritability metaanalysis/polderman2013-supp2.xlsx", sheet = "Supplementary Table 33", range = "A2:C328")
# Set to directory containing all MaTCH ICF/ICD10 subchapter (by country) tables 
setwd("~/Desktop/Muthukrishna_Lab/sociality-IQ/DATA/MaTCH all traits no constraints")

# Extract all country codes included in all country-based tables available for download from MaTCH  
filenames <- list.files()
country <- as.character(matrix(NA, 1, length(filenames)))
for (i in 1:length(filenames)){
  file <- read.csv(filenames[i], sep=";")
  if (dim(file)[1]>0){
    country[i] <- as.character(read.csv(filenames[i], sep=";")$Code) }}
country <- unique(country)
country <- as.factor(country[!is.na(country)])

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
setwd(original_wd) # read files already so set working directory back to what it was initially 

dat.pre <-  rbind(dat.pre.h2_all, dat.pre.h2_ss, dat.pre.h2_m, dat.pre.h2_f)
dat.pre <- cbind(rep(c("h2_all","h2_ss","h2_m","h2_f"), each= dim(dat.pre.h2_all)[1])  , dat.pre)
colnames(dat.pre)[1] <- "h2_type"
 
# Remove columns containing no data in any of the h2 types - left with 121 traits (first column is just countries)
columns.nodata <- matrix(0,0,0)
for (i in 1:dim(dat.pre)[2]){
  if (sum(!is.na(dat.pre[,i])) == 0){ columns.nodata <- c(columns.nodata, i)}}
dat.pre <- dat.pre[,-columns.nodata]

# All variables that appear in "recode_wvs.do" 
include_vars <-  c("A001","A002","A003","A004","A005","A006","A007","A025","A026","A029","A030","A032","A034","A035","A038","A039","A040","A041","A042","A165","A168","A057","A058","A059","A060","A061","A062","B001","B002","B003","B008","B009","A169",
"A064","A065","A066","A067","A068","A069","A070","A071","A072","A073","A074","A075","A076","A077","A079","A081","A082","A083","A084","A085","A086","A087","A088","A089","A090","A091","A092","A093","A094","C001","C002","C006","A170","A173",
"C008","C009","C010","C011","C012","C013","C014","C015","C016","C017","C018","C019","C020","C021","C036","C037","C038","C039","C040","C041","C059","C060","C061","X007","X011","D017","D018","D019","D022","D023","D054","D055","D056","D057","D058","D059","D060",
"E001","E002","E003","E004","E005","E006","E012","E014","E015","E016","E018","E019","E022","E023","E025","E026","E027","E028","E029","E033","E034","E035","E036","E037","E039","E046","E143","E112","E114","E115","E116","E117","E110",
"E120","E121","E122","E123","E124","E125","E128","E129","E135","E136","E137","E138","E139","F001","F022","F025","F028","F034","F035","F036","F037","F038","F050","F051","F052","F053","F054","F063","F064","F065","F066","F102","F103",
"F104","F105","F114","F115","F116","F117","F118","F119","F120","F121","F122","F123","G001","G002","G006","E150","G015","G016","E179","E180","E182")

wvs <- wvs[,(names(wvs) %in% include_vars | names(wvs) %in% c("S002", "S003", "X048WVS"))] #year, country, region
wvs <- wvs[,c(which(names(wvs) %in% c("S002", "S003", "X048WVS")), which(!(names(wvs) %in% c("S002", "S003", "X048WVS"))))] # put year country region in front
colnames(wvs)[which(names(wvs) == "S002")] <- "wave"
colnames(wvs)[which(names(wvs) == "S003")] <- "pop"
colnames(wvs)[which(names(wvs) == "X048WVS")] <- "region"
col_idx <- c(which(names(wvs)=="pop"), which(names(wvs)=="region"), which(names(wvs)=="wave")) # To move these variables to the front
N.extravars <- length(col_idx)
##wvs <- wvs[, c(col_idx, (1:ncol(wvs))[-col_idx]) ]
wvs <- wvs[,c( col_idx, setdiff(1:length(names(wvs)),col_idx) ) ]

wvs <- wvs[str_detect(wvs$wave,"^2005") | str_detect(wvs$wave,"^2010"),] # use just these WVS waves
wvs[,(N.extravars+1):dim(wvs)[2]] <- wvs[,match(include_vars, names(wvs))] # re-order everything except for year, country, region so that it matches the order of variables in "include_vars", which is also the order in "recode_wvs.do" 
wvs.variablenames <- paste0("v",c( 1:7,10:64,75:143,160:219)) # assigning variable names according to Stata script "recode_wvs.do" 
names(wvs)[(N.extravars+1):dim(wvs)[2]] <- wvs.variablenames
wvs.nominalvars <- c("v11","v33","v34","v35","v81","v82","v101","v103","v117","v118","v119","v120","v121","v122",
  "v174","v175","v176","v177","v178","v180","v181","v215","v216","v217","v218","v219")
# Note: v10 can be treated as binary because nobody responded with "neither", which we would have considered nominal
wvs <- wvs[,-which(names(wvs) %in% wvs.nominalvars)] # remove nominal variables

# ===== Recoding of factors (necessary for normalization of variables) =====
# Get rid of factor level 10 "neither" since it is not used (and since it prevents us from using the variable as a binary)
levels(wvs$v10)[7:8] <- levels(wvs$v10)[7]
# v75: rearrange order of responses from "agree, disagree, neither" so that "neither" is in the middle.
levels(wvs$v75)[7:8] <- c(levels(wvs$v75)[8], levels(wvs$v75)[7])
v75_resp7 <-  which(wvs$v75==levels(wvs$v75)[7])
v75_resp8 <-  which(wvs$v75==levels(wvs$v75)[8])
wvs$v75[v75_resp7] <- levels(wvs$v75)[8]  
wvs$v75[v75_resp8] <- levels(wvs$v75)[7]
# v76: rearrange order of responses from "agree, disagree, neither" so that "neither" is in the middle.
levels(wvs$v76)[7:8] <- c(levels(wvs$v76)[8], levels(wvs$v76)[7])
v76_resp7 <-  which(wvs$v76==levels(wvs$v76)[7])
v76_resp8 <-  which(wvs$v76==levels(wvs$v76)[8])
wvs$v76[v76_resp7] <- levels(wvs$v76)[8]  
wvs$v76[v76_resp8] <- levels(wvs$v76)[7]
# v102: rearrange order of responses from "follow instrucitons, must be convinced, depends" so that "depends" in in the middle
levels(wvs$v102)[7:8] <- c(levels(wvs$v102)[8], levels(wvs$v102)[7])
v102_resp7 <-  which(wvs$v102==levels(wvs$v102)[7])
v102_resp8 <-  which(wvs$v102==levels(wvs$v102)[8])
wvs$v102[v102_resp7] <- levels(wvs$v102)[8]  
wvs$v102[v102_resp8] <- levels(wvs$v102)[7]
# v109: rearrange order of responses from "disapprove, approve, depends" so that "depends" in in the middle
levels(wvs$v109)[7:8] <- c(levels(wvs$v109)[8], levels(wvs$v109)[7])
v109_resp7 <-  which(wvs$v109==levels(wvs$v109)[7])
v109_resp8 <-  which(wvs$v109==levels(wvs$v109)[8])
wvs$v109[v109_resp7] <- levels(wvs$v109)[8]  
wvs$v109[v109_resp8] <- levels(wvs$v109)[7]
# v129: rearrange order of responses from "will help, will harm, some of each" so that "some of each" in the middle
levels(wvs$v129)[7:8] <- c(levels(wvs$v129)[8], levels(wvs$v129)[7])
v129_resp7 <-  which(wvs$v129==levels(wvs$v129)[7])
v129_resp8 <-  which(wvs$v129==levels(wvs$v129)[8])
wvs$v129[v129_resp7] <- levels(wvs$v129)[8]  
wvs$v129[v129_resp8] <- levels(wvs$v129)[7]
# v182, need to collapse "Only on special holy days/Christmas/Easter days" and "Other specific holy days" together
levels(wvs$v182) <- levels(wvs$v182)[c(1:8,10,10,11:13)] #merge the two responses

# ===== Min-max normalization =====
for (i in (N.extravars+1):dim(wvs)[2] ){  
  wvs[,i] <- as.integer(wvs[,i]) # first convert variable from factor to integer
  wvs.norm <- wvs[,i]
  wvs.norm[which(wvs.norm<6)] <- NA # all negative-valued responses get turned into NA
  wvs[,i] <- (wvs.norm - min(wvs.norm, na.rm=T))/ (max(wvs.norm, na.rm=T) - min(wvs.norm, na.rm=T)) # perform min-max normalisation
}
wvs.colsums <- colSums(wvs[,(N.extravars+1):dim(wvs)[2]], na.rm=T) # remove variables that were all NA, i.e., negative values in the original WVS 
wvs <- wvs[, -as.vector(which(wvs.colsums==0)+3)] # adding 3 because index is from 4: above, to exclude pop, location, region 
pops <- unique(wvs$pop)

# Get modified ('non-allelised')  Muthukrishna cultural variance
cultvar.matrix <- matrix(NA, length(pops), dim(wvs)[2]-N.extravars) #rows are populations, columns are WVS questions
for (col in 1:dim(cultvar.matrix)[2]){
  print(col)
  cultvar.matrix[,col] <-  sapply(1:length(pops), function(y){ var(wvs[which(wvs$pop==pops[y]), col+N.extravars ], na.rm=T)})
  }
cultvar.means <- rowMeans(cultvar.matrix, na.rm=T)
var.societies <- data.frame(pops, "CountryISO"= rep(NA, length(cultvar.means)) , cultvar.means) 
var.societies <- var.societies[order(var.societies$pops),]
# attaching the ISO codes only to countries for which we have heritability data 
var.societies$CountryISO[c( 5,11,13,20,21,29,33,34,37,47,50,57,64,65,73,74)] <-  
  c("AU","CA","CN","FI","FR","IN","IT","JP","KR","NL","NO","RU","ES","SE","GB","US")
var.societies$CountryISO <- as.factor(var.societies$CountryISO)
var.societies.reduced <- var.societies[var.societies$CountryISO %in% country,] # only inlcude countries for which we have heritability data

# Merge Uz 2015 tightness-looseness data with MaTCH data
dat <- merge(uz[uz$CountryISO %in% country,], dat.pre, by.x="CountryISO", by.y="country", all.y = T)
dat <- data.frame(dat[,1:2], matrix(NA,dim(dat)[1],2), dat$CTL_DS, rep(NA,dim(dat)[1]), dat$CTL_DG, rep(NA,dim(dat)[1]), dat$CTL_C, rep(NA,dim(dat)[1]), dat[,-c(1:5)])
colnames(dat)[3:10] <- c("var_muth","var_muth.ascribed","var_uzDS","var_uzDS.ascribed","var_uzDG","var_uzDG.ascribed","var_uzC","var_uzC.ascribed")
dat$var_muth <- var.societies.reduced$cultvar.means[match(dat$CountryISO, var.societies.reduced$CountryISO)]

# organising the data
dat$var_uzDS[which(dat$var_uzDS=="n/a")] <- NA # Convert Uz "n/a"s into proper NAs
dat$var_uzDG[which(dat$var_uzDG=="n/a")] <- NA
dat$var_uzC[which(dat$var_uzC=="n/a")] <- NA
dat$var_uzDS <- as.numeric(as.character(dat$var_uzDS))# CTL values imported as factors so converting into numeric 
dat$var_uzDG <- as.numeric(as.character(dat$var_uzDG))# CTL values imported as factors so converting into numeric 
dat$var_uzC <- as.numeric(as.character(dat$var_uzC))# CTL values imported as factors so converting into numeric 
dat <- dat[order(as.character(dat$CountryISO)),] # reorder by country ISO code

# Which countries are missing from Muthukrishna cultural variance?
print(paste0("Countries that we have heritability data for, but which are missing from Muthukrishna cultural variance: ", as.character(unique(dat$CountryISO[which(is.na(dat$var_muth))])) ))  
print(paste0("Countries that we have heritability data for, but which are missing from Uz cultural variance: ", as.character(unique(dat$CountryISO[which(is.na(dat$var_uzDG) | is.na(dat$var_uzDS))]))))  
dat$var_muth.ascribed <- dat$var_muth
dat$var_uzDS.ascribed <- dat$var_uzDS
dat$var_uzDG.ascribed <- dat$var_uzDG
dat$var_uzC.ascribed <- dat$var_uzC

# Manually entering countries to use for ascription of Muthukrishna cultrual variance 
dat$var_muth.ascribed[which(dat$Country=="Belgium")] <- mean(var.societies$cultvar.means[which(var.societies$pops %in% c("France","Germany","Netherlands"))])
dat$var_muth.ascribed[which(dat$Country=="Denmark")] <- mean(var.societies$cultvar.means[which(var.societies$pops %in% c("Sweden","Germany"))])
# --> Skipping ascription for Czech Republic because heritability from there isn't included in our list of culturally sensitive traits

# Manually entering countries to use for ascription of Uz cultrual variance 
dat$var_uzDS.ascribed[which(dat$CountryISO=="AU")] <- mean(dat$var_uzDS[which(dat$Country=="Great Britain")[1]])
dat$var_uzDG.ascribed[which(dat$CountryISO=="AU")] <- mean(dat$var_uzDG[which(dat$Country=="Great Britain")[1]])
dat$var_uzC.ascribed[which(dat$CountryISO=="AU")] <- mean(dat$var_uzC[which(dat$Country=="Great Britain")[1]])
dat$var_uzDS.ascribed[which(dat$CountryISO=="NO")] <- mean(dat$var_uzDS[c(which(dat$Country=="Denmark")[1], which(dat$Country=="Sweden")[1])])
dat$var_uzDG.ascribed[which(dat$CountryISO=="NO")] <- mean(dat$var_uzDG[c(which(dat$Country=="Denmark")[1], which(dat$Country=="Sweden")[1])])
dat$var_uzC.ascribed[which(dat$CountryISO=="NO")] <- mean(dat$var_uzC[c(which(dat$Country=="Denmark")[1], which(dat$Country=="Sweden")[1])])

# subset of V-Dem dataset, for 2006 (median year of all papers included in MaTCH meta-analysis) and for countries remaining from MaTCH data
vdem.sub <- rbind(vdem[which(vdem$country_name == "Australia" & vdem$year == 2006),], 
                  vdem[which(vdem$country_name == "Belgium" & vdem$year == 2006),], 
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
                  vdem[which(vdem$country_name == "Norway" & vdem$year == 2006),], 
                  vdem[which(vdem$country_name == "Russia" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "Sweden" & vdem$year == 2006),],
                  vdem[which(vdem$country_name == "United States of America" & vdem$year == 2006),]
) 
# further subsetting, taking only variables of interest from full V-Dem dataset
col.select <- c(  
which(colnames(vdem.sub) == "country_name"),
which(colnames(vdem.sub) == "e_migdppc"), #GDP per capita
which(colnames(vdem.sub) == "e_peaveduc"), # years education for citizens older than 15 (missing for Croatia? but Croatia data not in MaTCH anyway)
which(colnames(vdem.sub) == "e_wb_pop"), # population
which(colnames(vdem.sub) == "e_regionpol") # region (politico-geographic). Less sub-divisions than 'geographic region' (e_regiongeo)
)

# Merge MaTCH and Uz data with V-Dem data and clean up data matrix
vdem.selected <- vdem.sub[,col.select]
vdem.selected <- vdem.selected[rep(1:nrow(vdem.selected), each=4),] # multiplying V-Dem rows for each of the four h2 types.
dat <- cbind(dat, vdem.selected)
dat <- dat[,-which(names(dat)== "CountryISO")] # remove columns for "CountryISO" and "Country", as we use "country_name" in the final data frame
dat <- dat[,-which(names(dat)== "Country")]
colnames(dat)[which(colnames(dat)=="country_name")] <- "country" # rename variable
###dat.long <- melt(dat, na.rm=T, id.vars=c("h2_type","country_name", "var_uzDS", "var_uzDG", "var_uzC", "e_migdppc", "e_peaveduc", "e_wb_pop", "e_regionpol")) # Use melt() to reorganise matrix
# ^ 535 obs, 14 countries, 112 unique subchapters, 
dat.long <- melt(dat, na.rm=T, id.vars=c("h2_type","country", "var_muth","var_muth.ascribed","var_uzDS","var_uzDS.ascribed", "var_uzDG","var_uzDG.ascribed", "var_uzC", "var_uzC.ascribed", "e_migdppc", "e_peaveduc", "e_wb_pop", "e_regionpol")) # Use melt() to reorganise matrix
dat.long[,c(dim(dat.long)[2],dim(dat.long)[2]-1)] <- dat.long[,c(dim(dat.long)[2]-1,dim(dat.long)[2])] # simply putting trait trait name after heritability
colnames(dat.long)[c(dim(dat.long)[2]-1,dim(dat.long)[2])] <- c("heritability","trait_string") # renaming "value" and "variable"

# trait labels in full verbal form (not just strings extracted from file names), from Polderman et al 2015 supplementary table 33. Using these to fetch higher-order categories as well the same table.  
supp32traitlabels <-  c("Acquired Absence of Organs, Not Elsewhere Classified","All-Cause Mortality","Allergy, Unspecified","Asthma","Atopic Dermatitis","Attention Functions","Basic Interpersonal Interactions","Bipolar Affective Disorder","Blood Pressure Functions","Blood Vessel Functions","Calculation Functions","Complex Interpersonal Interactions","Conduct Disorders","Coxarthrosis [Arthrosis of Hip]",
"Cystic Fibrosis","Dementia In Alzheimer Disease","Depressive Episode","Diseases of Pulp and Periapical Tissues","Disorders of Social Functioning with Onset Specific to Childhood and Adolescence","Dissocial Personality Disorder","Dorsalgia","Eating Disorders","Education","Emotional Disorders with Onset Specific to Childhood",
"Emotionally Unstable Personality Disorder","Endocrine Gland Functions","Energy and Drive Functions","Exercise Tolerance Functions","Experience of Self and Time Functions","Expression","Food","Function of Brain","Gene Expression","General Metabolic Functions","Gestational [Pregnancy-Induced] Hypertension Without Significant Proteinuria","Global Psychosocial Functions",
"Gout","Habit and Impulse Disorders","Haematological System Functions","Hearing Functions","Heart Functions","Height","Higher-Level Cognitive Functions","Hyperkinetic Disorders","Hyperplasia of Prostate","Immunological System Functions","Individual Attitudes of Strangers","Informal Social Relationships","Intellectual Functions",
"Intimate Relationships","Irritable Bowel Syndrome","Labour and Delivery Complicated By Umbilical Cord Complications","Leiomyoma of Uterus","Looking After One's Health","Memory Functions","Menstruation Functions","Mental and Behavioural Disorders Due to Multiple Drug Use and Use of Other Psychoactive Substances",
"Mental and Behavioural Disorders Due to Use of Alcohol","Mental and Behavioural Disorders Due to Use of Cannabinoids","Mental and Behavioural Disorders Due to Use of Cocaine","Mental and Behavioural Disorders Due to Use of Hallucinogens","Mental and Behavioural Disorders Due to Use of Opioids","Mental and Behavioural Disorders Due to Use of Sedatives Or Hypnotics",
"Mental and Behavioural Disorders Due to Use of Tobacco","Mental Functions, Unspecified","Migraine","Mild Mental Retardation","Mood [Affective] Disorders","Mortality From Heart Disease","Muscle Power Functions","Nonsuppurative Otitis Media","Obsessive-Compulsive Disorder","Osteoporosis In Diseases Classified Elsewhere","Other Anxiety Disorders",
"Other Functions of the Skin","Other Nontoxic Goitre","Other Symptoms and Signs Involving the Urinary System","Pain and Other Conditions Associated with Female Genital Organs and Menstrual Cycle","Parkinson Disease","Pervasive Developmental Disorders","Phobic Anxiety Disorders","Potential Health Hazards Related to Socioeconomic and Psychosocial Circumstances",
"Problems Related to Upbringing","Procreation Functions","Protective Functions of the Skin","Psychological and Behavioural Disorders Associated with Sexual Development and Orientation","Psychomotor Functions","Reaction to Severe Stress, and Adjustment Disorders","Recreation and Leisure","Recurrent Depressive Disorder","Religion and Spirituality","Schizophrenia","Schizotypal Disorder","Seeing Functions","Senile Cataract","Sensation of Pain","Sexual Functions","Sleep Disorders",
"Sleep Functions","Slow Fetal Growth and Fetal Malnutrition","Societal Attitudes","Specific Developmental Disorder of Motor Function","Specific Personality Disorders","Structure of Areas of Skin","Structure of Brain","Structure of Cardiovascular System","Structure of Eyeball","Structure of Head and Neck Region","Structure of Lower Extremity","Structure of Pelvic Region","Structure of Trunk",
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
  "attention_functions", "basic_interpersonal_interactions","calculation_functions","complex_interpersonal_interactions","education","global_psychosocial_functions","higher-level_cognitive_functions","individual_attitudes_of_strangers","informal_social_relationships","intellectual_functions",
  "looking_after_ones_health","memory_functions","mild_mental_retardation","potential_health_hazards_related_to_socioeconomic_and_psychosocial_circumstances","problems_related_to_upbringing","psychomotor_functions","religion_and_spirituality","societal_attitudes")
  )
# Standardize the variables
dat.long$var_muth <- scale(dat.long$var_muth)
dat.long$var_muth.ascribed <- scale(dat.long$var_muth.ascribed)
dat.long$var_uzDS <- scale(dat.long$var_uzDS)
dat.long$var_uzDG <- scale(dat.long$var_uzDG)
dat.long$var_uzC <- scale(dat.long$var_uzC)
dat.long$e_migdppc <- scale(dat.long$e_migdppc) 
dat.long$e_peaveduc <- scale(dat.long$e_peaveduc) 
dat.long$e_wb_pop <- scale(dat.long$e_wb_pop)
# Adding to the above list all traits within Psychiatric domain except for 'sleep functions'
culturaltraits_withpsychiatric <- union(culturaltraits_nonpsychiatric, setdiff(which(dat.long$trait_domain=="Psychiatric"), which(dat.long$trait_subchapter=="Sleep Functions")))
dat.core <- dat.long[culturaltraits_nonpsychiatric,]
dat.expanded <- dat.long[culturaltraits_withpsychiatric,]

#------
write.csv(dat.core, "heritability-dat-core.csv")
write.csv(dat.expanded, "heritability-dat-expanded.csv")

cor.test(dat.core$var_muth, dat.core$heritability)
cor.test(dat.expanded$var_muth, dat.expanded$heritability)
cor.test(dat.core$var_muth.ascribed, dat.core$heritability)
cor.test(dat.expanded$var_muth.ascribed, dat.expanded$heritability)
cor.test(dat.core$var_uzDG, dat.core$heritability)
cor.test(dat.expanded$var_uzDG, dat.expanded$heritability)
cor.test(dat.core$var_uzDS, dat.core$heritability)
cor.test(dat.expanded$var_uzDS, dat.expanded$heritability)

qplot(var_uzDS, heritability, geom = c("abline"), data=dat.core)

plot1 <- ggplot(data = dat.core, aes(x = var_muth, heritability))
plot2 <- ggplot(data = dat.expanded, aes(x = var_muth, heritability))
plot3 <- ggplot(data = dat.core, aes(x = var_muth.ascribed, heritability))
plot4 <- ggplot(data = dat.expanded, aes(x = var_muth.ascribed, heritability))
plot5 <- ggplot(data = dat.core, aes(x = var_uzDG, heritability))
plot6 <- ggplot(data = dat.expanded, aes(x = var_uzDG, heritability))
plot7 <- ggplot(data = dat.core, aes(x = var_uzDS, heritability))
plot8 <- ggplot(data = dat.expanded, aes(x = var_uzDS, heritability))

plot1 + geom_point() +  geom_smooth(method='lm')
plot2 + geom_point() +  geom_smooth(method='lm')
plot3 + geom_point() +  geom_smooth(method='lm')
plot4 + geom_point() +  geom_smooth(method='lm')
plot5 + geom_point() +  geom_smooth(method='lm')
plot6 + geom_point() +  geom_smooth(method='lm')
plot7 + geom_point() +  geom_smooth(method='lm')
plot8 + geom_point() +  geom_smooth(method='lm')


model.01 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.02 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.03 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.04 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.05 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.06 <- lmer(heritability ~ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.07 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.08 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.09 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.10 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.11 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.12 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.13 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.14 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.15 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.16 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.17 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.18 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.19 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.20 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.21 <- lmer(heritability ~ var_muth.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.22 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.23 <- lmer(heritability ~ var_muth.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.24 <- lmer(heritability ~ var_muth.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.25 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.26 <- lmer(heritability ~ var_muth.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.27 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.28 <- lmer(heritability ~ var_muth.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.29 <- lmer(heritability ~ var_muth.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.30 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.31 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.32 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.33 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.34 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.35 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.36 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.37 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.38 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.39 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.40 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.41 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.42 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.43 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.44 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.45 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.46 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.47 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.48 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.49 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.50 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.51 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.52 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.53 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.54 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.55 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.56 <- lmer(heritability ~ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.57 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.58 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.59 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.60 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.61 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.62 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.63 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.64 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.65 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.66 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.67 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.68 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.69 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.70 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.71 <- lmer(heritability ~ var_muth.ascribed+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.72 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.73 <- lmer(heritability ~ var_muth.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.74 <- lmer(heritability ~ var_muth.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.75 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.76 <- lmer(heritability ~ var_muth.ascribed+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.77 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.78 <- lmer(heritability ~ var_muth.ascribed+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.79 <- lmer(heritability ~ var_muth.ascribed+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.80 <- lmer(heritability ~ var_muth.ascribed+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.81 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.82 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.83 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.84 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.85 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.86 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.87 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.88 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.89 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.90 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.91 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.92 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.93 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.94 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.95 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.96 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.97 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.98 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.99 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.100 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

#-----------

models <- paste0("model.", 1:100)
var.types <- c("Muth index", "Muth index (ascr.)", "Uz DG index", "Uz DS index") 
rsquared.list <- matrix(NA,100,3)

for (i in seq(1,91,by=10)){
  datascope <- ifelse(i<51, "Core",  "Expanded")
  vartype.id <- (((i+9)/10)-1) %% 5 
  if (i %in% c(1,51)){ tablerows <- c("GDP per cap.", "education", "population")
    }else{ tablerows <- c(var.types[vartype.id], "GDP per cap.", "education", "population") }
  for (j in seq(0,9)){rsquared.list[i+j,1:2] <-  round(r.squaredGLMM(get(model.list[i+j])), 4)} # get R-squared
  if (i %in% 11:41){ rsquared.list[i:(i+9),3] <- round(rsquared.list[i:(i+9),1] - rsquared.list[1:10,1],4) } # get r-squared gains compared to models without cultural variance
  if (i %in% 61:91){ rsquared.list[i:(i+9),3] <- round(rsquared.list[i:(i+9),1] - rsquared.list[51:60,1],4)}
  
  stargazer(get(models[i]),get(models[i+1]),get(models[i+2]),get(models[i+3]),get(models[i+4]),
            get(models[i+5]),get(models[i+6]),get(models[i+7]),get(models[i+8]),get(models[i+9]),
            column.labels = c(paste0(datascope," traits at Subchapter-level"), paste0(datascope," traits at Domain-level")), 
            column.separate = c(5, 5), title= paste0("Models ",i," to ", i+9, ": ", datascope, " trait list, ", var.types[vartype.id]), 
            covariate.labels = tablerows, ci = FALSE, align=TRUE,
            add.lines = list(c("R2 marginal", rsquared.list[i:(i+9),1] ),
                             c("R2 conditional", rsquared.list[i:(i+9),2] ),
                             c("R2 marginal gain", rsquared.list[i:(i+9),3] )),
            out= paste0("modelresults_",i,"_to_",i+9, ".html" ))
}

# stack two tables?

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


