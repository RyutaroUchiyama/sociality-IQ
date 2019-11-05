

# In UZ, No var_uzDG for China. Can we find WVS data for 
# If we're reconstructing Uz CTL variable with broader range of WVS surveys, 
#   should we try to match WVS years with years of paper? (probably futile)


#----
library(reshape2)
library(lme4)
library(readxl)

###library(car)
library(ggplot2)
###library(stargazer) 
library(stringr)
library(MuMIn)
library(lmerTest)#????

original_wd <- getwd()
# Import WVS data used for Muthukrishna et al. Cultural Distance analysis
wvs <- read.csv("/Users/ryutaro/Desktop/Muthukrishna_Lab/cultural distance material/wvs_allelesv02_merged.csv")
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

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# From here, code adapted from Muthukrishna et al. Cultural Distance script 
# (from Table1_Calculate_CFst_copy2.R)

# extract relevant subsets of data
include_vars <- as.character(dimensions[dimensions$CAT.All==1,]$V)
wvs <- wvs[,(names(wvs) %in% include_vars | names(wvs) %in% c("S002", "S003", "X048WVS"))]
wvs$pop <- wvs$S003
col_idx <- grep("pop", names(wvs)) #Move pop to the front
wvs <- wvs[, c(col_idx, (1:ncol(wvs))[-col_idx])]
wvs.rep_subset <- wvs[str_detect(wvs$S002,"^2005") | str_detect(wvs$S002,"^2010"),]
wvs.rep_subset[wvs.rep_subset<0] <- NA

loci <- ls(wvs.rep_subset)
loci <- loci[loci != "X048WVS"]
loci <- loci[loci != "pop"]
loci <- loci[loci != "Region"]
loci <- loci[loci != "S002"]
loci <- loci[loci != "S003"]
# what is the type of data, discrete (0) or continuous (1)
type <- c(0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	1, 	1, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0 )
names(type) <- loci # make sure to give names to elements of this vector

# ---  ---- --- --- -- 
# from (CultureFst.r)
pops = as.character(unique(wvs.rep_subset[,1]))
var.loci = function( d, l ){
  # d is the data matrix
  # l is the name of the trait	
  print(paste0("locus = ", l))
  
  # condition on whether it is a discrete or quantitative character
  if( type[[l]]==1 ){ # quantitative character
    print( paste( "q trait", l ) )
    # find total variance
    var.locus = sapply( 1:length(pops), function(y){ yo = pops[y]; var(  d[ d[,1]==yo, l], na.rm=T ) } ) 

  }else{ # discrete character
    # count the unique variants   
    polymorphs = as.character( unique(d[,l]) ); 
    polymorphs = polymorphs[is.na(polymorphs)==FALSE]
    
    # find the frequency of any variant per trait per population
    freq.within = lapply( pops, function(z) { 
      data = d[ d[,1]==z, ]
      freq = rep(0, length(polymorphs) )
      for( i in 1:length(polymorphs) ){
        freq[i] = sum( data[,l]==polymorphs[i], na.rm=T ) / 	sum( is.na(data[,l])==FALSE )
        }
      freq
      } )
    names(freq.within) = pops      
    
    if (length(polymorphs)>0){ 
    var.each.k <- matrix(NA, length(pops), length(polymorphs))
      for( k in 1:length(polymorphs)){
        var.each.k[,k] <- sapply(1:length(pops), function(y){ freq.within[[pops[y]]][k] * (1 - freq.within[[pops[y]]][k]) })
      }
      var.locus <- rowSums(var.each.k) # variance at the current locus
    }else{ # if there is no information about answer frequencies, then fill array with NA 
      var.locus <- matrix(NA, length(pops), 1)    
    }  
    var.locus
  } 
  var.locus
}

# across all traits
var.all.loci = sapply( loci, function(z){ var.loci(wvs.rep_subset, z) } )
var.societies <-  rowMeans(var.all.loci, na.rm=T)
var.societies <- data.frame(pops, "CountryISO"= rep(NA, length(var.societies)) , var.societies) 
var.societies <- var.societies[order(var.societies$pops),]

# attaching the ISO codes only to countries for which we have heritability data 
var.societies$CountryISO[c( 5,12,14,21,22,26,30,34,35,47,50,58,64,66,74)] <-  
  c("AU","CA","CN","FI","FR","GB","IN","IT","JP","NL","NO","RU","KR","SE","US")
var.societies$CountryISO <- as.factor(var.societies$CountryISO)
var.societies.reduced <- var.societies[var.societies$CountryISO %in% country,] # only inlcude countries for which we have heritability data

# Merge Uz 2015 tightness-looseness data with MaTCH data
dat <- merge(uz[uz$CountryISO %in% country,], dat.pre, by.x="CountryISO", by.y="country", all.y = T)
dat <- data.frame(dat[,1:2], matrix(NA,dim(dat)[1],2), dat$CTL_DS, rep(NA,dim(dat)[1]), dat$CTL_DG, rep(NA,dim(dat)[1]), dat$CTL_C, rep(NA,dim(dat)[1]), dat[,-c(1:5)])
colnames(dat)[3:10] <- c("var_muth","var_muth.ascribed","var_uzDS","var_uzDS.ascribed","var_uzDG","var_uzDG.ascribed","var_uzC","var_uzC.ascribed")
dat$var_muth <- var.societies.reduced$var.societies[match(dat$CountryISO, var.societies.reduced$CountryISO)]

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
dat$var_muth.ascribed[which(dat$Country=="Belgium")] <- mean(var.societies$var.societies[which(var.societies$pops %in% c("France","Germany","Netherlands"))])
dat$var_muth.ascribed[which(dat$Country=="Denmark")] <- mean(var.societies$var.societies[which(var.societies$pops %in% c("Sweden","Germany"))])
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
# Adding to the above list all traits within Psychiatric domain except for 'sleep functions'
culturaltraits_withpsychiatric <- union(culturaltraits_nonpsychiatric, setdiff(which(dat.long$trait_domain=="Psychiatric"), which(dat.long$trait_subchapter=="Sleep Functions")))
dat.core <- dat.long[culturaltraits_nonpsychiatric,]
dat.expanded <- dat.long[culturaltraits_withpsychiatric,]

#------
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
plot3 <- ggplot(data = dat.core, aes(x = var_uzDG, heritability))
plot4 <- ggplot(data = dat.expanded, aes(x = var_uzDG, heritability))
plot5 <- ggplot(data = dat.core, aes(x = var_uzDS, heritability))
plot6 <- ggplot(data = dat.expanded, aes(x = var_uzDS, heritability))

plot1 + geom_point() +  geom_smooth(method='lm')
plot2 + geom_point() +  geom_smooth(method='lm')
plot3 + geom_point() +  geom_smooth(method='lm')
plot4 + geom_point() +  geom_smooth(method='lm')
plot5 + geom_point() +  geom_smooth(method='lm')
plot6 + geom_point() +  geom_smooth(method='lm')


model.1 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.2 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.3 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.4 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.5 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.6 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.7 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.8 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.9 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.10 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.11 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.12 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.13 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.14 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.15 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.16 <- lmer(heritability ~ var_muth+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.17 <- lmer(heritability ~ var_muth+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.18 <- lmer(heritability ~ var_muth+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.19 <- lmer(heritability ~ var_muth+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.20 <- lmer(heritability ~ var_muth+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.21 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.22 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.23 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.24 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.25 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.26 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.27 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.28 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.29 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.30 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)

model.31 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.32 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.33 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.34 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.35 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.36 <- lmer(heritability ~ var_uzDG+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.37 <- lmer(heritability ~ var_uzDG+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.38 <- lmer(heritability ~ var_uzDG+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.39 <- lmer(heritability ~ var_uzDG+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.40 <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

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

model.51 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.52 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.53 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.54 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.55 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.expanded)
model.56 <- lmer(heritability ~ var_uzDS+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.57 <- lmer(heritability ~ var_uzDS+ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.58 <- lmer(heritability ~ var_uzDS+ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.59 <- lmer(heritability ~ var_uzDS+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)
model.60 <- lmer(heritability ~ var_uzDS+ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.expanded)

model.61 <- lmer(heritability ~ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.62 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.63 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.64 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.65 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_subchapter)+ (1|country), data=dat.core)
model.66 <- lmer(heritability ~ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.67 <- lmer(heritability ~ e_migdppc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.68 <- lmer(heritability ~ e_peaveduc+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.69 <- lmer(heritability ~ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)
model.70 <- lmer(heritability ~ e_migdppc+ e_peaveduc+ e_wb_pop+ (1|h2_type)+ (1|trait_domain)+ (1|country), data=dat.core)


model.list <- paste0("model.", 1:70)
R2.list <- t(sapply(model.list, function(x)r.squaredGLMM(get(x))))
R2.list <- data.frame(round(R2.list, 3))
colnames(R2.list) <- c("marginal.R2","conditional.R2")
R2.list$ marginal.R2gain <- rep(NA, dim(R2.list)[1])
R2.list$conditional.R2gain <- rep(NA, dim(R2.list)[1])
##R2.list$marginal.R2gain <- R2.list$marginal.R2

R2.list$conditional.R2[1] - R2.list$conditional.R2[61]
R2.list$conditional.R2[1] / R2.list$conditional.R2[61]
R2.list$conditional.R2[5] - R2.list$conditional.R2[65]
R2.list$conditional.R2[5] / R2.list$conditional.R2[65]




baseline <- t(sapply(paste0("model.", 61:65), function(x)r.squaredGLMM(get(x))))
rbind(rep(baseline,12), nrow=nrow(baseline), bycol=T)


matrix( rep( t( baseline ) , 10 ) , ncol = ncol(m) , byrow = TRUE )

rsquared.list$test <- 
# subchapter higher r2 than domain. Core list higher than expanded list. uzDG generally higher than uzDS, but uzDS often higher|Domain



#******
colnames(dat.long)
table(dat.cultural_withpsychiatric$country_name)
table(dat.cultural$trait_chapter)


table(dat.cultural_psychiatric$country_name)

hist(unique(dat.cultural$var_uzDG), breaks=20)


table(dat.long$trait_domain)





dat.long.1 <- dat.long[culturaltraits_nonpsychiatric,]
dat.long.2 <-dat.long[culturaltraits_withpsychiatric,]

tabulate(dat.long.1$e_regionpol)
tabulate(dat.long.2$e_regionpol)

levels(as.factor(dat.long.2$country_name))
tabulate(as.factor(dat.long.2$country_name))



  length(culturaltraits_withpsychiatric)
  
  x <- dat.long[culturaltraits_nonpsychiatric,]
  length(unique(x$country_name))
  
    
    
)



trait_hierarchy2 <- unique( dat.long[,14:12] )
trait_hierarchy2 <- trait_hierarchy2[order(trait_hierarchy2$trait_subchapter) , ]
trait_hierarchy2 <- trait_hierarchy2[order(trait_hierarchy2$trait_chapter) , ]
trait_hierarchy2 <- trait_hierarchy2[order(trait_hierarchy2$trait_domain) , ]
write.csv(trait_hierarchy2, "actual_traits_foranalysis.csv")

trait_hierarchy3 <- unique( dat.long[,12:11] )
trait_hierarchy3 <- trait_hierarchy3[order(trait_hierarchy3$trait_subchapter) , ]
write.csv(trait_hierarchy3, "subchapter_and_trait.csv")




lmm <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ trait+ (1+ trait+ var_uzDG| e_regionpol))



lmm <- lmer(heritability ~ var_uzDG+ e_migdppc+ e_peaveduc+ e_wb_pop+ trait+ )




