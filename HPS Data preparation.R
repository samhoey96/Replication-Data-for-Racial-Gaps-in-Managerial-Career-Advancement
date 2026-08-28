# preamble ----------------------------------------------------------------

#clear environment
rm(list = ls())

# Set this to the folder that contains this script.
# Everything else is resolved relative to it, so the package runs anywhere.
pkgwd <- getwd()   # or: pkgwd <- "/path/to/Replication Package"

# data/source holds the raw inputs this script reads; data/ holds the environments it writes.
# output/ holds the tables and figures the results code writes.
datawd          <- file.path(pkgwd, "data")
inputwd         <- file.path(datawd, "source")
outputwdfiles   <- file.path(pkgwd, "output")
outputwdresults <- outputwdfiles
dir.create(datawd,        showWarnings = FALSE, recursive = TRUE)
dir.create(outputwdfiles, showWarnings = FALSE, recursive = TRUE)

#load packages
library(collapse)
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lfe)
library(xtable)
library(RColorBrewer)
library(survival)
library(Matching)
library(nnet)

setwd(inputwd)

starttime <- Sys.time()

#load data
load("playercareerlevel250423.Rdata")
load("coachdata030523.Rdata")
load("clubdata250423.Rdata")
load("natframe080523.Rdata")
load("leaguedata300623.Rdata")
load("combinedresponses.Rdata")
load("playerids.Rdata")
countries_app <- read.csv("all_clublinks_countries.csv")

#load excel files
jobdescrip <- read.csv("jobdescrip.csv", sep = ";")

# The BFP/DDS race-coding workbooks, playerlevelmergeDS.csv, playerdes.Rdata, clubkeynew.rtf
# and the England betting-odds workbook are no longer read. They fed the dsdata race coding,
# whose only route into the analysis (the block that appended it to responsesperid) is
# commented out, and the odds data was superseded by allgames plus the fixture scrape. Race
# now comes solely from the survey responses in responsesperid.

#unique all data sets to fix duplicate scraping where necessary
playercareerlevel <- unique(playercareerlevel)
coachdata <- unique(coachdata)

#define color palettes
cbp5 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999")

#cbp5 <- c("#039597", "#1ad4d9", "#2c4854", "#2296d6", "#d6e9e9")

cbp6 <- c("#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#CC79A7", "#999999")   

cbp8 <- c("#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999" )

cbp8 <-  brewer.pal(8,"Blues")[3:8]

cbp7 <- c("#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#999999" )

#Custom custom_join function with an alert for increased observations
custom_join <- function(df, df_alt) {
  # Step 1: Check initial number of observations
  initial_observations <- nrow(df)
  
  # Step 2: Perform the custom_join operation
  df <- plyr::join(df, df_alt)
  
  # Step 3: Check the number of observations after custom_joining
  final_observations <- nrow(df)
  
  # Step 4: Show an alert if the number of observations increases
  if (final_observations > initial_observations) {
    message("Alert: Number of observations increased after the custom_join!")
    break
  }
  
  # Step 5: Show an alert if custom_join succesful
  if (final_observations == initial_observations) {
    message("custom_join succesful")
  }
  
  return(df)
}

# definitions ex ante
joborder <- c("coached", "Manager", "Assistant_manager", "Youth_team_management", "Job_outside_UK", "Other_staff") 
insiders <- c("England")
ukcountries <- c("England", "Scotland", "Wales", "Northern Ireland", "Jersey", "Guernsey", "United Kingdom")
internationalleagues <- c(
  "International Friendlies",
  "World Cup qualification",
  "World Cup",
  "UEFA Euro qualifying",
  "EURO",
  "UEFA Nations League",
  "Africa Cup of Nations qualification",
  "Africa Cup of Nations",
  "Confederations Cup",
  "Copa América",
  "Copa América Centenario Play-In",
  "Asian Cup qualification",
  "AFC Asian Cup",
  "CONCACAF Nations League A",
  "CONCACAF Nations League B",
  "CONCACAF Nations League C",
  "CONCACAF Nations League Qualifikation",
  "CONCACAF Cup",
  "Gold Cup 1991",
  "Gold Cup 1993",
  "Gold Cup 1996",
  "Gold Cup 1998",
  "Gold Cup 2000",
  "Gold Cup 2002",
  "Gold Cup 2003",
  "Gold Cup 2005",
  "Gold Cup 2007",
  "Gold Cup 2009",
  "Gold Cup 2011",
  "Gold Cup 2013",
  "Gold Cup 2015",
  "Gold Cup 2017",
  "Gold Cup 2019",
  "Gold Cup 2021",
  "Gold Cup Qualifikation",
  "South Asian Football Federation Championship",
  "East Asian Football Championship",
  "Southeast Asian Championship 2018",
  "Caribbean Cup",
  "Caribbean Cup Qualifikation",
  "African Nations Championship"
)

#change working directory to output working directory
setwd(outputwdresults)

#load fonts
# library(extrafont)
# font_import()
# y

#set sizes
legsize <- 24/1.5
axislabsize <- 24/1.5
textsize <- 28/1.5
ticksize <- 20/1.5
labsize <- 10/1.5
geomtextsize <- 4
legendtextsize <- 20/1.5
facetsize <- 20/1.6

#set race majority cutoff percentage
cutthres <- 0.5

rhs_vars <- c(
  "as.factor(black)",
  "I(sum_minplayed_tier1/60)", "sum_goals_tier1", "sum_assists_tier1",
  "I(sum_minplayed_tier2/60)", "sum_goals_tier2", "sum_assists_tier2",
  "I(sum_minplayed_lowtier/60)", "sum_goals_lowtier", "sum_assists_lowtier",
  "attacker", "midfielder", "defender",
  "as.factor(yob)", "ukplus", "foreign",
  "natteam", "natapps"   # <- added here
)


# DATA PREP 1: country classification into groups --------------------------------------------------

# Given list of countries
countries <- c(
  "Albania", "Algeria", "Angola", "Antigua and Barbuda",
  "Argentina", "Armenia", "Australia", "Austria",
  "Bangladesh", "Barbados", "Belarus", "Belgium",
  "Benin", "Bermuda", "Bolivia", "Bosnia-Herzegovina",
  "Brazil", "British Virgin Islands", "Bulgaria", "Burundi",
  "Cameroon", "Canada", "Cape Verde", "Central African Republic",
  "Chile", "China", "Colombia", "Comoros",
  "Congo", "Costa Rica", "Cote d'Ivoire", "Croatia",
  "Curacao", "Cyprus", "Czech Republic", "Denmark",
  "Dominica", "DR Congo", "Ecuador", "Egypt",
  "England", "Equatorial Guinea", "Estonia", "Eswatini",
  "Faroe Islands", "Finland", "France", "French Guiana",
  "Gabon", "Georgia", "Germany", "Ghana",
  "Gibraltar", "Greece", "Grenada", "Guadeloupe",
  "Guatemala", "Guernsey", "Guinea", "Guinea-Bissau",
  "Guyana", "Haiti", "Honduras", "Hongkong",
  "Hungary", "Iceland", "India", "Iran",
  "Iraq", "Ireland", "Israel", "Italy",
  "Jamaica", "Japan", "Jersey", "Jordan",
  "Jugoslawien (SFR)", "Kenya", "Korea, South", "Kosovo",
  "Latvia", "Lebanon", "Liberia", "Lithuania",
  "Luxembourg", "Malawi", "Mali", "Malta",
  "Martinique", "Mauritius", "Mexico", "Moldova",
  "Montserrat", "Morocco", "Mozambique", "Netherlands",
  "Neukaledonien", "New Zealand", "Nigeria", "North Macedonia",
  "Northern Ireland", "Norway", "Oman", "Pakistan",
  "Paraguay", "Peru", "Philippines", "Poland",
  "Portugal", "Réunion", "Romania", "Russia",
  "Samoa", "Sao Tome and Principe", "Saudi Arabia", "Scotland",
  "Senegal", "Serbia", "Seychelles", "Sierra Leone",
  "Singapore", "Slovakia", "Slovenia", "Somalia",
  "South Africa", "Spain", "Sri Lanka", "St. Kitts & Nevis",
  "St. Lucia", "St. Vincent & Grenadinen", "Suriname", "Sweden",
  "Switzerland", "Tanzania", "The Gambia", "Togo",
  "Trinidad and Tobago", "Tunisia", "Turkey", "Uganda",
  "Ukraine", "United Kingdom", "United States", "Uruguay",
  "Venezuela", "Wales", "Yemen", "Zambia",
  "Zimbabwe"
)

# Create a new variable containing country categories
country_categories <- rep(NA, length(countries))  # Initialize the new variable with NA values

# Assign country categories based on conditions
country_categories[grep("England|Wales|Scotland|Jersey|Guernsey|United Kingdom", countries, ignore.case = TRUE)] <- "UK"
country_categories[grep("Iceland|Croatia|Italy|Luxembourg|Moldova|Serbia|Malta|Slovenia|Ukraine|Turkey|Switzerland|Finland|Germany|Armenia|Belarus|France|Norway|Ireland|Portugal|Poland|Belgium|Romania|Netherlands|Denmark|Cyprus|Sweden|Austria|Latvia|Estonia|Lithuania|Greece|Hungary|Slovakia|Czech Republic|Spain|Albania|Bosnia-Herzegovina|North Macedonia", countries, ignore.case = TRUE)] <- "Europe"
country_categories[grep("Zambia|Tanzania|Central African Republic|Congo|South Africa|Senegal|Cote d'Ivoire|Zimbabwe|Egypt|Morocco|Algeria|Nigeria|Ghana|Sierra Leone|Liberia|Cameroon|Angola|The Gambia|Guinea|Togo|Mali|DR Congo|Cape Verde|Tunisia|Mauritius|Uganda|Burundi|Seychelles|Mozambique|Kenya|Malawi|Somalia", countries, ignore.case = TRUE)] <- "Africa"
country_categories[grep("Northern Ireland", countries, ignore.case = TRUE)] <- "UK"
country_categories[is.na(country_categories)] <- "ROW"  # Assign "Rest of World" to the remaining countries

# put result in data frame
countryclass <- data.frame(nationality = countries, natcat = country_categories)

# DATA PREP 2: compiling survey responses ---------------------------------

surveyresponses <- combined_responses
# table(totalvotes$totalvotes)
table(combined_responses$ethnicity)
table(combined_responses$origin)
table(combined_responses$soccerfan)
table(combined_responses$age)
rm(combined_responses)

#remove unecessary information and plyr::join in link identifier
surveyresponses <- plyr::join(surveyresponses, ids)
responsesperid <- surveyresponses[,c("id", "link", "Response")]

#add bfp codings
# dsdata <- plyr::join(dsdata, ids)
# dsdata <- dsdata[,c("id", "link", "black")]
# rm(ids)
# dsdata$Response <- NA
# dsdata$Response <- ifelse(dsdata$black == 0, "Any other background/ethnicity", dsdata$Response)
# dsdata$Response <- ifelse(dsdata$black == 1, "Black / African / Caribbean / African American/ Mixed Black Heritage/ Any other Black background", dsdata$Response)
# dsdata$Response <- ifelse(dsdata$black == 2, "Don't know", dsdata$Response)
# dsdata <- dsdata[,-match("black", colnames(dsdata))]
# dsdata$id <- as.character(dsdata$id)
# dsdata <- subset.data.frame(dsdata, !is.na(dsdata$Response))
# responsesperid <- rbind(responsesperid, dsdata)

#count categorizations per id
responsesperid <- responsesperid %>%
  group_by(id, link, Response) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Response, values_from = count, values_fill = 0)
responsesperid$Other <- responsesperid$Other + responsesperid$`Any other background/ethnicity`
responsesperid <- responsesperid[,-match("Any other background/ethnicity", colnames(responsesperid))]
responsesperid <- subset.data.frame(responsesperid, !is.na(id))
colnames(responsesperid) <- c("id", "link", "other", "dn", "black")
responsesperid$total <- responsesperid$other + responsesperid$black + responsesperid$dn
responsesperid$id <- as.numeric(responsesperid$id)
responsesperid <- subset.data.frame(responsesperid, !is.na(id))
table(responsesperid$total)

#make disagreement variable
responsesperid$disagree <- ifelse(responsesperid$other > 0 & responsesperid$black > 0, 1, 0)
responsesperid$disagree <- ifelse(responsesperid$total < 2, NA, responsesperid$disagree)

#pick majority vote if possible
responsesperid$othershare <- responsesperid$other/responsesperid$total
responsesperid$blackshare <- responsesperid$black/responsesperid$total
responsesperid$othermaj <- ifelse(responsesperid$othershare > cutthres,1,0)
responsesperid$blackmaj <- ifelse(responsesperid$blackshare > cutthres,1,0)
responsesperid$nomaj <- ifelse(responsesperid$othermaj == 0 & responsesperid$blackmaj == 0, 1, 0)

#statistics
table(responsesperid$othermaj)
table(responsesperid$blackmaj)
table(responsesperid$nomaj)
table(responsesperid$disagree)
table(responsesperid$total)

# DATA PREP 3: compiling player season level data -------------------------------------------------------

#make relevant vars numeric
playercareerlevel$appearances <- as.numeric(playercareerlevel$appearances)
playercareerlevel$goals <- as.numeric(playercareerlevel$goals)
playercareerlevel$assists <- as.numeric(playercareerlevel$assists)
playercareerlevel$minplayed <- as.numeric(str_remove_all(playercareerlevel$minplayed, "\\D"))

#fix season coding
playercareerlevel$season <- ifelse(str_detect(playercareerlevel$season, "/"), substr(playercareerlevel$season, 4,5), playercareerlevel$season)
playercareerlevel$season <- ifelse(nchar(playercareerlevel$season) == 2 & as.numeric(playercareerlevel$season) < 25, paste0("20", playercareerlevel$season), ifelse(nchar(playercareerlevel$season) == 2 & as.numeric(playercareerlevel$season) >= 25, paste0("19", playercareerlevel$season), playercareerlevel$season))
playercareerlevel$season <- as.numeric(playercareerlevel$season)
unique(playercareerlevel$season)

#find first and last season or each player in the data
minmax <- collap(playercareerlevel, season ~ link + dob + playername, FUN = c(min, max), na.rm = T)
minmax <- minmax %>%
  mutate_all(~ replace(., is.infinite(.), NA))

#merge into playercareerlevel
playercareerlevel <- custom_join(playercareerlevel, minmax)

#fix season variable
clubdata$season <- ifelse(str_detect(clubdata$season, "/"), substr(clubdata$season, 4,5), clubdata$season)
clubdata$season <- ifelse(nchar(clubdata$season) == 2 & as.numeric(clubdata$season) < 25, 
                          paste0("20", clubdata$season), ifelse(nchar(clubdata$season) == 2 & as.numeric(clubdata$season) >= 25, 
                                                                paste0("19", clubdata$season), 
                                                                clubdata$season))
clubdata$season <- as.numeric(clubdata$season)
test <- subset.data.frame(clubdata, is.na(season))
rm(test)

#make numeric tier variable
table(clubdata$tier)
clubdata$tier_alt <- NA
clubdata$tier_alt <- ifelse(clubdata$tier == "First Tier", 1, clubdata$tier_alt)
clubdata$tier_alt <- ifelse(clubdata$tier == "Second Tier", 2, clubdata$tier_alt)
clubdata$tier_alt <- ifelse(clubdata$tier == "Third Tier", 3, clubdata$tier_alt)
clubdata$tier_alt <- ifelse(clubdata$tier == "Fourth Tier", 4, clubdata$tier_alt)
clubdata$tier_alt <- ifelse(clubdata$tier == "Fifth Tier", 5, clubdata$tier_alt)
clubdata$tier_alt <- ifelse(clubdata$tier == "Sixth Tier", 6, clubdata$tier_alt)

#some fixes for both clubdata and playercareerlevel
complink <- clubdata$complink
complink <- str_split(complink, "/saison_id/")
clubdata$complink <- unlist(lapply(complink, function(x) x[1]))
clublink <- playercareerlevel$clublink_alt
clublink <- str_split(clublink, "/")
playercareerlevel$clublink <- unlist(lapply(clublink, function(x) x[2]))
rm(clublink, complink)

#countries
countries <- leaguedata[,c(1,3)]
colnames(countries)[2] <- "leaguelink"
countries$dups <- ifelse(countries$leaguelink %in% countries$leaguelink[duplicated(countries$leaguelink)],1,0)
countries <- countries[,c(1,2)]

playercareerlevel <- custom_join(playercareerlevel, countries)
rm(countries)

#tiers by competition link
complinks <- unique(clubdata[, c("complink", "tier", "tier_alt")])
complinks$dups <- ifelse(complinks$complink %in% complinks$complink[duplicated(complinks$complink)], 1, 0)
test <- subset.data.frame(complinks, dups == 1)
rm(test) 

#keep only lowest (in numbers) tier
complinks <- complinks %>% 
  group_by(complink) %>% 
  filter(tier == min(tier))

playercareerlevel <- custom_join(playercareerlevel, complinks)
test <- playercareerlevel[,c("season", "link", "clublink", "complink", "country", "tier")]
rm(test)

#test which competitions did not get a tier assigned
test <- unique(subset.data.frame(playercareerlevel, is.na(tier) & country == "England"))
unique(test$competition) #they all seem to not be top 4 tiers
rm(test)

#international caps
playercareerlevel$natcomp <- ifelse(playercareerlevel$competition %in% internationalleagues, 1, 0)
playercareerlevel <- playercareerlevel %>%
  group_by(link) %>%
  mutate(natapps = sum(appearances * natcomp, na.rm = TRUE))
playercareerlevel$natteam <- ifelse(playercareerlevel$natapps > 0, 1, 0)
playercareerlevel$england <- ifelse(playercareerlevel$club %in% ukcountries, 1, 0)
playercareerlevel <- playercareerlevel %>%
  group_by(link) %>%
  mutate(natappseng = sum(appearances * england, na.rm = TRUE))
playercareerlevel$natteameng <- ifelse(playercareerlevel$natappseng > 0, 1, 0)

# DATA PREP 4: compiling coach data ---------------------------------------------------------

test <- subset.data.frame(coachdata, is.na(func)) #all of them truly have no descriptions
rm(test)

#manual edits
coachdata$func <- str_remove(coachdata$func, "Al-Ahli \\(UAE\\)")
coachdata$club <- ifelse(str_detect(coachdata$clubfunc, "Al-Ahli \\(UAE\\)"), "Al-Ahli (UAE)", coachdata$club)
coachdata$func <- str_remove(coachdata$func, "São Bento \\(SP\\)")
coachdata$club <- ifelse(str_detect(coachdata$clubfunc, "São Bento \\(SP\\)"), "São Bento (SP)", coachdata$club)
coachdata$func <- str_remove(coachdata$func, "San Martín \\(SJ\\)")
coachdata$club <- ifelse(str_detect(coachdata$clubfunc, "San Martín \\(SJ\\)"), "San Martín (SJ)", coachdata$club)
coachdata$func <- str_remove(coachdata$func, "Víkingur Ó\\.")
coachdata$club <- ifelse(str_detect(coachdata$clubfunc, "Víkingur Ó\\."), "Víkingur Ó.", coachdata$club)
coachdata$func <- str_remove(coachdata$func, "University \\(J\\)")
coachdata$club <- ifelse(str_detect(coachdata$clubfunc, "University \\(J\\)"), "University (J)", coachdata$club)

#save and edit old job descriptions
# functions <- functions[,-1]
# colnames(functions) <- c("func", "freq", "funccat")
# coachdata$func <- trimws(coachdata$func)

# coachdata <- custom_join(coachdata, functions)
# test <- subset.data.frame(coachdata, is.na(func))

# jobdescrip <- unique(coachdata[,c("func", "funccat")])
#write.csv(jobdescrip, "jobdescrip.csv")

#load and merge new job descriptions
jobdescrip <- jobdescrip[,-1]
coachdata$func <- trimws(coachdata$func)
coachdata <- custom_join(coachdata, jobdescrip)
coachdata$funccat <- ifelse(coachdata$funccat %in% c("Medical staff", "Other support staff", "Player coach", "Executive function", "Scout"), "Other staff", coachdata$funccat)
coachdata$funccat <- ifelse(coachdata$funccat %in% c("Lower Manager"), "Assistant manager", coachdata$funccat)
table(coachdata$funccat)
rm(jobdescrip)

#recode coachdata start and end seasons
coachdata$startspell <- ifelse(str_detect(coachdata$start, "/"), substr(coachdata$start, 4,5), coachdata$start)
coachdata$startspell <- ifelse(nchar(coachdata$startspell) == 2 & as.numeric(coachdata$startspell) < 25, paste0("20", coachdata$startspell), ifelse(nchar(coachdata$startspell) == 2 & as.numeric(coachdata$startspell) >= 25, paste0("19", coachdata$startspell), coachdata$startspell))
coachdata$startspell <- as.numeric(coachdata$startspell)
table(coachdata$startspell)
coachdata$endspell <- ifelse(str_detect(coachdata$end, "/"), substr(coachdata$end, 4,5), coachdata$end)
coachdata$endspell <- ifelse(substr(coachdata$end,5,5) == "/", substr(coachdata$end, 6, 7), coachdata$endspell)
coachdata$endspell <- ifelse(nchar(coachdata$endspell) == 2 & as.numeric(coachdata$endspell) < 25, paste0("20", coachdata$endspell), ifelse(nchar(coachdata$endspell) == 2 & as.numeric(coachdata$endspell) >= 25, paste0("19", coachdata$endspell), coachdata$endspell))
coachdata$endspell <- as.numeric(coachdata$endspell)
table(coachdata$endspell)

#some tests on startspell and endspell
teststart <- subset.data.frame(coachdata, is.na(startspell)) #no start dates available for these guys
testend <- subset.data.frame(coachdata, is.na(endspell)) #all not finished yet it seems
sum(str_detect(testend$end, "expected"), na.rm = T)
sum(is.na(testend$end))
sum(str_detect(testend$end, "expected"), na.rm = T) + sum(is.na(testend$end))
test <- subset.data.frame(testend, !(is.na(testend$end) | str_detect(testend$end, "expected")))
rm(teststart, testend, test)

#figuring out which club is which
splitout <- str_split(coachdata$clublink, "-")
splitout <- sapply(splitout, tail, 1)
splitout <- subset(splitout, !is.na(splitout))
splitout <- data.frame(table(splitout))

splitout <- str_split(coachdata$clublink, "-")
splitout <- sapply(splitout, tail, 1)
coachdata$lastterm  <- sapply(splitout, tail, 1)
coachdata$youth <- ifelse(coachdata$lastterm %in% c(unique(splitout[grep("^u\\d{2}$", splitout)]), "jugend", "res.", "Res.", "reserves", "ii", "b", "c", "2"), 1, 0)
coachdata$funccat <- ifelse(coachdata$youth == 1 & coachdata$funccat %in% c("Manager", "Lower Manager"), "Youth team management/development", coachdata$funccat)

#short last terms
shortlast <- sort(unique(subset(coachdata$lastterm, nchar(coachdata$lastterm) < 5 & coachdata$youth == 0)))
shortlast
rm(shortlast)

#merge in first and last season as player
coachdata <- custom_join(coachdata, minmax)

#fix coach licence
licenses <- str_split(coachdata$coachlicence, "Preferred formation")
coachdata$coachlicence <-  sapply(licenses, head, 1)
licenses <- str_split(coachdata$coachlicence, "Website")
coachdata$coachlicence <-  sapply(licenses, head, 1)
table(coachdata$coachlicence)
test <- subset.data.frame(coachdata, is.na(coachlicence)) #coach licence only missing when no history available it seems
rm(test)
coachdata$coachlicence <- ifelse(coachdata$coachlicence == "missing", NA, coachdata$coachlicence)

#add countries
clubs <- unique(clubdata[,c("clublink", "complink")])
countries <- leaguedata[,c(1,3)]
colnames(countries)[2] <- "complink"
clubs <- custom_join(clubs, countries)
countries <- unique(clubs[,c("clublink", "country")])
special_info_vector <- c(unique(coachdata$lastterm[grep("^u\\d{2}$", coachdata$lastterm)]), "jugend", "res.", "Res.", "reserves", "ii", "b", "c", "2", "uefa")
for (special_info in special_info_vector){
  pattern <- paste0("-", special_info, "$")
  coachdata$clublink <- sub(pattern, "", coachdata$clublink)
}
countries$dups <- ifelse(countries$clublink %in% countries$clublink[duplicated(countries$clublink)],1,0)
dups <- subset.data.frame(countries, dups == 1)
unique(dups$country)
countries <- subset.data.frame(countries, !duplicated(clublink))
countries <- countries[,c(1,2)]

coachdata <- custom_join(coachdata, countries_app)
test <- subset.data.frame(coachdata, is.na(country)) #no league information available
rm(test)

coachdata$funccat_alt <- coachdata$funccat
coachdata$funccat <- ifelse(!(coachdata$country %in% c("England")) & coachdata$funccat != "No Job" & coachdata$funccat != "Other staff", "Job outside UK", coachdata$funccat)
coachdata$funccat_bn <- coachdata$funccat_alt
coachdata$funccat_bn <- ifelse(coachdata$funccat_bn == "Manager" & coachdata$country == "England", "Manager England", coachdata$funccat_bn)
coachdata$funccat_bn <- ifelse(coachdata$funccat_bn == "Manager" & coachdata$country != "England", "Manager ROW", coachdata$funccat_bn)
coachdata$funccat_bn <- ifelse(coachdata$funccat_bn == "Assistant manager" & coachdata$country == "England", "Assistant manager England", coachdata$funccat_bn)
coachdata$funccat_bn <- ifelse(coachdata$funccat_bn == "Assistant manager" & coachdata$country != "England", "Assistant manager ROW", coachdata$funccat_bn)
coachdata$funccat_bn <- ifelse(coachdata$funccat_bn == "Youth team management/development" & coachdata$country == "England", "Youth England", coachdata$funccat_bn)
coachdata$funccat_bn <- ifelse(coachdata$funccat_bn == "Youth team management/development" & coachdata$country != "England", "Youth ROW", coachdata$funccat_bn)
table(coachdata$funccat)
table(coachdata$funccat_alt)
table(coachdata$funccat_bn)

# #who has done which job?
coachdata$jobdone <- 1
jobdone <- unique(coachdata[,c("link", "funccat", "jobdone")])
jobdone <- spread(jobdone, "funccat", "jobdone")
jobdone[is.na(jobdone)] <- 0
jobdone <- jobdone[,-c(match("<NA>", colnames(jobdone)))]

#who has done which job alt?
coachdata$jobdone_alt <- 1
jobdone_alt <- unique(coachdata[,c("link", "funccat_alt", "jobdone_alt")])
jobdone_alt <- spread(jobdone_alt, "funccat_alt", "jobdone_alt")
jobdone_alt[is.na(jobdone_alt)] <- 0
jobdone_alt <- jobdone_alt[,-c(match("<NA>", colnames(jobdone_alt)))]

#who has done which job bn?
coachdata$jobdone_bn <- 1
jobdone_bn <- unique(coachdata[,c("link", "funccat_bn", "jobdone_bn")])
jobdone_bn <- spread(jobdone_bn, "funccat_bn", "jobdone_bn")
jobdone_bn[is.na(jobdone_bn)] <- 0
jobdone_bn <- jobdone_bn[,-c(match("<NA>", colnames(jobdone_bn)))]

#subset on spells where we know exact start and enddate
coachdata <- subset.data.frame(coachdata, !is.na(startspell) & nchar(start) >= 5 & nchar(end) >= 5)

#look at obs that were dropped
rmtest <- subset.data.frame(coachdata, !(!is.na(startspell) & nchar(start) >= 5 & nchar(end) >= 5))
rmlinks <- unique(rmtest$link)
rm(rmtest)

#make charcounter for start and end
coachdata$startchar <- nchar(coachdata$start)
coachdata$endchar <- nchar(coachdata$end)

#adapt when nchar is lower than 5
coachdata$start <- ifelse(coachdata$startchar <= 5, paste0(coachdata$start, " (Jul 1, ", coachdata$startspell-1, ")"), coachdata$start)
coachdata$end <- ifelse(coachdata$endchar <= 5, paste0(coachdata$end, " (Jun 30, ", coachdata$endspell, ")"), coachdata$end)

#startdate
#extract the date using regular expressions
startdate <- regmatches(coachdata$start, regexpr("\\((.*?)\\)", coachdata$start))
startdate <- gsub("\\(|\\)", "", startdate)

#convert the date string to a Date object
startdate <-  as.Date(startdate, format = "%b %d, %Y")
coachdata$startdate <- startdate
#check
test <- coachdata[,c("start", "startdate")] #looks good
rm(test)

#end date
coachdata$end <- ifelse(str_detect(coachdata$end, "expected"), "(Jan 1, 2024)", coachdata$end)
#extract the date using regular expressions
enddate <- regmatches(coachdata$end, regexpr("\\((.*?)\\)", coachdata$end))
enddate <- gsub("\\(|\\)", "", enddate)

#convert the date string to a Date object
enddate <-  as.Date(enddate, format = "%b %d, %Y")
coachdata$enddate <- enddate
#check
test <- coachdata[,c("end", "enddate")] #looks good
rm(test)

#convert start and enddate to numeric 
coachdata$startdatenum <- as.numeric(coachdata$startdate)
coachdata$enddatenum <- as.numeric(coachdata$enddate)

#sort the dataframe by coach identifier and start date
coachdata <- coachdata %>% 
  arrange(link, startdatenum)

#add countries
clubs <- unique(clubdata[,c("clublink", "complink")])
countries <- leaguedata[,c(1,3)]
colnames(countries)[2] <- "complink"
clubs <- custom_join(clubs, countries)
countries <- unique(clubs[,c("clublink", "country")])
special_info_vector <- c(unique(coachdata$lastterm[grep("^u\\d{2}$", coachdata$lastterm)]), "jugend", "res.", "Res.", "reserves", "ii", "b", "c", "2")
for (special_info in special_info_vector){
  pattern <- paste0("-", special_info, "$")
  coachdata$clublink <- sub(pattern, "", coachdata$clublink)
}
countries$dups <- ifelse(countries$clublink %in% countries$clublink[duplicated(countries$clublink)],1,0)
dups <- subset.data.frame(countries, dups == 1)
unique(dups$country)
countries <- subset.data.frame(countries, !duplicated(clublink))
countries <- countries[,c(1,2)]

coachdata <- custom_join(coachdata, countries_app)
test <- subset.data.frame(coachdata, is.na(country)) #no league information available
rm(test)

#licence dummy
table(coachdata$coachlicence)
coachdata$licencedum <- ifelse(!is.na(coachdata$coachlicence), 1, 0)

#funccat table
funccattable <- data.frame(table(coachdata$func, coachdata$funccat))
funccattable <- spread(funccattable, Var2, Freq)
funccattable <- funccattable[,c(1,3,5,6,2,4)]
funccattable$Total <- round(rowSums(funccattable[,c(2:6)]),0)
funccattable <- funccattable[order(-funccattable$Total),]
funccattable$Total <- str_remove(as.character(funccattable$Total), ".00")
colnames(funccattable)[1] <- "Job Description"


print(xtable(funccattable), include.rownames = F)

for(i in 2:7){
  funccattable[,i] <- as.numeric(funccattable[,i])
}
colSums(funccattable[,c(2:7)])

# DATA PREP 5: compile player level data for those with tier 1-4 experience England -----------------------------------------------

test <- subset.data.frame(playercareerlevel, country == "England" & tier %in% c("First Tier", "Second Tier", "Third Tier", "Fourth Tier", "Fifth Tier", "Sixth Tier"))
test <- subset.data.frame(playercareerlevel, country == "England" & is.na(tier))
table(test$leaguelink) #none of them are non cup/qualification/knock out leagues.
table(test$competition) #none of them are non cup/qualification/knock out leagues.
rm(test)

#subset on correct leagues and tiers and only playerseasons with more than 0 appearances
combinedplayercareer <- subset.data.frame(playercareerlevel, country == "England" & tier %in% c("First Tier", "Second Tier", "Third Tier", "Fourth Tier", "Fifth Tier", "Sixth Tier") & appearances > 0)

#define career statistics within player-tiers
combinedplayercareer <- combinedplayercareer %>%
  group_by(link) %>%
  mutate(
    sum_appearances_tier1 = sum(appearances[tier_alt == 1], na.rm = TRUE),
    sum_goals_tier1 = sum(goals[tier_alt == 1], na.rm = TRUE),
    sum_assists_tier1 = sum(assists[tier_alt == 1], na.rm = TRUE),
    sum_minplayed_tier1 = sum(minplayed[tier_alt == 1], na.rm = TRUE),
    
    sum_appearances_tier2 = sum(appearances[tier_alt == 2], na.rm = TRUE),
    sum_goals_tier2 = sum(goals[tier_alt == 2], na.rm = TRUE),
    sum_assists_tier2 = sum(assists[tier_alt == 2], na.rm = TRUE),
    sum_minplayed_tier2 = sum(minplayed[tier_alt == 2], na.rm = TRUE),
    
    sum_appearances_tier3 = sum(appearances[tier_alt == 3], na.rm = TRUE),
    sum_goals_tier3 = sum(goals[tier_alt == 3], na.rm = TRUE),
    sum_assists_tier3 = sum(assists[tier_alt == 3], na.rm = TRUE),
    sum_minplayed_tier3 = sum(minplayed[tier_alt == 3], na.rm = TRUE),
    
    sum_appearances_tier4 = sum(appearances[tier_alt == 4], na.rm = TRUE),
    sum_goals_tier4 = sum(goals[tier_alt == 4], na.rm = TRUE),
    sum_assists_tier4 = sum(assists[tier_alt == 4], na.rm = TRUE),
    sum_minplayed_tier4 = sum(minplayed[tier_alt == 4], na.rm = TRUE),
    
    sum_appearances_tier5 = sum(appearances[tier_alt == 5], na.rm = TRUE),
    sum_goals_tier5 = sum(goals[tier_alt == 5], na.rm = TRUE),
    sum_assists_tier5 = sum(assists[tier_alt == 5], na.rm = TRUE),
    sum_minplayed_tier5 = sum(minplayed[tier_alt == 5], na.rm = TRUE),
    
    sum_appearances_tier6 = sum(appearances[tier_alt == 6], na.rm = TRUE),
    sum_goals_tier6 = sum(goals[tier_alt == 6], na.rm = TRUE),
    sum_assists_tier6 = sum(assists[tier_alt == 6], na.rm = TRUE),
    sum_minplayed_tier6 = sum(minplayed[tier_alt == 6], na.rm = TRUE)
  )

#create lower tier summed stats
combinedplayercareer$sum_appearances_lowtier <- combinedplayercareer$sum_appearances_tier3 +   combinedplayercareer$sum_appearances_tier4 +  combinedplayercareer$sum_appearances_tier5 +  combinedplayercareer$sum_appearances_tier6
combinedplayercareer$sum_goals_lowtier <- combinedplayercareer$sum_goals_tier3 +   combinedplayercareer$sum_goals_tier4 +  combinedplayercareer$sum_goals_tier5 +  combinedplayercareer$sum_goals_tier6
combinedplayercareer$sum_assists_lowtier <- combinedplayercareer$sum_assists_tier3 +   combinedplayercareer$sum_assists_tier4 +  combinedplayercareer$sum_assists_tier5 +  combinedplayercareer$sum_assists_tier6
combinedplayercareer$sum_minplayed_lowtier <- combinedplayercareer$sum_minplayed_tier3 +   combinedplayercareer$sum_minplayed_tier4 +  combinedplayercareer$sum_minplayed_tier5 +  combinedplayercareer$sum_minplayed_tier6


combinedplayercareer <- as.data.frame(combinedplayercareer)

#collapse to player level
playerlevel <- unique(combinedplayercareer[,c("playername", "link", "dob", "position", 
                                              "sum_appearances_tier1", "sum_goals_tier1", "sum_assists_tier1", "sum_minplayed_tier1", 
                                              "sum_appearances_tier2", "sum_goals_tier2", "sum_assists_tier2", "sum_minplayed_tier2",
                                              "sum_appearances_tier3", "sum_goals_tier3", "sum_assists_tier3", "sum_minplayed_tier3",
                                              "sum_appearances_tier4", "sum_goals_tier4", "sum_assists_tier4", "sum_minplayed_tier4",
                                              "sum_appearances_tier5", "sum_goals_tier5", "sum_assists_tier5", "sum_minplayed_tier5",
                                              "sum_appearances_tier6", "sum_goals_tier6", "sum_assists_tier6", "sum_minplayed_tier6",
                                              "sum_appearances_lowtier", "sum_goals_lowtier", "sum_assists_lowtier", "sum_minplayed_lowtier", 
                                              "natteam", "natapps", "natteameng", "natappseng")]
)
#merge in first and last season as a player
playerlevel <- custom_join(playerlevel, minmax)

#played in tier
playerlevel$onetier <- ifelse(playerlevel$sum_appearances_tier1 > 0, 1, 0)
playerlevel$twotier <- ifelse(playerlevel$sum_appearances_tier2 > 0, 1, 0)
playerlevel$threetier <- ifelse(playerlevel$sum_appearances_tier3 > 0, 1, 0)
playerlevel$fourtier <- ifelse(playerlevel$sum_appearances_tier4 > 0, 1, 0)
# playerlevel$fifthtier <- ifelse(playerlevel$sum_appearances_tier5 > 0, 1, 0)
# playerlevel$sixthtier <- ifelse(playerlevel$sum_appearances_tier6 > 0, 1, 0)

#recode position
sum(is.na(playerlevel$position))
position <- sort(unique(playerlevel$position))
position_alt <- c("A", "M", "M", "D", "A", "D", "M", "G", "M", "A", "D", "M", "M", "A", "D", "A", "D")
positions <- data.frame(position, position_alt) #check this
playerlevel <- custom_join(playerlevel, positions)
playerlevel$attacker <- ifelse(playerlevel$position_alt == "A", 1, 0)
playerlevel$defender <- ifelse(playerlevel$position_alt == "D", 1, 0)
playerlevel$midfielder <- ifelse(playerlevel$position_alt == "M", 1, 0)
playerlevel$goalkeeper  <- ifelse(playerlevel$position_alt == "G", 1, 0)
rm(position, positions, position_alt)

#add in coach jobs
# colnames(jobdone) <- str_replace_all(colnames(jobdone), "\\ ", "_")
# colnames(jobdone) <- str_replace_all(colnames(jobdone), "/development", "")
# playerlevel <- custom_join(playerlevel, jobdone)
# jobs <- colnames(jobdone)[2:length(colnames(jobdone))]
# playerlevel[jobs] <-
#   lapply(playerlevel[jobs], function(x) ifelse(is.na(x), 0, x))

colnames(jobdone_alt) <- str_replace_all(colnames(jobdone_alt), "\\ ", "_")
colnames(jobdone_alt) <- str_replace_all(colnames(jobdone_alt), "/development", "")
playerlevel <- custom_join(playerlevel, jobdone_alt)
jobs <- colnames(jobdone_alt)[2:length(colnames(jobdone_alt))]
playerlevel[jobs] <-
  lapply(playerlevel[jobs], function(x) ifelse(is.na(x), 0, x))

colnames(jobdone_bn) <- str_replace_all(colnames(jobdone_bn), "\\ ", "_")
colnames(jobdone_bn) <- str_replace_all(colnames(jobdone_bn), "/development", "")
playerlevel <- custom_join(playerlevel, jobdone_bn)
jobs <- colnames(jobdone_bn)[2:length(colnames(jobdone_bn))]
playerlevel[jobs] <-
  lapply(playerlevel[jobs], function(x) ifelse(is.na(x), 0, x))

#more performance statistics
playerlevel$avg_goalspergame_tier1 <- ifelse(playerlevel$sum_appearances_tier1 > 0, playerlevel$sum_goals_tier1/playerlevel$sum_appearances_tier1, NA)
playerlevel$avg_assistspergame_tier1 <- ifelse(playerlevel$sum_appearances_tier1 > 0,playerlevel$sum_assists_tier1/playerlevel$sum_appearances_tier1, NA)
playerlevel$avg_minplayedpergame_tier1 <- ifelse(playerlevel$sum_appearances_tier1 > 0,playerlevel$sum_minplayed_tier1/playerlevel$sum_appearances_tier1, NA)
playerlevel$avg_goalspergame_tier2 <- ifelse(playerlevel$sum_appearances_tier2 > 0,playerlevel$sum_goals_tier2/playerlevel$sum_appearances_tier2, NA)
playerlevel$avg_assistspergame_tier2 <- ifelse(playerlevel$sum_appearances_tier2 > 0,playerlevel$sum_assists_tier2/playerlevel$sum_appearances_tier2, NA)
playerlevel$avg_minplayedpergame_tier2 <- ifelse(playerlevel$sum_appearances_tier2 > 0,playerlevel$sum_minplayed_tier2/playerlevel$sum_appearances_tier2, NA)
playerlevel$avg_goalspergame_tier3 <- ifelse(playerlevel$sum_appearances_tier3 > 0,playerlevel$sum_goals_tier/playerlevel$sum_appearances_tier3, NA)
playerlevel$avg_assistspergame_tier3 <- ifelse(playerlevel$sum_appearances_tier3 > 0,playerlevel$sum_assists_tier3/playerlevel$sum_appearances_tier3, NA)
playerlevel$avg_minplayedpergame_tier3 <- ifelse(playerlevel$sum_appearances_tier3 > 0,playerlevel$sum_minplayed_tier3/playerlevel$sum_appearances_tier3, NA)
playerlevel$avg_goalspergame_tier4 <- ifelse(playerlevel$sum_appearances_tier4 > 0,playerlevel$sum_goals_tier4/playerlevel$sum_appearances_tier4, NA)
playerlevel$avg_assistspergame_tier4 <- ifelse(playerlevel$sum_appearances_tier4 > 0,playerlevel$sum_assists_tier4/playerlevel$sum_appearances_tier4, NA)
playerlevel$avg_minplayedpergame_tier4 <- ifelse(playerlevel$sum_appearances_tier4 > 0,playerlevel$sum_minplayed_tier4/playerlevel$sum_appearances_tier4, NA)
playerlevel$avg_goalspergame_tier5 <- ifelse(playerlevel$sum_appearances_tier5 > 0,playerlevel$sum_goals_tier5/playerlevel$sum_appearances_tier5, NA)
playerlevel$avg_assistspergame_tier5 <- ifelse(playerlevel$sum_appearances_tier5 > 0,playerlevel$sum_assists_tier5/playerlevel$sum_appearances_tier5, NA)
playerlevel$avg_minplayedpergame_tier5 <- ifelse(playerlevel$sum_appearances_tier5 > 0,playerlevel$sum_minplayed_tier5/playerlevel$sum_appearances_tier5, NA)
playerlevel$avg_goalspergame_tier6 <- ifelse(playerlevel$sum_appearances_tier6 > 0,playerlevel$sum_goals_tier6/playerlevel$sum_appearances_tier6, NA)
playerlevel$avg_assistspergame_tier6 <- ifelse(playerlevel$sum_appearances_tier6 > 0,playerlevel$sum_assists_tier6/playerlevel$sum_appearances_tier6, NA)
playerlevel$avg_minplayedpergame_tier6 <- ifelse(playerlevel$sum_appearances_tier6 > 0,playerlevel$sum_minplayed_tier6/playerlevel$sum_appearances_tier6, NA)
playerlevel$avg_goalspergame_lowtier <- ifelse(playerlevel$sum_appearances_lowtier > 0,playerlevel$sum_goals_lowtier/playerlevel$sum_appearances_lowtier, NA)
playerlevel$avg_assistspergame_lowtier <- ifelse(playerlevel$sum_appearances_lowtier > 0,playerlevel$sum_assists_lowtier/playerlevel$sum_appearances_lowtier, NA)
playerlevel$avg_minplayedpergame_lowtier <- ifelse(playerlevel$sum_appearances_lowtier > 0,playerlevel$sum_minplayed_lowtier/playerlevel$sum_appearances_lowtier, NA)

#add in nationalities
playerlevel <- custom_join(playerlevel, natframe)
playerlevel <- custom_join(playerlevel, countryclass)
colnames(countryclass) <- c("nationality2", "natcat2")
playerlevel <- custom_join(playerlevel, countryclass)
playerlevel$natcat2 <- ifelse(is.na(playerlevel$nationality2), NA, playerlevel$natcat2)

#spread the nationalities across columns
natcats <- unique(playerlevel[,c("link", "natcat")])
natcats$dum <- 1
natcats <- spread(natcats, natcat, -c(1))
natcats[is.na(natcats)] <- 0
playerlevel <- custom_join(playerlevel, natcats)

#adapt for second nationalities
playerlevel$natcat2 <- ifelse(is.na(playerlevel$natcat2), "missing", playerlevel$natcat2)
playerlevel$Africa <- ifelse(playerlevel$natcat2 == "Africa", 1, playerlevel$Africa)
playerlevel$UK <- ifelse(playerlevel$natcat2 == "UK", 1, playerlevel$UK)
playerlevel$ROW <- ifelse(playerlevel$natcat2 == "ROW", 1, playerlevel$ROW)
playerlevel$Europe <- ifelse(playerlevel$natcat2 == "Europe", 1, playerlevel$Europe)

# DATA PREP 6: compile previous and next job data ----------------------------------------------

#create data about previous and next jobs
jobevol <- unique(coachdata[,c("playername", "link", "club", "clublink", "country", "funccat", "startdate", "startdatenum", "enddate", "enddatenum", "max.season")])
jobevol$max.season <- as.numeric(as.Date(paste0(jobevol$max.season, "-06-30")))

#who has done which job?
jobevol$jobdone <- 1
jobdone <- unique(jobevol[,c("link", "funccat", "jobdone")])
jobdone <- spread(jobdone, "funccat", "jobdone")
jobdone[is.na(jobdone)] <- 0

#add no job category right after player career
nojob <- subset.data.frame(playerlevel)[,c("playername", "link", "max.season")]
nojob$max.season <- as.numeric(as.Date(paste0(nojob$max.season, "-06-30")))
nojob$funccat <- "No Job"
nojob$startdatenum <- nojob$max.season 
jobevol <- rbind.fill(jobevol, nojob)
jobevol <- jobevol[order(jobevol$playername, jobevol$startdatenum),]

#remove if no job at start of career if started job right at end of playing career
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(minstartdatenum = min(startdatenum)) %>%
  ungroup
jobevol <- subset.data.frame(jobevol, !(startdatenum > minstartdatenum & funccat == "No Job"))

#add jobsdone
jobevol <- custom_join(jobevol, jobdone)
jobs <- colnames(jobdone)[2:length(colnames(jobdone))]
jobevol[jobs] <-
  lapply(jobevol[jobs], function(x) ifelse(is.na(x), 0, x))
jobevol$totaljobs <- rowSums(jobevol[jobs])

#set enddate to 2024-01-01 if had no jobs ever
jobevol$enddatenum <- ifelse(jobevol$totaljobs == 0, as.numeric(as.Date("2024-01-01")), jobevol$enddatenum)

#add no job at end if applicable
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(maxend = max(enddatenum, na.rm =T))
jobevol <- as.data.frame(jobevol)
nojob <- unique(subset.data.frame(jobevol, totaljobs >0)[,c("playername", "link", "max.season", "maxend")])
nojob$funccat <- "No Job"
nojob$startdatenum <- nojob$maxend + 1
nojob$enddatenum <- as.numeric(as.Date("2024-01-01"))
jobevol <- rbind.fill(jobevol, nojob)
jobevol <- jobevol[order(jobevol$playername, jobevol$startdatenum),]

#remove players without playing history
jobevol <- subset.data.frame(jobevol, !is.na(max.season))

#recode job categories
jobevol$funccat <- ifelse(!(jobevol$country %in% c("England")) & jobevol$funccat != "No Job", "Job outside UK", jobevol$funccat)

#previous jobs
jobevol <- jobevol %>% 
  group_by(link) %>% 
  mutate(previous_endspell = coalesce(lag(enddatenum), max.season),
         previous_funccat = coalesce(lag(funccat), "No Job"),
         previous_club = coalesce(lag(club), "No Job"),
         previous_clublink = coalesce(lag(clublink), "No Job"),
         previous_country = coalesce(lag(country), "No Job")) %>% 
  ungroup()

#make adaptations if start coach career before end of player career
jobevol$previous_endspell <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat)
jobevol$previous_club <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_country)

jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(count = n())

jobevol$previous_endspell <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat)
jobevol$previous_club <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_country)

#next jobs
jobevol <- jobevol %>% 
  group_by(link) %>% 
  mutate(next_startspell = coalesce(lead(startdatenum), enddatenum),
         next_funccat = coalesce(lead(funccat), "No Job"),
         next_club = coalesce(lead(club), "No Job"),
         next_clublink = coalesce(lead(clublink), "No Job"),
         next_country = coalesce(lead(country), "No Job")) %>% 
  ungroup()

jobevol$enddatenum <- ifelse(is.na(jobevol$enddatenum), jobevol$next_startspell - 1, jobevol$enddatenum)
jobevol$next_startspell <- ifelse(jobevol$enddatenum== 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_startspell)
jobevol$next_funccat <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat)
jobevol$next_club <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_club)
jobevol$next_clublink <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_clublink)
jobevol$next_country <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_country)

#remove last job = no job if they still have a job
jobevol <- subset.data.frame(jobevol, !(startdatenum == 19724 & enddatenum == 19723))

#remove overlap between jobs
jobevol$enddatenum <- ifelse(!is.na(jobevol$next_startspell), ifelse(jobevol$next_startspell <= jobevol$enddatenum, jobevol$next_startspell -1, jobevol$enddatenum), jobevol$enddatenum)

#remove jobs where startdate is after enddate
jobevol <- subset.data.frame(jobevol, startdatenum < enddatenum)

#redo previous jobs
jobevol <- jobevol %>% 
  group_by(link) %>% 
  mutate(previous_endspell = coalesce(lag(enddatenum), max.season),
         previous_funccat = coalesce(lag(funccat), "No Job"),
         previous_club = coalesce(lag(club), "No Job"),
         previous_clublink = coalesce(lag(clublink), "No Job"),
         previous_country = coalesce(lag(country), "No Job")) %>% 
  ungroup()

#make adaptations if start coach career before end of player career
jobevol$previous_endspell <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat)
jobevol$previous_club <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_country)

jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(count = n())

jobevol$previous_endspell <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat)
jobevol$previous_club <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_country)

#redo next jobs
jobevol <- jobevol %>% 
  group_by(link) %>% 
  mutate(next_startspell = coalesce(lead(startdatenum), enddatenum),
         next_funccat = coalesce(lead(funccat), "No Job"),
         next_club = coalesce(lead(club), "No Job"),
         next_clublink = coalesce(lead(clublink), "No Job"),
         next_country = coalesce(lead(country), "No Job")) %>% 
  ungroup()

jobevol$enddatenum <- ifelse(is.na(jobevol$enddatenum), jobevol$next_startspell - 1, jobevol$enddatenum)
jobevol$next_startspell <- ifelse(jobevol$enddatenum== 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_startspell)
jobevol$next_funccat <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat)
jobevol$next_club <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_club)
jobevol$next_clublink <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_clublink)
jobevol$next_country <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_country)

#remove vars
jobs <- subset(unique(jobevol$funccat), unique(jobevol$funccat) != "No Job" & unique(jobevol$funccat) != "Job outside UK")
jobevol <- dplyr::select(jobevol, -c(count, jobs, totaljobs))

#add jobsdone back again
jobevol <- custom_join(jobevol, jobdone)
jobs <- colnames(jobdone)[2:length(colnames(jobdone))]
jobevol[jobs] <-
  lapply(jobevol[jobs], function(x) ifelse(is.na(x), 0, x))
jobevol$totaljobs <- rowSums(jobevol[jobs])

#add a sequence identifier column within each coach's link
jobevol <- jobevol %>% 
  group_by(link) %>% 
  mutate(spellid = cumsum(startdatenum > lag(startdatenum, default = first(startdatenum))) + 1) %>% 
  ungroup()

#data subsets
#must have played between 1990 and 2010 in tier 1 or tier 2 football (excluding cup games)
length(unique(jobevol$link))
links <- unique(subset(combinedplayercareer$link, combinedplayercareer$season >= 1990 & combinedplayercareer$season <= 2010))
jobevol <- subset.data.frame(jobevol, link %in% links)
length(unique(jobevol$link))
#must have finished career before 1 july 2013
jobevol <- subset.data.frame(jobevol, max.season <= as.numeric(as.Date("2013-06-30"))) ###Important
length(unique(jobevol$link))
#subset on those that have race well defined (coded & certain)
# links <- unique(subset(playercareerlevel$link, playercareerlevel$black %in% c(1,0)))
# jobevol <- subset.data.frame(jobevol, link %in% links)
# length(unique(jobevol$link))

#add tiers
tiers <- clubdata[,c("season", "clublink", "tier_alt")]
#keep only lowest (in numbers) tier
tiers <- tiers %>% 
  group_by(clublink, season) %>% 
  filter(tier_alt == min(tier_alt))
tiers <- unique(tiers)
tiers <- as.data.frame(tiers)

#complete missing ones using average tier in data
tieravg <- collap(tiers, tier_alt ~ clublink, FUN = mean)
tieravg$tierround <- round(tieravg$tier_alt)

#tier at end current club
colnames(tiers)[c(1,3)] <- c("endseason", "tier_end")
jobevol$endseason <- ifelse(as.numeric(substr(jobevol$enddate,6,7)) < 7, as.numeric(substr(jobevol$enddate, 1,4)) -1, as.numeric(substr(jobevol$enddate, 1,4)))
jobevol <- custom_join(jobevol, tiers)
jobevol$tier_end <- ifelse(jobevol$funccat == "No Job", "No Job", jobevol$tier_end)
jobevol$tier_end <- ifelse(!(jobevol$country %in% c("England")) & jobevol$funccat != "No Job", "Job outside UK", jobevol$tier_end)
#merge in averaged tiers as alternative
colnames(tieravg)[c(2,3)] <- c("tier_endavg", "tier_endround")
jobevol <- custom_join(jobevol, tieravg)
jobevol$tier_end <- ifelse(is.na(jobevol$tier_end), jobevol$tier_endround, jobevol$tier_end)

#tier start next club
colnames(tiers)[c(1,2,3)] <- c("next_startseason", "next_clublink", "tier_next")
jobevol$next_startspelldate <- as.Date(jobevol$next_startspell, origin='1970-01-01')
jobevol$next_startseason <- ifelse(as.numeric(substr(jobevol$next_startspelldate,6,7)) < 7, as.numeric(substr(jobevol$next_startspelldate, 1,4)) -1, as.numeric(substr(jobevol$next_startspelldate, 1,4)))
jobevol <- custom_join(jobevol, tiers)
jobevol$tier_next <- ifelse(jobevol$next_funccat == "No Job", "No Job", jobevol$tier_next)
jobevol$tier_next <- ifelse(!(jobevol$next_country %in% c("England")) & jobevol$next_funccat != "No Job", "Job outside UK", jobevol$tier_next)
#merge in averaged tiers as alternative
colnames(tieravg) <- c("next_clublink", "tier_nextavg", "tier_nextround")
jobevol <- custom_join(jobevol, tieravg)
jobevol$tier_next <- ifelse(!is.na(jobevol$next_startspelldate) & is.na(jobevol$tier_next), jobevol$tier_nextround, jobevol$tier_next)

#tiers worked per played
jobevol$counter <- 1
tiercount <- collap(jobevol, counter ~ link + tier_next, FUN = sum)
tiercount <- subset.data.frame(tiercount, !is.na(tier_next) & tier_next != "No Job")
tiercount <- spread(tiercount, tier_next, -c(1))
tiercount[is.na(tiercount)] <- 0
tiercount[,c(2:8)][tiercount[,c(2:8)] > 1] <- 1
colnames(tiercount)[2:7] <- paste0("tier", colnames(tiercount)[2:7])

#add some more personal characteristics
jobevol <- custom_join(jobevol, unique(playerlevel[,c("link", "dob")]))
jobevol$yob <- substr(jobevol$dob, nchar(jobevol$dob)-3, nchar(jobevol$dob))
jobevol <- custom_join(jobevol, unique(playerlevel)[,c("link", "natcat")])

#remove errors
jobevol <- subset.data.frame(jobevol, !(link == "https://www.transfermarkt.com/marc-bridge-wilkinson/profil/spieler/39231"  & club == "Liverpool FC YL"))

#make days elapsed since end of playing career
jobevol$days_elapsed <- jobevol$startdatenum - jobevol$max.season

#make state number
jobevol$state <- NA
jobevol$state <- ifelse(jobevol$funccat == "No Job", 1, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Other staff", 2, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Job outside UK", 3, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Youth team management/development", 4, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Assistant manager", 5, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Manager", 6, jobevol$state)


# DATA PREP 7: compile previous and next job data (with funccat_alt and funccat_bn) ----------

#create data about previous and next jobs
jobevol <- unique(coachdata[,c("playername", "link", "club", "clublink", "country", "funccat", "funccat_alt", "funccat_bn", "startdate", "startdatenum", "enddate", "enddatenum", "max.season")])
jobevol$max.season <- as.numeric(as.Date(paste0(jobevol$max.season, "-06-30")))

#who has done which job?
jobevol$jobdone <- 1
jobdone <- unique(jobevol[,c("link", "funccat", "jobdone")])
jobdone <- spread(jobdone, "funccat", "jobdone")
jobdone[is.na(jobdone)] <- 0

#add no job category right after player career
nojob <- subset.data.frame(playerlevel)[,c("playername", "link", "max.season")]
nojob$max.season <- as.numeric(as.Date(paste0(nojob$max.season, "-06-30")))
nojob$funccat <- "No Job"
nojob$funccat_alt <- "No Job"
nojob$funccat_bn <- "No Job"
nojob$startdatenum <- nojob$max.season
jobevol <- rbind.fill(jobevol, nojob)
jobevol <- jobevol[order(jobevol$playername, jobevol$startdatenum),]

#remove if no job at start of career if started job right at end of playing career
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(minstartdatenum = min(startdatenum)) %>%
  ungroup
jobevol <- subset.data.frame(jobevol, !(startdatenum > minstartdatenum & funccat == "No Job"))

#add jobsdone
jobevol <- custom_join(jobevol, jobdone)
jobs <- colnames(jobdone)[2:length(colnames(jobdone))]
jobevol[jobs] <-
  lapply(jobevol[jobs], function(x) ifelse(is.na(x), 0, x))
jobevol$totaljobs <- rowSums(jobevol[jobs])

#set enddate to 2024-01-01 if had no jobs ever
jobevol$enddatenum <- ifelse(jobevol$totaljobs == 0, as.numeric(as.Date("2024-01-01")), jobevol$enddatenum)

#add no job at end if applicable
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(maxend = max(enddatenum, na.rm =T))
jobevol <- as.data.frame(jobevol)
nojob <- unique(subset.data.frame(jobevol, totaljobs >0)[,c("playername", "link", "max.season", "maxend")])
nojob$funccat <- "No Job"
nojob$funccat_alt <- "No Job"
nojob$funccat_bn <- "No Job"
nojob$startdatenum <- nojob$maxend + 1
nojob$enddatenum <- as.numeric(as.Date("2024-01-01"))
jobevol <- rbind.fill(jobevol, nojob)
jobevol <- jobevol[order(jobevol$playername, jobevol$startdatenum),]

#remove players without playing history
jobevol <- subset.data.frame(jobevol, !is.na(max.season))

#recode job categories
jobevol$funccat <- ifelse(!(jobevol$country %in% c("England")) & jobevol$funccat != "No Job", "Job outside UK", jobevol$funccat)

#previous jobs
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(previous_endspell = coalesce(lag(enddatenum), max.season),
         previous_funccat = coalesce(lag(funccat), "No Job"),
         previous_funccat_alt = coalesce(lag(funccat_alt), "No Job"),
         previous_funccat_bn = coalesce(lag(funccat_bn), "No Job"),
         previous_club = coalesce(lag(club), "No Job"),
         previous_clublink = coalesce(lag(clublink), "No Job"),
         previous_country = coalesce(lag(country), "No Job")) %>%
  ungroup()

#make adaptations if start coach career before end of player career
jobevol$previous_endspell <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat)
jobevol$previous_funccat_alt <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat_alt)
jobevol$previous_funccat_bn <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat_bn)
jobevol$previous_club <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_country)

jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(count = n())

jobevol$previous_endspell <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat)
jobevol$previous_funccat_alt <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat_alt)
jobevol$previous_funccat_bn <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat_bn)
jobevol$previous_club <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_country)

#next jobs
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(next_startspell = coalesce(lead(startdatenum), enddatenum),
         next_funccat = coalesce(lead(funccat), "No Job"),
         next_funccat_alt = coalesce(lead(funccat_alt), "No Job"),
         next_funccat_bn = coalesce(lead(funccat_bn), "No Job"),
         next_club = coalesce(lead(club), "No Job"),
         next_clublink = coalesce(lead(clublink), "No Job"),
         next_country = coalesce(lead(country), "No Job")) %>%
  ungroup()

jobevol$enddatenum <- ifelse(is.na(jobevol$enddatenum), jobevol$next_startspell - 1, jobevol$enddatenum)
jobevol$next_startspell <- ifelse(jobevol$enddatenum== 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_startspell)
jobevol$next_funccat <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat)
jobevol$next_funccat_alt <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat_alt)
jobevol$next_funccat_bn <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat_bn)
jobevol$next_club <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_club)
jobevol$next_clublink <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_clublink)
jobevol$next_country <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_country)

#remove last job = no job if they still have a job
jobevol <- subset.data.frame(jobevol, !(startdatenum == 19724 & enddatenum == 19723))

#remove overlap between jobs
jobevol$enddatenum <- ifelse(!is.na(jobevol$next_startspell), ifelse(jobevol$next_startspell <= jobevol$enddatenum, jobevol$next_startspell -1, jobevol$enddatenum), jobevol$enddatenum)

#remove jobs where startdate is after enddate
jobevol <- subset.data.frame(jobevol, startdatenum < enddatenum)

#redo previous jobs
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(previous_endspell = coalesce(lag(enddatenum), max.season),
         previous_funccat = coalesce(lag(funccat), "No Job"),
         previous_funccat_alt = coalesce(lag(funccat_alt), "No Job"),
         previous_funccat_bn = coalesce(lag(funccat_bn), "No Job"),
         previous_club = coalesce(lag(club), "No Job"),
         previous_clublink = coalesce(lag(clublink), "No Job"),
         previous_country = coalesce(lag(country), "No Job")) %>%
  ungroup()

#make adaptations if start coach career before end of player career
jobevol$previous_endspell <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat)
jobevol$previous_funccat_alt <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat_alt)
jobevol$previous_funccat_bn <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_funccat_bn)
jobevol$previous_club <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$startdatenum <= jobevol$max.season, NA, jobevol$previous_country)

jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(count = n())

jobevol$previous_endspell <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_endspell)
jobevol$previous_funccat <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat)
jobevol$previous_funccat_alt <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat_alt)
jobevol$previous_funccat_bn <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_funccat_bn)
jobevol$previous_club <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_club)
jobevol$previous_clublink <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_clublink)
jobevol$previous_country <- ifelse(jobevol$count == 1 | (jobevol$funccat == "No Job" & jobevol$previous_funccat == "No Job"), NA, jobevol$previous_country)

#redo next jobs
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(next_startspell = coalesce(lead(startdatenum), enddatenum),
         next_funccat = coalesce(lead(funccat), "No Job"),
         next_funccat_alt = coalesce(lead(funccat_alt), "No Job"),
         next_funccat_bn = coalesce(lead(funccat_bn), "No Job"),
         next_club = coalesce(lead(club), "No Job"),
         next_clublink = coalesce(lead(clublink), "No Job"),
         next_country = coalesce(lead(country), "No Job")) %>%
  ungroup()

jobevol$enddatenum <- ifelse(is.na(jobevol$enddatenum), jobevol$next_startspell - 1, jobevol$enddatenum)
jobevol$next_startspell <- ifelse(jobevol$enddatenum== 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_startspell)
jobevol$next_funccat <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat)
jobevol$next_funccat_alt <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat_alt)
jobevol$next_funccat_bn <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_funccat_bn)
jobevol$next_club <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_club)
jobevol$next_clublink <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_clublink)
jobevol$next_country <- ifelse(jobevol$enddatenum == 19723 | is.na(jobevol$enddatenum), NA, jobevol$next_country)

#remove vars
jobs <- subset(unique(jobevol$funccat), unique(jobevol$funccat) != "No Job" & unique(jobevol$funccat) != "Job outside UK")
jobevol <- dplyr::select(jobevol, -c(count, jobs, totaljobs))

#add jobsdone back again
jobevol <- custom_join(jobevol, jobdone)
jobs <- colnames(jobdone)[2:length(colnames(jobdone))]
jobevol[jobs] <-
  lapply(jobevol[jobs], function(x) ifelse(is.na(x), 0, x))
jobevol$totaljobs <- rowSums(jobevol[jobs])

#add a sequence identifier column within each coach's link
jobevol <- jobevol %>%
  group_by(link) %>%
  mutate(spellid = cumsum(startdatenum > lag(startdatenum, default = first(startdatenum))) + 1) %>%
  ungroup()

#data subsets
#must have played between 1990 and 2010 in tier 1 or tier 2 football (excluding cup games)
length(unique(jobevol$link))
links <- unique(subset(combinedplayercareer$link, combinedplayercareer$season >= 1990 & combinedplayercareer$season <= 2010))
jobevol <- subset.data.frame(jobevol, link %in% links)
length(unique(jobevol$link))
#must have finished career before 1 july 2013
jobevol <- subset.data.frame(jobevol, max.season <= as.numeric(as.Date("2013-06-30"))) ###Important
length(unique(jobevol$link))
#subset on those that have race well defined (coded & certain)
# links <- unique(subset(playercareerlevel$link, playercareerlevel$black %in% c(1,0)))
# jobevol <- subset.data.frame(jobevol, link %in% links)
# length(unique(jobevol$link))

#add tiers
tiers <- clubdata[,c("season", "clublink", "tier_alt")]
#keep only lowest (in numbers) tier
tiers <- tiers %>%
  group_by(clublink, season) %>%
  filter(tier_alt == min(tier_alt))
tiers <- unique(tiers)
tiers <- as.data.frame(tiers)

#complete missing ones using average tier in data
tieravg <- collap(tiers, tier_alt ~ clublink, FUN = mean)
tieravg$tierround <- round(tieravg$tier_alt)

#tier at end current club
colnames(tiers)[c(1,3)] <- c("endseason", "tier_end")
jobevol$endseason <- ifelse(as.numeric(substr(jobevol$enddate,6,7)) < 7, as.numeric(substr(jobevol$enddate, 1,4)) -1, as.numeric(substr(jobevol$enddate, 1,4)))
jobevol <- custom_join(jobevol, tiers)
jobevol$tier_end <- ifelse(jobevol$funccat == "No Job", "No Job", jobevol$tier_end)
jobevol$tier_end <- ifelse(!(jobevol$country %in% c("England")) & jobevol$funccat != "No Job", "Job outside UK", jobevol$tier_end)
#merge in averaged tiers as alternative
colnames(tieravg)[c(2,3)] <- c("tier_endavg", "tier_endround")
jobevol <- custom_join(jobevol, tieravg)
jobevol$tier_end <- ifelse(is.na(jobevol$tier_end), jobevol$tier_endround, jobevol$tier_end)

#tier start next club
colnames(tiers)[c(1,2,3)] <- c("next_startseason", "next_clublink", "tier_next")
jobevol$next_startspelldate <- as.Date(jobevol$next_startspell, origin='1970-01-01')
jobevol$next_startseason <- ifelse(as.numeric(substr(jobevol$next_startspelldate,6,7)) < 7, as.numeric(substr(jobevol$next_startspelldate, 1,4)) -1, as.numeric(substr(jobevol$next_startspelldate, 1,4)))
jobevol <- custom_join(jobevol, tiers)
jobevol$tier_next <- ifelse(jobevol$next_funccat == "No Job", "No Job", jobevol$tier_next)
jobevol$tier_next <- ifelse(!(jobevol$next_country %in% c("England")) & jobevol$next_funccat != "No Job", "Job outside UK", jobevol$tier_next)
#merge in averaged tiers as alternative
colnames(tieravg) <- c("next_clublink", "tier_nextavg", "tier_nextround")
jobevol <- custom_join(jobevol, tieravg)
jobevol$tier_next <- ifelse(!is.na(jobevol$next_startspelldate) & is.na(jobevol$tier_next), jobevol$tier_nextround, jobevol$tier_next)

#tiers worked per played
jobevol$counter <- 1
tiercount <- collap(jobevol, counter ~ link + tier_next, FUN = sum)
tiercount <- subset.data.frame(tiercount, !is.na(tier_next) & tier_next != "No Job")
tiercount <- spread(tiercount, tier_next, -c(1))
tiercount[is.na(tiercount)] <- 0
tiercount[,c(2:8)][tiercount[,c(2:8)] > 1] <- 1
colnames(tiercount)[2:7] <- paste0("tier", colnames(tiercount)[2:7])

#add some more personal characteristics
jobevol <- custom_join(jobevol, unique(playerlevel[,c("link", "dob")]))
jobevol$yob <- substr(jobevol$dob, nchar(jobevol$dob)-3, nchar(jobevol$dob))
jobevol <- custom_join(jobevol, unique(playerlevel)[,c("link", "natcat")])

#remove errors
jobevol <- subset.data.frame(jobevol, !(link == "https://www.transfermarkt.com/marc-bridge-wilkinson/profil/spieler/39231"  & club == "Liverpool FC YL"))

#make days elapsed since end of playing career
jobevol$days_elapsed <- jobevol$startdatenum - jobevol$max.season

#make state number (funccat)
jobevol$state <- NA
jobevol$state <- ifelse(jobevol$funccat == "No Job", 1, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Other staff", 2, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Job outside UK", 3, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Youth team management/development", 4, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Assistant manager", 5, jobevol$state)
jobevol$state <- ifelse(jobevol$funccat == "Manager", 6, jobevol$state)

#make state number for funccat_alt (no "Job outside UK" category)
jobevol$state_alt <- NA
jobevol$state_alt <- ifelse(jobevol$funccat_alt == "No Job", 1, jobevol$state_alt)
jobevol$state_alt <- ifelse(jobevol$funccat_alt == "Other staff", 2, jobevol$state_alt)
jobevol$state_alt <- ifelse(jobevol$funccat_alt == "Youth team management/development", 3, jobevol$state_alt)
jobevol$state_alt <- ifelse(jobevol$funccat_alt == "Assistant manager", 4, jobevol$state_alt)
jobevol$state_alt <- ifelse(jobevol$funccat_alt == "Manager", 5, jobevol$state_alt)

#make state number for funccat_bn (England/ROW split)
jobevol$state_bn <- NA
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "No Job", 1, jobevol$state_bn)
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "Other staff", 2, jobevol$state_bn)
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "Youth England", 3, jobevol$state_bn)
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "Youth ROW", 4, jobevol$state_bn)
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "Assistant manager England", 5, jobevol$state_bn)
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "Assistant manager ROW", 6, jobevol$state_bn)
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "Manager England", 7, jobevol$state_bn)
jobevol$state_bn <- ifelse(jobevol$funccat_bn == "Manager ROW", 8, jobevol$state_bn)


# DATA PREP 8: duration analysis preparation --------------------------------------------

#all job spells in coachdata (England and rest of world), matched to games from allgames.Rdata.
#games come from the TM scrape rather than the England-only odds file, so no betting odds are
#available and cumulative surprise is not constructed.

load(file.path(inputwd, "allgames.Rdata"))

#collapse the player-appearance level scrape to one row per team per match
data.table::setDT(gameframes)
gameteam <- unique(gameframes[, .(league, season, id, date, homescore, awayscore,
                                  hometeamlink, awayteamlink, home)])
rm(gameframes)

#TM writes dates as "Sun, 8/24/03", "Sat, 23/08/03", "Sat, 31/07/1999" or "Sat, Feb 10, 1973"
datestr <- sub("^[A-Za-z]{3}, ", "", trimws(gameteam$date))
gdate <- as.Date(rep(NA_character_, length(datestr)))
sel <- grepl("^[A-Za-z]{3} [0-9]{1,2}, [0-9]{4}$", datestr)
gdate[sel] <- as.Date(datestr[sel], format = "%b %d, %Y")
sel <- grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$", datestr)
gdate[sel] <- as.Date(datestr[sel], format = "%m/%d/%Y")
sel <- sel & is.na(gdate)
gdate[sel] <- as.Date(datestr[sel], format = "%d/%m/%Y")
sel <- grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2}$", datestr)
gdate[sel] <- as.Date(datestr[sel], format = "%m/%d/%y")
sel <- sel & is.na(gdate)
gdate[sel] <- as.Date(datestr[sel], format = "%d/%m/%y")
gameteam[, gdate := gdate]
rm(datestr, gdate, sel)

#clubs are matched on the numeric TM club id, so no club-name key is needed
gameteam[, teamlink := data.table::fifelse(home == 1, hometeamlink, awayteamlink)]
gameteam[, `:=`(clubid  = str_extract(teamlink, "(?<=verein/)[0-9]+"),
                datenum = as.numeric(gdate),
                sc      = suppressWarnings(as.numeric(data.table::fifelse(home == 1, homescore, awayscore))),
                op_sc   = suppressWarnings(as.numeric(data.table::fifelse(home == 1, awayscore, homescore))))]

#matchid is stored as the bare numeric Transfermarkt match id rather than the full report
#URL, so it is directly comparable with the id the fixture scrape extracts and the two
#sources can be deduplicated against each other
gamedata_all <- unique(gameteam[!is.na(clubid) & !is.na(datenum) & !is.na(sc) & !is.na(op_sc),
                                .(clubid, datenum, date = gdate, div = league, season,
                                  matchid = str_extract(id, "[0-9]+$"), sc, op_sc)])
rm(gameteam)

#add the targeted fixture scrape (Scraper files/spellgamescrape.R), which covers club-seasons
#the league scrape never reached: MLS, Wales, Australia, Norway, Ireland, Sweden, Denmark, the
#lower Scottish tiers, the academy leagues, and seasons in which a club otherwise present in
#allgames played outside the covered divisions. div holds competition names for these rows
#rather than league codes.
#
#The scrape fetches every season of a qualifying spell, including seasons allgames already
#covers, so the two sources overlap and the same game can arrive twice with a different div
#label (league code vs competition name). Deduplication is therefore on club and Transfermarkt
#match id rather than on the whole row, which would leave both copies in and double-count
#games. allgames is kept where both have a game, since its div codes are the cleaner ones.
#
#set use_scraped_games to FALSE before running this section to reproduce the earlier
#allgames-only analysis; it defaults to TRUE so the normal pipeline is unaffected
if (!exists("use_scraped_games")) use_scraped_games <- TRUE
scrapedfile <- file.path(inputwd, "gamedata_scraped.Rdata")
if (use_scraped_games && file.exists(scrapedfile)) {
  load(scrapedfile)   # gamedata_scraped
  gamedata_scraped <- data.table::as.data.table(gamedata_scraped)
  gamedata_scraped[, `:=`(clubid = as.character(clubid), matchid = as.character(matchid),
                          div = as.character(div), season = as.numeric(season))]
  gamedata_all[, matchid := as.character(matchid)]
  gamedata_all[, src := 1L]   # allgames first, so it wins ties on dedup
  gamedata_scraped[, src := 2L]
  gamedata_all <- rbind(
    gamedata_all,
    gamedata_scraped[, .(clubid, datenum, date, div, season, matchid, sc, op_sc, src)])
  #games with no usable match id fall back to club and date as the identity key
  gamedata_all[, dedupkey := data.table::fifelse(
    is.na(matchid) | matchid == "",
    paste0("d_", clubid, "_", datenum),
    paste0("m_", clubid, "_", matchid))]
  data.table::setorder(gamedata_all, dedupkey, src)
  n_before <- nrow(gamedata_all)
  gamedata_all <- gamedata_all[!duplicated(dedupkey)]
  cat("dropped", n_before - nrow(gamedata_all), "duplicate games across the two sources\n")
  gamedata_all[, c("src", "dedupkey") := NULL]
  rm(gamedata_scraped)
  cat("game pool after adding scrape:", nrow(gamedata_all), "rows,",
      length(unique(gamedata_all$clubid)), "clubs\n")
} else {
  cat("running on allgames only:", nrow(gamedata_all), "rows,",
      length(unique(gamedata_all$clubid)), "clubs\n")
}

#prepare the job spells: every spell, not just England, and keyed on a unique spell id because
#the same coach can hold two roles at one club from the same start date
coachdata_alt <- coachdata
coachdata_alt$jobspellid <- seq_len(nrow(coachdata_alt))
coachdata_alt$clubid <- str_extract(coachdata_alt$clublink_alt, "(?<=verein/)[0-9]+")
spellframe <- data.table::as.data.table(
  subset.data.frame(coachdata_alt, !is.na(clubid) & !is.na(startdatenum) & !is.na(enddatenum)))

#make data frame containing all games for each job spell
masterdf <- gamedata_all[spellframe,
                         on = .(clubid, datenum >= startdatenum, datenum <= enddatenum),
                         allow.cartesian = TRUE, nomatch = 0L,
                         .(jobspellid, link, playername, dob, club, clubid, clublink, country,
                           funccat, funccat_alt, funccat_bn, youth, coachlicence, licencedum,
                           startdatenum, enddatenum, matchid, date, datenum = x.datenum,
                           div, season, sc, op_sc)]
save(masterdf, file = file.path(datawd, "durationdata.Rdata"))

load(file.path(datawd, "durationdata.Rdata"))
masterdfreg <- data.table::as.data.table(masterdf)

#make points variable
masterdfreg[, points := data.table::fifelse(sc < op_sc, 0, data.table::fifelse(sc == op_sc, 1, 3))]

#make running points performance variable
data.table::setorder(masterdfreg, jobspellid, datenum, matchid)
masterdfreg[, cumulpoints := cumsum(points), by = jobspellid]

#cumulative games in current job
masterdfreg[, cumulgames := seq_len(.N), by = jobspellid]
masterdfreg[, cumulpoints_last10 := data.table::fifelse(
  cumulgames == 10, cumulpoints, cumulpoints - data.table::shift(cumulpoints, 10L)),
  by = jobspellid]

#cumulative games since start of coaching career
data.table::setorder(masterdfreg, link, datenum, matchid)
masterdfreg[, cumulgamestot := seq_len(.N), by = link]

#cumulative games since start of coaching career at start of job
masterdfreg[, cumulgamesstart := min(cumulgamestot) - 1, by = jobspellid]

#make indicator for leaving job
masterdfreg[, leftjob := data.table::fifelse(cumulgames == max(cumulgames) & enddatenum != 19723, 1, 0),
            by = jobspellid]

regframe <- as.data.frame(masterdfreg)

# Order the funccat categories
funccat_order <- rev(c("Other staff", "Youth team management/development", "Assistant manager", "Manager"))

# Convert funccat to a factor with the specified order. funccat_alt is used because funccat
# collapses every non-England spell into "Job outside UK", which the England/ROW split needs kept apart.
regframe$funccat <- factor(regframe$funccat_alt, levels = funccat_order)
regframe <- within(regframe, funccat <- relevel(funccat, ref = 4))
table(regframe$funccat)

#flag whether the spell itself is at an English club
regframe$engspell <- ifelse(!is.na(regframe$country) & regframe$country == "England", 1, 0)
table(regframe$engspell, regframe$funccat)

#No divisional tier variable is built. It could be read off the trailing digit of the league
#codes in allgames, but the fixture scrape brings div in as competition names: 423 distinct
#values, only 48 ending in a digit, with sponsor renames splitting single leagues (Scotland's
#top flight appears as Premiership, Betway Premiership and Scottish Premiership). A reliable
#tier would have to come from the tier column in clubdata, merged on club and season.

#merge in player career details
regframe <- plyr::join(regframe, dplyr::select(playerlevel, -c(min.season, max.season)))

regframe$cumulgames <- regframe$cumulgames - 1
regframe$cumulgames_next <- regframe$cumulgames + 1

#nationality dummies, built exactly as elsewhere in the paper: the three categories are
#exhaustive, so regressions include ukplus and foreign and leave UK-only as the reference
regframe$dualnational <- ifelse(regframe$natcat2 == "missing", 0, 1)
regframe$ukonly  <- ifelse(regframe$natcat == "UK" & regframe$dualnational == 0, 1, 0)
regframe$ukplus  <- ifelse((regframe$natcat == "UK"  & regframe$dualnational == 1) |
                             (regframe$natcat2 == "UK" & regframe$dualnational == 1), 1, 0)
regframe$foreign <- ifelse(!(regframe$natcat == "UK" | regframe$natcat2 == "UK"), 1, 0)

#drop national team spells. These are not club jobs, and youth national sides do carry
#fixtures on Transfermarkt (U20 Elite League, OFC championships, the Algarve Tournament),
#so without this they leak into the sample rather than simply failing to match.
regframe <- subset.data.frame(regframe, is.na(country) | country != "international")

#add race
regframe <- plyr::join(regframe, responsesperid)
regframe$black <- ifelse(regframe$blackmaj == 1, 1, 0)

#filter out players with unclear race
regframe <- subset.data.frame(regframe, nomaj == 0)

#only spells longer than 10 games. Set require_ten_games to FALSE before running this
#section to keep shorter spells, which is possible for specifications with no performance
#control; it defaults to TRUE so the normal pipeline is unaffected.
if (!exists("require_ten_games")) require_ten_games <- TRUE
if (require_ten_games) regframe <- subset.data.frame(regframe, !is.na(cumulpoints_last10))

#add next job data
regframe <- plyr::join(regframe, jobevol[c("link", "startdatenum", "enddatenum", "next_funccat")])

# SAVE ENVIRONMENT DATA PREP 1-8 --------------------------------------------------------


#doubtful if should keep: jobdate, jobdate2
keepdata <- c("counter_simact", "countlist_simeao", "countlist_simpao", "countlist_simom", "countlist_simoma",
              "countlist_simact", "countlist_simact_all", "coachdata", "date_sequence_sim", "jobevol", "natframe", 
              "playercareerlevel", "playerlevel", "plotframe", "plotframe_share", "plotframe_share_simact", "plotframe_share_simeao", "plotframe_share_simpao",
              "plotframe_share_simom", "plotframe_share_simoma", "regframe", "responsesperid", "set1_days", "statelist_simact", "tiercount")

keepvecs <- c("axislabsize", "cbp5", "cbp6", "cbp8", "custom_join", "cutthres", "facetsize", "geomtextsize",
              "inputwd", "datawd", "pkgwd", "insiders", "joborder", "jobs",
              "jobstates", "labsize", "legendtextsize", "legsize", "links", "outputwdfiles",
              "outputwdresults", "rmlinks", "textsize", "ticksize", "ukcountries", "rhs_vars")

regressions <- c("simreg_multi", "startregmulti")

keeps <- c(keepdata, keepvecs)
removes <- setdiff(ls(), keeps)
rm(list=setdiff(ls(), keeps))

# Identifier columns that the code above needs while it runs but that nothing reads
# afterwards: not the results script, and not the simulation half of this script. They are
# dropped here rather than in the anonymisation step, because the anonymisation step works on
# the source files and these columns are built in between, so stripping them there would
# either break this script or be undone by the next rebuild. Only `link` survives, since it is
# the join key, and `dob`, which supplies year of birth.
strip_ids <- list(
  regframe          = c("matchid", "clubid", "clublink", "playername", "nationality_alt"),
  coachdata         = c("clublink", "clublink_alt", "playername"),
  jobevol           = c("clublink", "previous_clublink", "next_clublink", "spellid",
                        "playername"),
  playerlevel       = c("playername", "nationality_alt"),
  natframe          = c("playername", "nationality_alt"),
  playercareerlevel = c("complink", "leaguelink", "clublink_alt", "clublink", "playername")
)
for (.o in names(strip_ids)) {
  if (!exists(.o)) next
  .x <- get(.o)
  .drop <- intersect(names(.x), strip_ids[[.o]])
  if (length(.drop)) {
    assign(.o, .x[, setdiff(names(.x), .drop), drop = FALSE])
    cat("stripped from", .o, ":", paste(.drop, collapse = ", "), "\n")
  }
}
rm(list = intersect(ls(all.names = TRUE), c(".o", ".x", ".drop")))

setwd(datawd)
save.image(file.path(datawd, "analysis_environment.RData"))

cat("\n>>> analysis_environment.RData written -",
    round(file.size(file.path(datawd, "analysis_environment.RData")) / 1024^2, 1), "MB\n")
cat(">>> next: run HPS Simulation preparation.R once for each uksim setting\n")

# Leave the working directory where it started. Each script sets pkgwd from getwd(), so a
# script that ends inside data/ would make the next one in the sequence resolve every path
# one level down and create stray data/ and output/ folders there.
setwd(pkgwd)
