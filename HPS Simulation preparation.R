# preamble ----------------------------------------------------------------
# HPS Simulation preparation.R
#
# Builds the simulation environments that the simulation section of HPS Racial Gaps.R loads.
# Run "HPS Data preparation.R" first: this script picks up from the analysis environment that
# script writes, data/analysis_environment.RData, and adds nothing to it.
#
# The results code loads three simulation environments, one per sample, so this script has to
# be run three times, once for each setting of uksim below:
#
#   uksim = 0   all former players           -> data/data_environment_all_<reps>.RData
#   uksim = 1   UK nationals only            -> data/data_environment_uk_<reps>.RData
#   uksim = 2   propensity-score matched     -> data/data_environment_ms_<reps>.RData
#
# <reps> is the replication count set below, so a 1-replication test run writes its own files
# and cannot overwrite the 1000-replication ones the package ships with.
#
# Each run fits the multinomial transition model and the 1000-replication simulation, and
# takes several hours. The three environments ship with the package, so the results code can
# be run without repeating them.
#
# These sections were DATA PREP 9-17 of HPS Data preparation.R before being split out; they
# are numbered SIM PREP 1-9 here.

#clear environment
rm(list = ls())

#set seed
set.seed(123)

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

# Set this to the folder that contains this script.
# Everything else is resolved relative to it, so the package runs anywhere.
pkgwd <- getwd()   # or: pkgwd <- "/path/to/Replication Package"

# analysis_environment.RData is a save.image(), so loading it also restores the working
# directories of the machine that built it. Load first, then re-derive every path from the
# current location, so this script runs from wherever the package sits.
load(file.path(pkgwd, "data", "analysis_environment.RData"))

pkgwd           <- getwd()
datawd          <- file.path(pkgwd, "data")
inputwd         <- datawd
outputwdfiles   <- file.path(pkgwd, "output")
outputwdresults <- outputwdfiles
setwd(datawd)

#replications for simulation
reps <- 10

#Which sample to build (see the header). The restriction itself is applied at the end of
#SIM PREP 1, so everything downstream is fitted on the selected sample.
uksim <- 0

# SIM PREP 1: multi stage analysis preparation --------------------------------------------------

#data set up

plotframe <- jobevol
plotframe <- subset.data.frame(plotframe, !(link %in% rmlinks)) #remove people whose job start and enddates have been miscoded
nall <- length(unique(plotframe$link))

#plyr::join in race codings
plotframe <- plyr::join(plotframe, responsesperid)
plotframe$black <- ifelse(plotframe$blackmaj == 1, 1, 0)

#filter out players with unclear race
plotframe <- subset.data.frame(plotframe, nomaj == 0)

#define race counts
races <- table(unique(plotframe[,c("link", "black")])$black)
sum(races)

#must have finished 10 years ago and must have race well defined
plotframe <- subset.data.frame(plotframe, max.season <= as.numeric(as.Date("2013-06-30")) & black %in% c(0,1))

plotframe$next_funccat_bn <- ifelse(plotframe$next_startspell - plotframe$enddatenum > 1, "No Job", plotframe$next_funccat_bn)

# Generate a sequence of dates to represent the missing days
date_sequence_sim <- plotframe %>%
  group_by(link) %>%
  summarize(min_date = min(startdatenum), max_date = as.numeric(as.Date("2024-01-01"))) %>%
  mutate(date = purrr::map2(min_date, max_date, seq)) %>%
  unnest(date)

#merge back in the jobs
jobdate <- unique(plotframe[,c("link", "startdatenum", "funccat_bn", "black")])
colnames(jobdate)[2] <- "date"
jobdate2 <- unique(plotframe[,c("link", "enddatenum", "next_funccat_bn", "black")])
jobdate2$enddatenum <- jobdate2$enddatenum + 1
colnames(jobdate2)[c(2,3)] <- c("date", "funccat_bn")
jobdate <- unique(rbind(jobdate, jobdate2))
date_sequence_sim <- plyr::join(date_sequence_sim, jobdate)
date_sequence_sim <- date_sequence_sim %>%
  fill(funccat_bn, black)

# Calculate days elapsed since the first startdatenum for each player
date_sequence_sim$days_elapsed <- date_sequence_sim$date - date_sequence_sim$min_date

# #subset on first 10 years
date_sequence_sim <- subset.data.frame(date_sequence_sim, days_elapsed <= 3650)

#join in playerstats
date_sequence_sim <- join(date_sequence_sim, playerlevel)

#years elapsed
date_sequence_sim$yearselapsed <- floor(date_sequence_sim$days_elapsed/365)


date_sequence_sim$id <- group(date_sequence_sim$link)

date_sequence_sim$yob <- substr(date_sequence_sim$dob, nchar(date_sequence_sim$dob)-3, nchar(date_sequence_sim$dob))

date_sequence_sim <- subset.data.frame(date_sequence_sim, days_elapsed %in% seq(1,22222,14))

#previous job (previos day)
date_sequence_sim <- date_sequence_sim %>%
  arrange(id, days_elapsed) %>%
  group_by(id) %>%
  mutate(funccat_bn_lag = lag(funccat_bn, n = 1, default = NA))

test <- subset.data.frame(date_sequence_sim, yearselapsed == 4 & black == 0 & funccat_bn_lag == "Manager England")
table(date_sequence_sim$funccat_bn)

date_sequence_sim$uk <- ifelse(date_sequence_sim$nationality %in% ukcountries | date_sequence_sim$nationality2 %in% ukcountries | date_sequence_sim$nationality3 %in% ukcountries, 1 ,0)

#dual national
date_sequence_sim$dualnational <- ifelse(date_sequence_sim$natcat2 == "missing", 0, 1)
table(date_sequence_sim$dualnational)

#ukness
date_sequence_sim$ukonly <- ifelse(date_sequence_sim$natcat == "UK" & date_sequence_sim$dualnational == 0, 1, 0)
date_sequence_sim$ukplus <- ifelse((date_sequence_sim$natcat == "UK" & date_sequence_sim$dualnational == 1) | (date_sequence_sim$natcat2 == "UK" & date_sequence_sim$dualnational == 1), 1, 0)
date_sequence_sim$foreign <- ifelse(!(date_sequence_sim$natcat == "UK" | date_sequence_sim$natcat2 == "UK"), 1, 0)

#sample for the simulation, set by uksim at the top of this script. The simulation section of
#HPS Racial Gaps.R loops over uksim 0, 1 and 2 and loads a separate environment for each, so
#this script has to be run three times, once per setting. Everything from here on -- the
#multinomial model, the simulations and the shares they produce -- is fitted on whichever
#sample is selected here.
#The sample is recorded as a set of links, because SIM PREP 3 rebuilds the career panel from
#jobevol to compute the actual prevalence series. Restricting only date_sequence_sim here would
#leave that series on the full sample, so plotframe_share would disagree with everything else
#in the environment.
sim_links <- NULL

if(uksim == 1){
  #UK nationals only
  sim_links <- unique(date_sequence_sim$link[date_sequence_sim$foreign == 0])
}

if(uksim == 2){
  #propensity-score matched sample. Same probit and same 1:1 match without replacement as
  #DATA PREP 2 of HPS Racial Gaps.R, rebuilt from playerlevel here so that this script does
  #not depend on the results script having been run first. The seed and the row order of
  #psmframe both matter: Match breaks ties at random, so changing either changes the match.
  psmframe <- plyr::join(playerlevel, responsesperid)
  psmframe$black <- ifelse(psmframe$blackmaj == 1, 1, 0)
  psmframe <- subset.data.frame(psmframe, nomaj == 0)
  psmframe <- subset.data.frame(psmframe, link %in% unique(jobevol$link))

  coachedflag <- unique(data.frame(link = coachdata$link, coached = 1))
  psmframe <- custom_join(psmframe, coachedflag)
  psmframe$coached <- ifelse(is.na(psmframe$coached), 0, psmframe$coached)

  psmframe$dualnational <- ifelse(psmframe$natcat2 == "missing", 0, 1)
  psmframe$ukonly <- ifelse(psmframe$natcat == "UK" & psmframe$dualnational == 0, 1, 0)
  psmframe$ukplus <- ifelse((psmframe$natcat == "UK" & psmframe$dualnational == 1) | (psmframe$natcat2 == "UK" & psmframe$dualnational == 1), 1, 0)
  psmframe$foreign <- ifelse(!(psmframe$natcat == "UK" | psmframe$natcat2 == "UK"), 1, 0)
  psmframe$yob <- substr(psmframe$dob, nchar(psmframe$dob)-3, nchar(psmframe$dob))

  psmprobit <- glm(black ~ I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
                     I(sum_minplayed_tier2/60) + sum_goals_tier2 + sum_assists_tier2 +
                     I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
                     attacker + midfielder + defender + ukonly + ukplus + foreign +
                     natteam + natapps + as.factor(yob),
                   family = binomial(link = "probit"), data = psmframe)

  #Seed immediately before the call, not just in the preamble: Match breaks ties at random,
  #and anything upstream that draws from the RNG shifts the state and changes which controls
  #are picked. Roughly a fifth of the matched set moves between seeds. DATA PREP 2 of
  #HPS Racial Gaps.R seeds at the same point with the same value.
  set.seed(123)
  psmatch <- Match(Y = psmframe$coached, Tr = psmframe$black, X = psmprobit$fitted,
                   M = 1, ties = TRUE, replace = FALSE)
  matchedlinks <- psmframe$link[c(psmatch$index.control, psmatch$index.treated)]

  sim_links <- unique(matchedlinks)
  rm(psmframe, coachedflag, psmprobit, psmatch)
}

if (!is.null(sim_links))
  date_sequence_sim <- subset.data.frame(date_sequence_sim, link %in% sim_links)

cat("simulation sample: uksim =", uksim, "-",
    length(unique(date_sequence_sim$link)), "players,",
    nrow(date_sequence_sim), "rows\n")


#regression
# define states and set factor levels for outcome (reference = No Job)
jobstates <- c("Manager England", "Manager ROW", "Assistant manager England",
               "Assistant manager ROW", "Youth England", "Youth ROW", "Other staff", "No Job")
date_sequence_sim$funccat_bn <- factor(date_sequence_sim$funccat_bn, levels = jobstates)

simreg_multi <- multinom(
  funccat_bn ~
    as.factor(black)*days_elapsed*as.factor(funccat_bn_lag) +
    I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
    I(sum_minplayed_tier2/60) + sum_goals_tier2  + sum_assists_tier2 +
    I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
    attacker + midfielder + defender + as.numeric(yob) + ukplus + foreign + natteam + natapps,
  data  = date_sequence_sim,
  decay = 1e-3,
  maxit = 1000
)
# SIM PREP 2: simulate initial state -------------------------------------

#grab starting data
start <- subset.data.frame(date_sequence_sim, days_elapsed == 15)


#regression
start$funccat_bn_lag <- factor(start$funccat_bn_lag, levels = jobstates)
startregmulti <- multinom(
  funccat_bn_lag ~
    as.factor(black) +
    I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
    I(sum_minplayed_tier2/60) + sum_goals_tier2  + sum_assists_tier2 +
    I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
    attacker + midfielder + defender + as.numeric(yob) + ukplus + foreign + natteam + natapps,
  data  = start,
  decay = 1e-3,
  maxit = 1000
)

# SIM PREP 3: Actual prevalence -------------------------------------------------------

plotframe <- jobevol
plotframe <- subset.data.frame(plotframe, !(link %in% rmlinks)) #remove people whose job start and enddates have been miscoded
nall <- length(unique(plotframe$link))

#plyr::join in race codings
plotframe <- plyr::join(plotframe, responsesperid)
plotframe$black <- ifelse(plotframe$blackmaj == 1, 1, 0)

#filter out players with unclear race
plotframe <- subset.data.frame(plotframe, nomaj == 0)

#define race counts
races <- table(unique(plotframe[,c("link", "black")])$black)
sum(races)

#must have finished 10 years ago and must have race well defined
plotframe <- subset.data.frame(plotframe, max.season <= as.numeric(as.Date("2013-06-30")) & black %in% c(0,1))

#same sample restriction as SIM PREP 1, so the actual series is computed on the sample the
#simulations are fitted to rather than on all former players
if (!is.null(sim_links))
  plotframe <- subset.data.frame(plotframe, link %in% sim_links)

plotframe$next_funccat_bn <- ifelse(plotframe$next_startspell - plotframe$enddatenum > 1, "No Job", plotframe$next_funccat_bn)

# Generate a sequence of dates to represent the missing days
date_sequence <- plotframe %>%
  group_by(link) %>%
  summarize(min_date = min(startdatenum), max_date = as.numeric(as.Date("2024-01-01"))) %>%
  mutate(date = purrr::map2(min_date, max_date, seq)) %>%
  unnest(date)

#merge back in the jobs
jobdate <- unique(plotframe[,c("link", "startdatenum", "funccat_bn", "black")])
colnames(jobdate)[2] <- "date"
jobdate2 <- unique(plotframe[,c("link", "enddatenum", "next_funccat_bn", "black")])
jobdate2$enddatenum <- jobdate2$enddatenum + 1
colnames(jobdate2)[c(2,3)] <- c("date", "funccat_bn")
jobdate <- unique(rbind(jobdate, jobdate2))
date_sequence <- plyr::join(date_sequence, jobdate)
date_sequence <- date_sequence %>%
  fill(funccat_bn, black)

# Calculate days elapsed since the first startdatenum for each player
date_sequence$days_elapsed <- date_sequence$date - date_sequence$min_date

#subset on first 10 years
date_sequence <- subset.data.frame(date_sequence, days_elapsed <= 3650)

# Calculate total observations per race/spellnr/funccat_bn category
plotframe_count <- date_sequence %>%
  group_by(days_elapsed, black, funccat_bn) %>%
  count()

# Calculate shares by race/spellnr for each funccat_bn
plotframe_share <- plotframe_count %>%
  group_by(days_elapsed, black) %>%
  mutate(share = prop.table(n))

# Order the funccat_bn categories
funccat_bn_order <- rev(c("No Job", "Other staff", "Youth ROW", "Youth England", "Assistant manager ROW", "Assistant manager England", "Manager ROW", "Manager England"))

nblack <- paste0(length(unique(subset.data.frame(date_sequence, black == 1)$link)),
                 " (",
                 100*round(length(unique(subset.data.frame(date_sequence, black == 1)$link))/
                             races[2], 2),
                 "% of black players in data)")

nwhite <- paste0(length(unique(subset.data.frame(date_sequence, black == 0)$link)),
                 " (",
                 100*round(length(unique(subset.data.frame(date_sequence, black == 0)$link))/
                             races[1], 2),
                 "% of non-black players in data)")

# Convert funccat_bn to a factor with the specified order
plotframe_share$funccat_bn <- factor(plotframe_share$funccat_bn, levels = funccat_bn_order)

plotframe_share$type <- "Actual"
# SIM PREP 4: Simulate actual prevalence of job states SIMACT --------------------------------

#set seeds for replicability
seeds <- seq(1,reps,1)

#captures simulations here
statelist <- list()
countlist <- list()

starter <- subset.data.frame(date_sequence_sim, days_elapsed == 15)

for(h in 1:reps){
  
  #set seed
  set.seed(seeds[h])
  
  #grab starting data
  start <- starter
  table(start$funccat_bn)
  
  #make alternative black variable
  start$black_alt <- start$black
  
  #initialize objects
  states <- data.frame()
  counts <- data.frame()
  for(i in seq(15, 3650, 14)){
    print(paste("SIMACT", h,i))
    start$days_elapsed <- i
    
    #simulate actual
    #probabilities
    probs <- rbind(predict(simreg_multi, newdata = start, type = "probs"))[, jobstates]
    start$newstate <- apply(probs, 1, function(x) sample(jobstates, size = 1, prob = x))
    
    states_alt <- start[,c("link", "black", "black_alt", "days_elapsed", "newstate")]
    states <- states[,setdiff(colnames(states), "counter")]
    states <- rbind(states, states_alt)
    
    #count days spend in each category up until this point
    states$black <- states$black_alt
    states$counter <- 1
    
    pcollap <- collap(states, counter ~ link + newstate + black, FUN = sum)
    
    #add back categories with zero count
    all_combinations <- expand.grid(link = unique(pcollap$link),
                                    newstate = unique(pcollap$newstate))
    
    result <- merge(all_combinations, pcollap, by = c("link", "newstate"), all.x = TRUE)
    
    #fill missing race
    race <- subset.data.frame(unique(result[,c("link", "black")]), !is.na(black))
    result <- plyr::join(result[setdiff(colnames(result), "black")], race)
    
    # Replace missing values of "counter" with zeros
    result[is.na(result)] <- 0
    
    result <- collap(result, counter ~ newstate + black, FUN = c(mean))
    
    result$days_elapsed <- i
    
    counts <- rbind(counts, result)
    
    #move to new state for next simulation
    start$funccat_bn_lag <- start$newstate
    
  }
  
  #set back to null if needed
  plotframe_count <- NULL
  plotframe_share_simact <- NULL
  counter_simact <- NULL
  
  states$black <- states$black_alt
  colnames(states)[5] <- "funccat_bn"
  
  # Calculate total observations per race/spellnr/funccat_bn category
  plotframe_count <- states %>%
    group_by(days_elapsed, black, funccat_bn) %>%
    count()
  
  # Calculate shares by race/spellnr for each funccat_bn
  plotframe_share_simact <- plotframe_count %>%
    group_by(days_elapsed, black) %>%
    mutate(share = prop.table(n))
  
  # Order the funccat_bn categories
  funccat_bn_order <- rev(c("No Job", "Other staff", "Youth ROW", "Youth England", "Assistant manager ROW", "Assistant manager England", "Manager ROW", "Manager England"))
  
  # Convert funccat_bn to a factor with the specified order
  plotframe_share_simact$funccat_bn <- factor(plotframe_share_simact$funccat_bn, levels = funccat_bn_order)
  
  plotframe_share_simact$type <- "Simulated Actual"
  
  statelist[[h]] <- plotframe_share_simact
  
  counter_simact <- collap(counts, counter ~ newstate + black + days_elapsed, FUN =  mean)
  
  countlist[[h]] <- counter_simact
  
}

plotframe_share_simact <- do.call(rbind, statelist)

plotframe_share_simact <- collap(plotframe_share_simact, share + n ~ black + days_elapsed + funccat_bn + type, FUN = mean)

statelist_simact <- statelist

countlist_simact <- do.call(rbind, countlist)

countlist_simact <- collap(countlist_simact, counter ~ newstate + black + days_elapsed, FUN = mean)

countlist_simact_all <- countlist

# SIM PREP 5: Simulate as if enter as other race SIMEAO --------------------------------------


#set seeds for replicability
seeds <- seq(1,reps,1)

#captures simulations here
statelist <- list()
countlist <- list()

starter <- subset.data.frame(date_sequence_sim, days_elapsed == 15)

for(h in 1:reps){
  
  #set seed
  set.seed(seeds[h])
  
  #grab starting data
  start <- starter
  table(start$funccat_bn)
  
  #make alternative black variable
  start$black_alt <- start$black
  
  #simulate enter as other race
  start$black <- start$black_alt
  start$black <- ifelse(start$black_alt == 0, 1, start$black)
  start$black <- ifelse(start$black_alt == 1, 0, start$black)
  
  #simulate initial state
  #probabilities
  startprobs <- rbind(predict(startregmulti, newdata = start, type = "probs"))[, jobstates]
  start$startstate <- apply(startprobs, 1, function(x) sample(jobstates, size = 1, prob = x))
  start$funccat_bn_lag <- start$startstate
  
  #initialize objects
  states <- data.frame()
  counts <- data.frame()
  for(i in seq(15, 3650, 14)){
    print(paste("SIMEAO", h,i))
    start$days_elapsed <- i
    
    # initialize ever_employed tracker once (only during first iteration)
    if (i == 15) {
      start$ever_employed <- FALSE
    }
    
    # determine who has ever had a job
    start$ever_employed <- ifelse(start$funccat_bn_lag != "No Job", TRUE, start$ever_employed)
    
    # allow race switch only if never employed
    start$black <- start$black_alt
    start$black <- ifelse(!start$ever_employed & start$black_alt == 0, 1, start$black)
    start$black <- ifelse(!start$ever_employed & start$black_alt == 1, 0, start$black)
    
    #probabilities
    probs <- rbind(predict(simreg_multi, newdata = start, type = "probs"))[, jobstates]
    start$newstate <- apply(probs, 1, function(x) sample(jobstates, size = 1, prob = x))
    
    states_alt <- start[,c("link", "black", "black_alt", "days_elapsed", "newstate")]
    states <- states[,setdiff(colnames(states), "counter")]
    states <- rbind(states, states_alt)
    
    #count days spend in each category up until this point
    states$black <- states$black_alt
    states$counter <- 1
    
    pcollap <- collap(states, counter ~ link + newstate + black, FUN = sum)
    
    #add back categories with zero count
    all_combinations <- expand.grid(link = unique(pcollap$link),
                                    newstate = unique(pcollap$newstate))
    
    result <- merge(all_combinations, pcollap, by = c("link", "newstate"), all.x = TRUE)
    
    #fill missing race
    race <- subset.data.frame(unique(result[,c("link", "black")]), !is.na(black))
    result <- plyr::join(result[setdiff(colnames(result), "black")], race)
    
    # Replace missing values of "counter" with zeros
    result[is.na(result)] <- 0
    
    result <- collap(result, counter ~ newstate + black, FUN = c(mean))
    
    result$days_elapsed <- i
    
    counts <- rbind(counts, result)
    
    #move to new state for next simulation
    start$funccat_bn_lag <- start$newstate
    
  }
  
  #set back to null if needed
  plotframe_count <- NULL
  plotframe_share_simeao <- NULL
  counter_simeao <- NULL
  
  states$black <- states$black_alt
  colnames(states)[5] <- "funccat_bn"
  
  # Calculate total observations per race/spellnr/funccat category
  plotframe_count <- states %>%
    group_by(days_elapsed, black, funccat_bn) %>%
    count()
  
  # Calculate shares by race/spellnr for each funccat
  plotframe_share_simeao <- plotframe_count %>%
    group_by(days_elapsed, black) %>%
    mutate(share = prop.table(n))
  
  # Order the funccat categories
  funccat_bn_order <- rev(c("No Job", "Other staff", "Youth ROW", "Youth England", "Assistant manager ROW", "Assistant manager England", "Manager ROW", "Manager England"))
  
  # Convert funccat_bn to a factor with the specified order
  plotframe_share_simeao$funccat_bn <- factor(plotframe_share_simeao$funccat_bn, levels = funccat_bn_order)
  
  plotframe_share_simeao$type <- "Simulated Actual"
  
  statelist[[h]] <- plotframe_share_simeao
  
  counter_simeao <- collap(counts, counter ~ newstate + black + days_elapsed, FUN =  mean)
  
  countlist[[h]] <- counter_simeao
  
}

plotframe_share_simeao <- do.call(rbind, statelist)

plotframe_share_simeao <- collap(plotframe_share_simeao, share + n ~ black + days_elapsed + funccat_bn + type, FUN = mean)

countlist_simeao <- do.call(rbind, countlist)

countlist_simeao <- collap(countlist_simeao, counter ~ newstate + black + days_elapsed, FUN = mean)

# SIM PREP 6: Simulate progress as other race SIMPAO --------------------------------------

#set seeds for replicability
seeds <- seq(1,reps,1)

#captures simulations here
statelist <- list()
countlist <- list()

starter <- subset.data.frame(date_sequence_sim, days_elapsed == 15)

for(h in 1:reps){
  
  #set seed
  set.seed(seeds[h])
  
  #grab starting data
  start <- starter
  table(start$funccat_bn)
  
  #make alternative black variable
  start$black_alt <- start$black
  
  #initialize objects
  states <- data.frame()
  counts <- data.frame()
  for(i in seq(15, 3650, 14)){
    print(paste("SIMPAO", h,i))
    start$days_elapsed <- i
    
    # initialize ever_employed tracker once (only during first iteration)
    if (i == 15) {
      start$ever_employed <- FALSE
    }
    
    # determine who has ever had a job
    start$ever_employed <- ifelse(start$funccat_bn_lag != "No Job", TRUE, start$ever_employed)
    
    # allow race switch only if never employed
    start$black <- start$black_alt
    start$black <- ifelse(start$ever_employed & start$black_alt == 0, 1, start$black)
    start$black <- ifelse(start$ever_employed & start$black_alt == 1, 0, start$black)
    
    #probabilities
    probs <- rbind(predict(simreg_multi, newdata = start, type = "probs"))[, jobstates]
    start$newstate <- apply(probs, 1, function(x) sample(jobstates, size = 1, prob = x))
    start$funccat_bn_lag <- start$startstate
    
    states_alt <- start[,c("link", "black", "black_alt", "days_elapsed", "newstate")]
    states <- states[,setdiff(colnames(states), "counter")]
    states <- rbind(states, states_alt)
    
    #count days spend in each category up until this point
    states$black <- states$black_alt
    states$counter <- 1
    
    pcollap <- collap(states, counter ~ link + newstate + black, FUN = sum)
    
    #add back categories with zero count
    all_combinations <- expand.grid(link = unique(pcollap$link),
                                    newstate = unique(pcollap$newstate))
    
    result <- merge(all_combinations, pcollap, by = c("link", "newstate"), all.x = TRUE)
    
    #fill missing race
    race <- subset.data.frame(unique(result[,c("link", "black")]), !is.na(black))
    result <- plyr::join(result[setdiff(colnames(result), "black")], race)
    
    # Replace missing values of "counter" with zeros
    result[is.na(result)] <- 0
    
    result <- collap(result, counter ~ newstate + black, FUN = c(mean))
    
    result$days_elapsed <- i
    
    counts <- rbind(counts, result)
    
    #move to new state for next simulation
    start$funccat_bn_lag <- start$newstate
    
  }
  
  #set back to null if needed
  plotframe_count <- NULL
  plotframe_share_simpao <- NULL
  counter_simpao <- NULL
  
  states$black <- states$black_alt
  colnames(states)[5] <- "funccat_bn"
  
  # Calculate total observations per race/spellnr/funccat category
  plotframe_count <- states %>%
    group_by(days_elapsed, black, funccat_bn) %>%
    count()
  
  # Calculate shares by race/spellnr for each funccat
  plotframe_share_simpao <- plotframe_count %>%
    group_by(days_elapsed, black) %>%
    mutate(share = prop.table(n))
  
  # Order the funccat categories
  funccat_bn_order <- rev(c("No Job", "Other staff", "Youth ROW", "Youth England", "Assistant manager ROW", "Assistant manager England", "Manager ROW", "Manager England"))
  
  # Convert funccat_bn to a factor with the specified order
  plotframe_share_simpao$funccat_bn <- factor(plotframe_share_simpao$funccat_bn, levels = funccat_bn_order)
  
  plotframe_share_simpao$type <- "Simulated Actual"
  
  statelist[[h]] <- plotframe_share_simpao
  
  counter_simpao <- collap(counts, counter ~ newstate + black + days_elapsed, FUN =  mean)
  
  countlist[[h]] <- counter_simpao
  
}

plotframe_share_simpao <- do.call(rbind, statelist)

plotframe_share_simpao <- collap(plotframe_share_simpao, share + n ~ black + days_elapsed + funccat_bn + type, FUN = mean)

countlist_simpao <- do.call(rbind, countlist)

countlist_simpao <- collap(countlist_simpao, counter ~ newstate + black + days_elapsed, FUN = mean)

# SIM PREP 7: Simulate matrix of other race SIMOM --------------------------------------


#set seeds for replicability
seeds <- seq(1,reps,1)

#captures simulations here
statelist <- list()
countlist <- list()

starter <- subset.data.frame(date_sequence_sim, days_elapsed == 15)

for(h in 1:reps){
  
  #set seed
  set.seed(seeds[h])
  
  #grab starting data
  start <- starter
  table(start$funccat_bn)
  
  #make alternative black variable
  start$black_alt <- start$black
  
  #simulate enter as other race
  start$black <- start$black_alt
  start$black <- ifelse(start$black_alt == 0, 1, start$black)
  start$black <- ifelse(start$black_alt == 1, 0, start$black)
  
  #simulate initial state
  #probabilities
  startprobs <- rbind(predict(startregmulti, newdata = start, type = "probs"))[, jobstates]
  start$startstate <- apply(startprobs, 1, function(x) sample(jobstates, size = 1, prob = x))
  start$funccat_bn_lag <- start$startstate
  
  #initialize objects
  states <- data.frame()
  counts <- data.frame()
  for(i in seq(15, 3650, 14)){
    print(paste("SIMOM", h,i))
    start$days_elapsed <- i
    
    #simulate enter as other race
    start$black <- start$black_alt
    start$black <- ifelse(start$black_alt == 0, 1, start$black)
    start$black <- ifelse(start$black_alt == 1, 0, start$black)
    
    #probabilities
    probs <- rbind(predict(simreg_multi, newdata = start, type = "probs"))[, jobstates]
    start$newstate <- apply(probs, 1, function(x) sample(jobstates, size = 1, prob = x))
    start$funccat_bn_lag <- start$startstate
    
    states_alt <- start[,c("link", "black", "black_alt", "days_elapsed", "newstate")]
    states <- states[,setdiff(colnames(states), "counter")]
    states <- rbind(states, states_alt)
    
    #count days spend in each category up until this point
    states$black <- states$black_alt
    states$counter <- 1
    
    pcollap <- collap(states, counter ~ link + newstate + black, FUN = sum)
    
    #add back categories with zero count
    all_combinations <- expand.grid(link = unique(pcollap$link),
                                    newstate = unique(pcollap$newstate))
    
    result <- merge(all_combinations, pcollap, by = c("link", "newstate"), all.x = TRUE)
    
    #fill missing race
    race <- subset.data.frame(unique(result[,c("link", "black")]), !is.na(black))
    result <- plyr::join(result[setdiff(colnames(result), "black")], race)
    
    # Replace missing values of "counter" with zeros
    result[is.na(result)] <- 0
    
    result <- collap(result, counter ~ newstate + black, FUN = c(mean))
    
    result$days_elapsed <- i
    
    counts <- rbind(counts, result)
    
    #move to new state for next simulation
    start$funccat_bn_lag <- start$newstate
    
  }
  
  #set back to null if needed
  plotframe_count <- NULL
  plotframe_share_simom <- NULL
  counter_simom <- NULL
  
  states$black <- states$black_alt
  colnames(states)[5] <- "funccat_bn"
  
  # Calculate total observations per race/spellnr/funccat category
  plotframe_count <- states %>%
    group_by(days_elapsed, black, funccat_bn) %>%
    count()
  
  # Calculate shares by race/spellnr for each funccat
  plotframe_share_simom <- plotframe_count %>%
    group_by(days_elapsed, black) %>%
    mutate(share = prop.table(n))
  
  # Order the funccat categories
  funccat_bn_order <- rev(c("No Job", "Other staff", "Youth ROW", "Youth England", "Assistant manager ROW", "Assistant manager England", "Manager ROW", "Manager England"))
  
  # Convert funccat_bn to a factor with the specified order
  plotframe_share_simom$funccat_bn <- factor(plotframe_share_simom$funccat_bn, levels = funccat_bn_order)
  
  plotframe_share_simom$type <- "Simulated Actual"
  
  statelist[[h]] <- plotframe_share_simom
  
  counter_simom <- collap(counts, counter ~ newstate + black + days_elapsed, FUN =  mean)
  
  countlist[[h]] <- counter_simom
  
}

plotframe_share_simom <- do.call(rbind, statelist)

plotframe_share_simom <- collap(plotframe_share_simom, share + n ~ black + days_elapsed + funccat_bn + type, FUN = mean)

countlist_simom <- do.call(rbind, countlist)

countlist_simom <- collap(countlist_simom, counter ~ newstate + black + days_elapsed, FUN = mean)

# SIM PREP 8: Simulate matrix of other race and initial states SIMOMA ------------------------
# 
# 
# #set seeds for replicability
# seeds <- seq(1,reps,1)
# 
# #captures simulations here
# statelist <- list()
# countlist <- list()
# 
# starter <- subset.data.frame(date_sequence_sim, days_elapsed == 15)
# 
# for(h in 1:reps){
#   
#   #set seed
#   set.seed(seeds[h])
#   
#   #grab starting data
#   start <- starter
#   table(start$funccat)
#   
#   #make alternative black variable
#   start$black_alt <- start$black
#   
#   #simulate enter as other race
#   start$black <- start$black_alt
#   start$black <- ifelse(start$black_alt == 0, 1, start$black)
#   start$black <- ifelse(start$black_alt == 1, 0, start$black)
#   
#   #simulate initial state
#   #probabilities
#   start$oneprobstart <- as.numeric(predict(startreg1, newdata = start, type = "response"))
#   start$twoprobstart <- as.numeric(predict(startreg2, newdata = start, type = "response"))
#   start$threeprobstart <- as.numeric(predict(startreg3, newdata = start, type = "response"))
#   start$fourprobstart <- as.numeric(predict(startreg4, newdata = start, type = "response"))
#   start$fiveprobstart <- as.numeric(predict(startreg5, newdata = start, type = "response"))
#   start$sixprobstart <- as.numeric(predict(startreg6, newdata = start, type = "response"))
#   start$totprobstart <- start$oneprobstart + start$twoprobstart + start$threeprobstart + start$fourprobstart + start$fiveprobstart + start$sixprobstart
#   
#   #rescale
#   start$oneprobstart <- start$oneprobstart/start$totprobstart
#   start$twoprobstart <- start$twoprobstart/start$totprobstart
#   start$threeprobstart <- start$threeprobstart/start$totprobstart
#   start$fourprobstart <- start$fourprobstart/start$totprobstart
#   start$fiveprobstart <- start$fiveprobstart/start$totprobstart
#   start$sixprobstart <- start$sixprobstart/start$totprobstart
#   
#   start$startstate <- apply(start[, c("oneprobstart", "twoprobstart", "threeprobstart", "fourprobstart", "fiveprobstart", "sixprobstart")], 1, function(x) sample(jobstates, size = 1, prob = x))
#   start$funccat_lag <- start$startstate
#   
#   #initialize objects
#   states <- data.frame() 
#   counts <- data.frame()
#   for(i in seq(15, 3650, 14)){
#     print(paste("SIMOMA", h,i))
#     start$days_elapsed <- i
#     
#     #simulate enter as other race
#     start$black <- start$black_alt
#     start$black <- ifelse(start$black_alt == 0, 1, start$black)
#     start$black <- ifelse(start$black_alt == 1, 0, start$black)
#     
#     #probabilities
#     start$oneprob <- as.numeric(predict(simreg1, newdata = start, type = "response"))
#     start$twoprob <- as.numeric(predict(simreg2, newdata = start, type = "response"))
#     start$threeprob <- as.numeric(predict(simreg3, newdata = start, type = "response"))
#     start$fourprob <- as.numeric(predict(simreg4, newdata = start, type = "response"))
#     start$fiveprob <- as.numeric(predict(simreg5, newdata = start, type = "response"))
#     start$sixprob <- as.numeric(predict(simreg6, newdata = start, type = "response"))
#     start$totprob <- start$oneprob + start$twoprob + start$threeprob + start$fourprob + start$fiveprob + start$sixprob
#     
#     #rescale
#     start$oneprob <- start$oneprob/start$totprob
#     start$twoprob <- start$twoprob/start$totprob
#     start$threeprob <- start$threeprob/start$totprob
#     start$fourprob <- start$fourprob/start$totprob
#     start$fiveprob <- start$fiveprob/start$totprob
#     start$sixprob <- start$sixprob/start$totprob
#     
#     start$newstate <- apply(start[, c("oneprob", "twoprob", "threeprob", "fourprob", "fiveprob", "sixprob")], 1, function(x) sample(jobstates, size = 1, prob = x))
#     
#     states_alt <- start[,c("link", "black", "black_alt", "days_elapsed", "newstate")]
#     states <- states[,setdiff(colnames(states), "counter")]
#     states <- rbind(states, states_alt)
#     
#     #count days spend in each category up until this point
#     states$black <- states$black_alt     
#     states$counter <- 1
#     
#     pcollap <- collap(states, counter ~ link + newstate + black, FUN = sum)
#     
#     #add back categories with zero count
#     all_combinations <- expand.grid(link = unique(pcollap$link),
#                                     newstate = unique(pcollap$newstate))
#     
#     result <- merge(all_combinations, pcollap, by = c("link", "newstate"), all.x = TRUE)
#     
#     #fill missing race
#     race <- subset.data.frame(unique(result[,c("link", "black")]), !is.na(black))
#     result <- plyr::join(result[setdiff(colnames(result), "black")], race)
#     
#     # Replace missing values of "counter" with zeros
#     result[is.na(result)] <- 0
#     
#     result <- collap(result, counter ~ newstate + black, FUN = c(mean))
#     
#     result$days_elapsed <- i
#     
#     counts <- rbind(counts, result)
#     
#     #move to new state for next simulation
#     start$funccat_lag <- start$newstate
#     
#   }
#   
#   #set back to null if needed
#   plotframe_count <- NULL
#   plotframe_share_simoma <- NULL
#   counter_simoma <- NULL
#   
#   states$black <- states$black_alt
#   colnames(states)[5] <- "funccat"
#   
#   # Calculate total observations per race/spellnr/funccat category
#   plotframe_count <- states %>%
#     group_by(days_elapsed, black, funccat) %>%
#     count()
#   
#   # Calculate shares by race/spellnr for each funccat
#   plotframe_share_simoma <- plotframe_count %>%
#     group_by(days_elapsed, black) %>%
#     mutate(share = prop.table(n))
#   
#   # Order the funccat categories
#   funccat_order <- rev(c("No Job", "Job outside England", "Other staff", "Youth team management/development", "Assistant manager", "Manager"))
#   
#   #plotframe_share_simoma <- subset.data.frame(plotframe_share_simoma, funccat != "No Job")
#   
#   plotframe_share_simoma$funccat <- ifelse(plotframe_share_simoma$funccat == "Job outside UK", "Job outside England", plotframe_share_simoma$funccat)
#   plotframe_share_simoma$funccat <- ifelse(plotframe_share_simoma$funccat == "Youth team management/development", "Youth", plotframe_share_simoma$funccat)
#   
#   # Order the funccat categories
#   funccat_order <- rev(c("No Job", "Job outside England", "Other staff", "Youth", "Assistant manager", "Manager"))
#   
#   # Convert funccat to a factor with the specified order
#   plotframe_share_simoma$funccat <- factor(plotframe_share_simoma$funccat, levels = funccat_order)
#   
#   plotframe_share_simoma$type <- "Simulated Actual"
#   
#   statelist[[h]] <- plotframe_share_simoma
#   
#   counter_simoma <- collap(counts, counter ~ newstate + black + days_elapsed, FUN =  mean)
#   
#   countlist[[h]] <- counter_simoma
#   
# }
# 
# plotframe_share_simoma <- do.call(rbind, statelist)
# 
# plotframe_share_simoma <- collap(plotframe_share_simoma, share + n ~ black + days_elapsed + funccat + type, FUN = mean)
# 
# countlist_simoma <- do.call(rbind, countlist)
# 
# countlist_simoma <- collap(countlist_simoma, counter ~ newstate + black + days_elapsed, FUN = mean)

# SIM PREP 9: Actual days spent -----------------------------------------

set1 <- data.frame()
for(i in seq(15, 3650, 14)){
  #collapse to player level
  dater <- subset.data.frame(date_sequence_sim, days_elapsed <= i)
  
  dater$counter <- 1
  
  pcollap <- collap(dater, counter ~ link + funccat_bn + black, FUN = sum)
  
  #add back categories with zero count
  all_combinations <- expand.grid(link = unique(pcollap$link),
                                  funccat_bn = unique(pcollap$funccat_bn))
  
  result <- merge(all_combinations, pcollap, by = c("link", "funccat_bn"), all.x = TRUE)
  
  # Replace missing values of "counter" with zeros
  result$counter[is.na(result$counter)] <- 0
  
  #fill missing black values
  race <- subset.data.frame(unique(result[,c("link", "black")]), !is.na(black))
  result <- plyr::join(result[setdiff(colnames(result), "black")], race)
  if(i == 3650){
    resultsave <- result
  }
  
  #collapse for chart
  set <- collap(result, counter ~ funccat_bn + black, FUN = c(mean, sd))
  set$days_elapsed <- i
  set1 <- rbind(set1, set)
}

set1_days <- set1

# SAVE SIMULATION ENVIRONMENT SIM PREP 1-9 --------------------------------------------------------

#doubtful if should keep: jobdate, jobdate2
keepdata <- c("counter_simact", "countlist_simeao", "countlist_simpao", "countlist_simom", "countlist_simoma",
              "countlist_simact", "countlist_simact_all", "coachdata", "date_sequence_sim", "jobevol", "natframe",
              "playercareerlevel", "playerlevel", "plotframe", "plotframe_share", "plotframe_share_simact", "plotframe_share_simeao", "plotframe_share_simpao",
              "plotframe_share_simom", "plotframe_share_simoma", "regframe", "responsesperid", "set1_days", "statelist_simact", "tiercount")

keepvecs <- c("axislabsize", "cbp5", "cbp6", "cbp8", "custom_join", "cutthres", "facetsize", "geomtextsize",
              "inputwd", "datawd", "pkgwd", "insiders", "joborder", "jobs",
              "jobstates", "labsize", "legendtextsize", "legsize", "links", "outputwdfiles",
              "outputwdresults", "reps", "rmlinks", "textsize", "ticksize", "ukcountries", "ticks", "uksim")

regressions <- c("simreg_multi", "startregmulti")

keeps <- c(keepdata, keepvecs)
removes <- setdiff(ls(), keeps)
rm(list=setdiff(ls(), keeps))
setwd(datawd)

# The file name carries both the sample it was built from and the number of replications, so
# the three uksim settings do not overwrite one another and a short test run cannot overwrite
# the full one. HPS Simulation results.R loads whichever combination it is asked for.
simsample <- c("all",   # uksim = 0, all former players
               "uk",    # uksim = 1, UK nationals only
               "ms")[uksim + 1]
simenv    <- sprintf("data_environment_%s_%d.RData", simsample, reps)

# Save only the objects the simulation sections of HPS Racial Gaps.R actually use, with
# unused columns dropped. A full save.image() here runs to roughly 130 MB, of which the
# simulations touch under a tenth: regframe alone accounts for 427 MB in memory and is never
# referenced after the simulation load, since every regression result comes from
# analysis_environment.RData instead. Columns are cut to those the simulation code names, and
# dplyr grouping indices are stripped because for these frames the attribute is larger than
# the data it describes.
# What HPS Simulation results.R actually reads: the five share series drive every figure and
# the decomposition, and statelist_simact supplies the 1000 replication summaries the 2.5th
# and 97.5th percentile bands are computed from. date_sequence_sim is read only by APPENDIX
# TABLE C1, which runs off the "all" environment, so the uk and ms files leave it out --
# it is the single largest object in them and nothing would ever read it. playerlevel is not
# read at all, and dropping it also stops the simulation load from silently overwriting the
# playerlevel that the results preamble took from analysis_environment.RData.
sim_keep <- c("plotframe_share", "plotframe_share_simact", "plotframe_share_simeao",
              "plotframe_share_simom", "plotframe_share_simpao",
              "statelist_simact")
if (uksim == 0) sim_keep <- c(sim_keep, "date_sequence_sim")

sim_dss_cols <- c("link", "dob", "yob", "days_elapsed", "black",
                  "funccat_bn", "funccat_bn_lag", "Manager",
                  "position", "attacker", "midfielder", "defender", "goalkeeper",
                  "nationality", "nationality2", "nationality3", "natteam", "natapps",
                  "uk", "UK", "ROW", "ukonly", "ukplus", "foreign",
                  "sum_goals_tier1", "sum_assists_tier1", "sum_minplayed_tier1",
                  "sum_goals_tier2", "sum_assists_tier2", "sum_minplayed_tier2",
                  "sum_goals_lowtier", "sum_assists_lowtier", "sum_minplayed_lowtier")

sim_plain <- function(x) { x <- as.data.frame(x); attr(x, "groups") <- NULL; x }

# NOTE ON CLASSES. The plotting code rbinds plotframe_share, which by then carries p25 and
# p975, with plotframe_share_simact, which does not. That only works because both are
# grouped_df: rbind on a grouped_df fills the missing columns with NA, whereas on a plain
# data.frame or a tbl_df it errors with "numbers of columns of arguments do not match". The
# five plotframe_share objects are therefore left in their original class. The remaining
# objects are never rbound against a differently-shaped frame, so they can be flattened,
# which is where the space is saved anyway.
if ("date_sequence_sim" %in% sim_keep) {
  date_sequence_sim <- sim_plain(date_sequence_sim)
  date_sequence_sim <- date_sequence_sim[, intersect(sim_dss_cols, names(date_sequence_sim)),
                                         drop = FALSE]
}
statelist_simact <- lapply(statelist_simact, sim_plain)

save(list = intersect(sim_keep, ls()), file = file.path(datawd, simenv), compress = "xz")
cat("wrote", simenv, "-",
    round(file.size(file.path(datawd, simenv)) / 1024^2, 1), "MB\n")

 

# Leave the working directory where it started. Each script sets pkgwd from getwd(), so a
# script that ends inside data/ would make the next one in the sequence resolve every path
# one level down and create stray data/ and output/ folders there.
setwd(pkgwd)
