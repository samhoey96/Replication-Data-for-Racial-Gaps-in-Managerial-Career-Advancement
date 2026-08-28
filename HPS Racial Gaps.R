# Preamble ----------------------------------------------------------------


rm(list = ls())

# packages
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

datawd <- file.path(pkgwd, "data")

# analysis_environment.RData is produced by "HPS Data preparation.R"; run that first.
load(file.path(datawd, "analysis_environment.RData"))
# enggames.Rdata is not loaded: its only object, `england`, is never referenced.

# That file is a save.image(), so the load also restored the working directories of the
# machine that built it. Re-derive them from the current location before anything uses them,
# or output/ resolves to a path that only exists on that machine.
pkgwd  <- getwd()
datawd <- file.path(pkgwd, "data")
setwd(datawd)

# In a non-interactive run, printing a plot opens R's default device and leaves an Rplots.pdf
# in the working directory. Every figure here is written explicitly with ggsave, so the default
# device is not needed.
if (!interactive()) pdf(NULL)


# Export working directories
# One flat output folder for both tables and figures.
outputwd        <- file.path(pkgwd, "output")
# output/ is split three ways: numbered figures, numbered tables, and everything else the
# code produces that the paper does not report.
outputwdfigures <- file.path(outputwd, "figures")
outputwdtables  <- file.path(outputwd, "tables")
outputwdother   <- file.path(outputwd, "other")
for (.d in c(outputwdfigures, outputwdtables, outputwdother))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
rm(.d)

# A figure keeps its number only if the paper reports it; the unnumbered sample variants go
# to output/other, so the numbered folders match the manuscript one-for-one.
fig_dir <- function(name) if (grepl("^figure", name)) outputwdfigures else outputwdother

# Shared table-output helpers

num_fmt <- function(x) {
  abs_x <- abs(x)
  val <- if (abs_x > 0 & abs_x < 0.001) formatC(abs_x, format = "f", digits = 4)
  else formatC(round(abs_x, 3), format = "f", digits = 3)
  if (x < 0) paste0("$-$", val) else val
}

rsq_fmt <- function(x) {
  if (round(x, 4) < 0.001) formatC(x, format = "f", digits = 4)
  else formatC(round(x, 3), format = "f", digits = 3)
}

r2_wide <- function(x) {
  if (x < 0.0001) formatC(x, format = "f", digits = 5)
  else if (x < 0.001) formatC(x, format = "f", digits = 4)
  else formatC(round(x, 3), format = "f", digits = 3)
}

# Cox model cell extractor (coxph summary columns)
get_cell <- function(m, varname) {
  s <- summary(m)$coefficients
  if (varname %in% rownames(s)) {
    b  <- s[varname, "coef"]
    se <- s[varname, "se(coef)"]
    p  <- s[varname, "Pr(>|z|)"]
    st <- ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
    list(coef = paste0(num_fmt(b), "$^{", st, "}$"), se = paste0("(", num_fmt(se), ")"))
  } else {
    list(coef = "", se = "")
  }
}

# var_row uses models_cox (set in FIGURE 7 section before calling)
var_row <- function(label, varname) {
  cells <- lapply(models_cox, get_cell, varname = varname)
  n <- length(cells)
  cat(label, " & ", paste(sapply(cells, `[[`, "coef"), collapse = " & "), " \\\\ \n")
  cat("  & ", paste(sapply(cells, `[[`, "se"), collapse = " & "), " \\\\ \n")
  cat(paste(rep(" &", n), collapse = ""), "\\\\ \n")
}

# lm/felm cell extractor (lm summary columns)
get_cell_lm <- function(m, varname) {
  s <- summary(m)$coefficients
  if (varname %in% rownames(s)) {
    b  <- s[varname, "Estimate"]
    se <- s[varname, "Std. Error"]
    p  <- s[varname, "Pr(>|t|)"]
    st <- ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
    list(coef = paste0(num_fmt(b), ifelse(st == "", "", paste0("$^{", st, "}$"))),
         se   = paste0("(", num_fmt(se), ")"))
  } else {
    list(coef = "", se = "")
  }
}

# var_row_lm uses lms_deg (set in APPENDIX TABLE A2 before calling)
var_row_lm <- function(label, varname) {
  cells <- lapply(lms_deg, get_cell_lm, varname = varname)
  cat(label, " & ", paste(sapply(cells, `[[`, "coef"), collapse = " & "), " \\\\ \n")
  cat("  & ", paste(sapply(cells, `[[`, "se"), collapse = " & "), " \\\\ \n")
  cat("  & & & \\\\ \n")
}

# Generic lm cell extractor and 5-column row (used in EXTRA TABLE 1)
cell <- function(m, varname) {
  s <- summary(m)$coefficients
  if (varname %in% rownames(s)) {
    b  <- s[varname, "Estimate"]
    se <- s[varname, "Std. Error"]
    p  <- s[varname, "Pr(>|t|)"]
    st <- ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
    list(coef = paste0(num_fmt(b), ifelse(st == "", "", paste0("$^{", st, "}$"))),
         se   = paste0("(", num_fmt(se), ")"))
  } else {
    list(coef = "", se = "")
  }
}
empty <- list(coef = "", se = "")
row5 <- function(label, c1, c2, c3, c4, c5) {
  cat(label, " & ", c1$coef, " & ", c2$coef, " & ", c3$coef, " & ", c4$coef, " & ", c5$coef, " \\\\ \n")
  cat(" & ", c1$se, " & ", c2$se, " & ", c3$se, " & ", c4$se, " & ", c5$se, " \\\\ \n")
  cat(" & & & & & \\\\ \n")
}
row5c <- function(label, c1, c2, c3, c4, c5) {
  cat(label, " & ", c1$coef, " & ", c2$coef, " & ", c3$coef, " & ", c4$coef, " & ", c5$coef, " \\\\ \n")
  cat(" & ", c1$se, " & ", c2$se, " & ", c3$se, " & ", c4$se, " & ", c5$se, " \\\\ \n")
}

# Poisson/IRR cell extractor (fixest fepois); returns same $coef/$se structure as cell()
cell_irr <- function(m, varname) {
  ct <- summary(m)$coeftable
  if (varname %in% rownames(ct)) {
    b  <- ct[varname, "Estimate"]
    se <- ct[varname, "Std. Error"]
    p_col <- if ("Pr(>|z|)" %in% colnames(ct)) "Pr(>|z|)" else "Pr(>|t|)"
    p  <- ct[varname, p_col]
    st <- ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
    irr <- exp(b)
    fmt <- function(x) formatC(round(x, 3), format = "f", digits = 3)
    list(coef = paste0(fmt(irr), ifelse(st == "", "", paste0("$^{", st, "}$"))),
         se   = paste0("[", fmt(exp(b - 1.96*se)), ", ", fmt(exp(b + 1.96*se)), "]"))
  } else {
    list(coef = "", se = "")
  }
}

# color palettes (cbp7 not in .RData)
cbp7 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999")

# country classification (not saved in .RData, needed for transition matrices)
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
  "Venezuela", "Wales", "Yemen", "Zambia", "Zimbabwe"
)
country_categories <- rep(NA, length(countries))
country_categories[grep("England|Wales|Scotland|Jersey|Guernsey|United Kingdom|Northern Ireland", countries, ignore.case = TRUE)] <- "UK"
country_categories[grep("Iceland|Croatia|Italy|Luxembourg|Moldova|Serbia|Malta|Slovenia|Ukraine|Turkey|Switzerland|Finland|Germany|Armenia|Belarus|France|Norway|Ireland|Portugal|Poland|Belgium|Romania|Netherlands|Denmark|Cyprus|Sweden|Austria|Latvia|Estonia|Lithuania|Greece|Hungary|Slovakia|Czech Republic|Spain|Albania|Bosnia-Herzegovina|North Macedonia", countries, ignore.case = TRUE)] <- "Europe"
country_categories[grep("Zambia|Tanzania|Central African Republic|Congo|South Africa|Senegal|Cote d'Ivoire|Zimbabwe|Egypt|Morocco|Algeria|Nigeria|Ghana|Sierra Leone|Liberia|Cameroon|Angola|The Gambia|Guinea|Togo|Mali|DR Congo|Cape Verde|Tunisia|Mauritius|Uganda|Burundi|Seychelles|Mozambique|Kenya|Malawi|Somalia", countries, ignore.case = TRUE)] <- "Africa"
country_categories[is.na(country_categories)] <- "ROW"
countryclass <- data.frame(nationality = countries, natcat = country_categories)
rm(countries, country_categories)

# DATA PREP 1: Summary Statistics Data Preparation --------------------------------------------------------

#plyr::join in race codings
playerlevel <- plyr::join(playerlevel, responsesperid)
playerlevel$black <- ifelse(playerlevel$blackmaj == 1, 1, 0)

#filter out players with unclear race
playerlevel <- subset.data.frame(playerlevel, nomaj == 0)

coachdata$coached <- 1
playerlevelstats <- subset.data.frame(playerlevel, link %in% unique(jobevol$link))
playerlevelstats <- custom_join(playerlevelstats, unique(coachdata[, c("link", "coached")]))
playerlevelstats$coached <- ifelse(is.na(playerlevelstats$coached), 0, playerlevelstats$coached)

#observations by tier/race
table(playerlevelstats$black)
table(playerlevelstats$onetier, playerlevelstats$black)[2,]
table(playerlevelstats$twotier, playerlevelstats$black)[2,]

#create played in tier variable
playerlevelstats$played_tier1 <- ifelse(playerlevelstats$sum_appearances_tier1 > 0,1,0)
playerlevelstats$played_tier2 <- ifelse(playerlevelstats$sum_appearances_tier2 > 0,1,0)
playerlevelstats$played_lowtier <- ifelse(playerlevelstats$sum_appearances_tier3 > 0 | playerlevelstats$sum_appearances_tier4 > 0 | playerlevelstats$sum_appearances_tier5 > 0 | playerlevelstats$sum_appearances_tier6 > 0,1,0)

#plyr::join in tier data
playerlevelstats <- custom_join(playerlevelstats, tiercount)
playerlevelstats[,(dim(playerlevelstats)[2]-6):dim(playerlevelstats)[2]][is.na(playerlevelstats[,(dim(playerlevelstats)[2]-6):dim(playerlevelstats)[2]])] <- 0
playerlevelstats$lowtier <- ifelse(playerlevelstats$tier3 == 1 | playerlevelstats$tier4 == 1 | playerlevelstats$tier5 == 1 | playerlevelstats$tier6 == 1, 1, 0)

#join in licence
licences <- unique(coachdata[,c("link", "licencedum", "coachlicence")])
licences <- licences %>%
  mutate(id = row_number()) %>%  # Add an ID column to uniquely identify rows
  pivot_wider(names_from = coachlicence, 
              values_from = coachlicence, 
              values_fill = list(coachlicence = 0), 
              values_fn = list(coachlicence = ~1)) %>%  # Convert to 1s
  dplyr::select(-id)  # Remove the ID column
colnames(licences) <- str_replace_all(colnames(licences), "\\ ", "_")
playerlevelstats <- custom_join(playerlevelstats, licences)

playerlevelstats$licencedum <- ifelse(is.na(playerlevelstats$licencedum), 0, playerlevelstats$licencedum)
# List of variables
variables <- c("UEFA_Pro_Licence", "UEFA_B_Licence", "UEFA_A_Licence", "Licence_trainer", "Goalkeeping_Coach_Licence", "A_Licence", "B_Licence")

# Loop through each variable and apply the ifelse condition
for (var in variables) {
  playerlevelstats[[var]] <- ifelse(is.na(playerlevelstats[[var]]), 0, playerlevelstats[[var]])
}
playerlevelstats$UEFA_A_Licence <- ifelse(playerlevelstats$A_Licence == 1, 1, playerlevelstats$UEFA_A_Licence )
playerlevelstats$UEFA_B_Licence <- ifelse(playerlevelstats$B_Licence == 1, 1, playerlevelstats$UEFA_B_Licence )
playerlevelstats$Licence_trainer <- ifelse(playerlevelstats$Goalkeeping_Coach_Licence == 1, 1, playerlevelstats$Licence_trainer)
playerlevelstats$Licence_trainer <- ifelse(playerlevelstats$UEFA_B_Licence == 1, 1, playerlevelstats$Licence_trainer)

#dual national
playerlevelstats$dualnational <- ifelse(playerlevelstats$natcat2 == "missing", 0, 1)
table(playerlevelstats$black, playerlevelstats$dualnational)

#ukness
playerlevelstats$ukonly <- ifelse(playerlevelstats$natcat == "UK" & playerlevelstats$dualnational == 0, 1, 0)
playerlevelstats$ukplus <- ifelse((playerlevelstats$natcat == "UK" & playerlevelstats$dualnational == 1) | (playerlevelstats$natcat2 == "UK" & playerlevelstats$dualnational == 1), 1, 0)
playerlevelstats$foreign <- ifelse(!(playerlevelstats$natcat == "UK" | playerlevelstats$natcat2 == "UK"), 1, 0)

#tiers coached in
managerlevel <- jobevol
managerlevel <- subset.data.frame(managerlevel, !(link %in% rmlinks))
managerlevel <- plyr::join(jobevol, responsesperid)
managerlevel$black <- ifelse(managerlevel$blackmaj == 1, 1, 0)
#filter out managers with unclear race
managerlevel <- subset.data.frame(managerlevel, nomaj == 0)
managerlevel <- unique(managerlevel[,c("link", "tier_end")])
managerlevel$job <- 1
managerlevel <- managerlevel %>%
  pivot_wider(
    names_from = tier_end,
    values_from = job,
    values_fill = 0  # Fill other values with 0
  )
colnames(managerlevel)[4:10] <- paste0("mantier", colnames(managerlevel)[4:10] )
managerlevel$lowtierman <- ifelse((managerlevel$mantier3 + managerlevel$mantier4 + managerlevel$mantier5 + managerlevel$mantier6) > 0, 1,0)
managerlevel <- managerlevel[,c(1,4,5,7:11)]
playerlevelstats <-  plyr::join(playerlevelstats, managerlevel)

# totals across ALL playing tiers (tier1 + tier2 + lowtier)
playerlevelstats$sum_appearances_all <- with(playerlevelstats,
                                             sum_appearances_tier1 + sum_appearances_tier2 + sum_appearances_lowtier
)
playerlevelstats$sum_goals_all <- with(playerlevelstats,
                                       sum_goals_tier1 + sum_goals_tier2 + sum_goals_lowtier
)
playerlevelstats$sum_assists_all <- with(playerlevelstats,
                                         sum_assists_tier1 + sum_assists_tier2 + sum_assists_lowtier
)
playerlevelstats$sum_minplayed_all <- with(playerlevelstats,
                                           sum_minplayed_tier1 + sum_minplayed_tier2 + sum_minplayed_lowtier
)

# (optional) overall per-game averages across all tiers
playerlevelstats$avg_goalspergame_all <- with(playerlevelstats,
                                              ifelse(sum_appearances_all > 0, sum_goals_all / sum_appearances_all, NA_real_)
)
playerlevelstats$avg_assistspergame_all <- with(playerlevelstats,
                                                ifelse(sum_appearances_all > 0, sum_assists_all / sum_appearances_all, NA_real_)
)
playerlevelstats$avg_minplayedpergame_all <- with(playerlevelstats,
                                                  ifelse(sum_appearances_all > 0, sum_minplayed_all / sum_appearances_all, NA_real_)
)

playerlevelstats$managerialeng <- ifelse(playerlevelstats$Manager_England == 1 | playerlevelstats$Assistant_manager_England == 1 | playerlevelstats$Youth_England == 1, 1, 0)
playerlevelstats$managerialrow <- ifelse(playerlevelstats$Manager_ROW== 1 | playerlevelstats$Assistant_manager_ROW == 1 | playerlevelstats$Youth_ROW == 1, 1, 0)
playerlevelstats$managerial <- ifelse(playerlevelstats$managerialeng == 1 | playerlevelstats$managerialrow == 1, 1, 0)


# -- Spell-level setup (managerlevel) --

managerlevel <- jobevol
managerlevel <- subset.data.frame(managerlevel, !(link %in% rmlinks))

managerlevel <- plyr::join(jobevol, responsesperid)
managerlevel$black <- ifelse(managerlevel$blackmaj == 1, 1, 0)

#filter out managers with unclear race
managerlevel <- subset.data.frame(managerlevel, nomaj == 0)

managerlevel$managerialeng <- ifelse(managerlevel$funccat == "Manager" | managerlevel$funccat == "Assistant manager" | managerlevel$funccat == "Youth team management/development", 1, 0)

#first job
firstjob <- subset.data.frame(managerlevel, funccat != "No Job")
firstjob <- firstjob %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
firstjob <- firstjob[,c("link", "black", "funccat")]
colnames(firstjob)[3] <- "funccatfirst"

# Pivot data to create indicators for the first job
firstjob$job <- 1
firstjob <- firstjob %>%
  pivot_wider(
    names_from = funccatfirst,
    values_from = job,
    values_fill = 0  # Fill other values with 0
  )
colnames(firstjob)[3:dim(firstjob)[2]] <- paste0("first", colnames(firstjob)[3:dim(firstjob)[2]])

managerlevel <- plyr::join(managerlevel, firstjob)


#time to first job
ttfj <- subset.data.frame(managerlevel, funccat != "No Job")
ttfj <- ttfj %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttfj$timetofirstjob <- ttfj$days_elapsed
ttfj <- ttfj[,c("link", "black", "timetofirstjob")]

managerlevel <- plyr::join(managerlevel, ttfj)

#time to manager
ttm <- subset.data.frame(managerlevel, funccat == "Manager")
ttm <- ttm %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttm$timetomanager <- ttm$days_elapsed
ttm <- ttm[,c("link", "black", "timetomanager")]

managerlevel <- plyr::join(managerlevel, ttm)

#time to assistant manager
ttam <- subset.data.frame(managerlevel, funccat == "Assistant manager")
ttam <- ttam %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttam$timetoassmanager <- ttam$days_elapsed
ttam <- ttam[,c("link", "black", "timetoassmanager")]

managerlevel <- plyr::join(managerlevel, ttam)

#time to youth
tty <- subset.data.frame(managerlevel, funccat == "Youth team management/development")
tty <- tty %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
tty$timetoyouth<- tty$days_elapsed
tty <- tty[,c("link", "black", "timetoyouth")]

managerlevel <- plyr::join(managerlevel, tty)

#time to other
tto <- subset.data.frame(managerlevel, funccat == "Other staff")
tto <- tto %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
tto$timetoother<- tto$days_elapsed
tto <- tto[,c("link", "black", "timetoother")]

managerlevel <- plyr::join(managerlevel, tto)

#time to outside uk
ttjo <- subset.data.frame(managerlevel, funccat == "Job outside UK")
ttjo <- ttjo %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttjo$timetoouk <- ttjo$days_elapsed
ttjo <- ttjo[,c("link", "black", "timetoouk")]

managerlevel <- plyr::join(managerlevel, ttjo)

#time to managerial england
ttme <- subset.data.frame(managerlevel, managerialeng == 1)
ttme <- ttme %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttme$timetomanagerialeng <- ttme$days_elapsed
ttme <- ttme[,c("link", "black", "timetomanagerialeng")]

managerlevel <- plyr::join(managerlevel, ttme)

#--- funccat_alt derivatives (Manager, Assistant manager, Youth - all countries) ---

#managerial indicator (funccat_alt)
managerlevel$managerial_alt <- ifelse(managerlevel$funccat_alt == "Manager" | managerlevel$funccat_alt == "Assistant manager" | managerlevel$funccat_alt == "Youth team management/development", 1, 0)

#first job (funccat_alt)
firstjob_alt <- subset.data.frame(managerlevel, funccat_alt != "No Job")
firstjob_alt <- firstjob_alt %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
firstjob_alt <- firstjob_alt[,c("link", "black", "funccat_alt")]
colnames(firstjob_alt)[3] <- "funccataltfirst"
firstjob_alt$job <- 1
firstjob_alt <- firstjob_alt %>%
  pivot_wider(
    names_from = funccataltfirst,
    values_from = job,
    values_fill = 0
  )
colnames(firstjob_alt)[3:dim(firstjob_alt)[2]] <- paste0("first_alt_", colnames(firstjob_alt)[3:dim(firstjob_alt)[2]])
managerlevel <- plyr::join(managerlevel, firstjob_alt)

#time to manager (funccat_alt - includes all countries)
ttm_alt <- subset.data.frame(managerlevel, funccat_alt == "Manager")
ttm_alt <- ttm_alt %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttm_alt$timetomanager_alt <- ttm_alt$days_elapsed
ttm_alt <- ttm_alt[,c("link", "black", "timetomanager_alt")]
managerlevel <- plyr::join(managerlevel, ttm_alt)

#time to assistant manager (funccat_alt - includes all countries)
ttam_alt <- subset.data.frame(managerlevel, funccat_alt == "Assistant manager")
ttam_alt <- ttam_alt %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttam_alt$timetoassmanager_alt <- ttam_alt$days_elapsed
ttam_alt <- ttam_alt[,c("link", "black", "timetoassmanager_alt")]
managerlevel <- plyr::join(managerlevel, ttam_alt)

#time to youth (funccat_alt - includes all countries)
tty_alt <- subset.data.frame(managerlevel, funccat_alt == "Youth team management/development")
tty_alt <- tty_alt %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
tty_alt$timetoyouth_alt <- tty_alt$days_elapsed
tty_alt <- tty_alt[,c("link", "black", "timetoyouth_alt")]
managerlevel <- plyr::join(managerlevel, tty_alt)

#time to managerial (funccat_alt - includes all countries)
ttme_alt <- subset.data.frame(managerlevel, managerial_alt == 1)
ttme_alt <- ttme_alt %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttme_alt$timetomanagerial_alt <- ttme_alt$days_elapsed
ttme_alt <- ttme_alt[,c("link", "black", "timetomanagerial_alt")]
managerlevel <- plyr::join(managerlevel, ttme_alt)

#--- funccat_bn derivatives (England/ROW split) ---

#first job (funccat_bn)
firstjob_bn <- subset.data.frame(managerlevel, funccat_bn != "No Job")
firstjob_bn <- firstjob_bn %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
firstjob_bn <- firstjob_bn[,c("link", "black", "funccat_bn")]
colnames(firstjob_bn)[3] <- "funccatbnfirst"
firstjob_bn$job <- 1
firstjob_bn <- firstjob_bn %>%
  pivot_wider(
    names_from = funccatbnfirst,
    values_from = job,
    values_fill = 0
  )
colnames(firstjob_bn)[3:dim(firstjob_bn)[2]] <- paste0("first_bn_", colnames(firstjob_bn)[3:dim(firstjob_bn)[2]])
managerlevel <- plyr::join(managerlevel, firstjob_bn)

#time to manager England (funccat_bn)
ttm_bn_eng <- subset.data.frame(managerlevel, funccat_bn == "Manager England")
ttm_bn_eng <- ttm_bn_eng %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttm_bn_eng$timetomanager_bn_eng <- ttm_bn_eng$days_elapsed
ttm_bn_eng <- ttm_bn_eng[,c("link", "black", "timetomanager_bn_eng")]
managerlevel <- plyr::join(managerlevel, ttm_bn_eng)

#time to manager ROW (funccat_bn)
ttm_bn_row <- subset.data.frame(managerlevel, funccat_bn == "Manager ROW")
ttm_bn_row <- ttm_bn_row %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttm_bn_row$timetomanager_bn_row <- ttm_bn_row$days_elapsed
ttm_bn_row <- ttm_bn_row[,c("link", "black", "timetomanager_bn_row")]
managerlevel <- plyr::join(managerlevel, ttm_bn_row)

#time to assistant manager England (funccat_bn)
ttam_bn_eng <- subset.data.frame(managerlevel, funccat_bn == "Assistant manager England")
ttam_bn_eng <- ttam_bn_eng %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttam_bn_eng$timetoassmanager_bn_eng <- ttam_bn_eng$days_elapsed
ttam_bn_eng <- ttam_bn_eng[,c("link", "black", "timetoassmanager_bn_eng")]
managerlevel <- plyr::join(managerlevel, ttam_bn_eng)

#time to assistant manager ROW (funccat_bn)
ttam_bn_row <- subset.data.frame(managerlevel, funccat_bn == "Assistant manager ROW")
ttam_bn_row <- ttam_bn_row %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
ttam_bn_row$timetoassmanager_bn_row <- ttam_bn_row$days_elapsed
ttam_bn_row <- ttam_bn_row[,c("link", "black", "timetoassmanager_bn_row")]
managerlevel <- plyr::join(managerlevel, ttam_bn_row)

#time to youth England (funccat_bn)
tty_bn_eng <- subset.data.frame(managerlevel, funccat_bn == "Youth England")
tty_bn_eng <- tty_bn_eng %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
tty_bn_eng$timetoyouth_bn_eng <- tty_bn_eng$days_elapsed
tty_bn_eng <- tty_bn_eng[,c("link", "black", "timetoyouth_bn_eng")]
managerlevel <- plyr::join(managerlevel, tty_bn_eng)

#time to youth ROW (funccat_bn)
tty_bn_row <- subset.data.frame(managerlevel, funccat_bn == "Youth ROW")
tty_bn_row <- tty_bn_row %>%
  group_by(link) %>%
  arrange(startdatenum) %>%
  slice(1)
tty_bn_row$timetoyouth_bn_row <- tty_bn_row$days_elapsed
tty_bn_row <- tty_bn_row[,c("link", "black", "timetoyouth_bn_row")]
managerlevel <- plyr::join(managerlevel, tty_bn_row)

#tiers worked


colnames(managerlevel) <- str_replace_all(colnames(managerlevel), "\\s", "_")
colnames(managerlevel) <- str_replace(colnames(managerlevel), "Youth_team_management/development", "Youth")

#job spell duration (funccat)
managerlevel$duration <- managerlevel$enddatenum - managerlevel$startdatenum
jobspelldur <- managerlevel[,c("link", "black", "funccat", "duration")]
jobspelldur$jobid <- seq(1, dim(jobspelldur)[1], 1)
jobspelldur <- jobspelldur %>%
  pivot_wider(
    names_from = funccat,
    values_from = duration
  )
colnames(jobspelldur)[4:dim(jobspelldur)[2]] <- paste0("dur", colnames(jobspelldur)[4:dim(jobspelldur)[2]])
colnames(jobspelldur) <- str_replace(colnames(jobspelldur), "Youth team management/development", "Youth")
colnames(jobspelldur) <- str_replace_all(colnames(jobspelldur), "\\s", "_")
jobspelldur$durany <- ifelse(is.na(jobspelldur$durNo_Job), rowSums(jobspelldur[,c(4,6,7,8,9)], na.rm =T), NA)

#job spell duration (funccat_alt)
jobspelldur_alt <- managerlevel[,c("link", "black", "funccat_alt", "duration")]
jobspelldur_alt$jobid <- seq(1, dim(jobspelldur_alt)[1], 1)
jobspelldur_alt <- jobspelldur_alt %>%
  pivot_wider(
    names_from = funccat_alt,
    values_from = duration
  )
colnames(jobspelldur_alt)[4:dim(jobspelldur_alt)[2]] <- paste0("dur_alt_", colnames(jobspelldur_alt)[4:dim(jobspelldur_alt)[2]])
colnames(jobspelldur_alt) <- str_replace(colnames(jobspelldur_alt), "Youth_team_management/development", "Youth")
colnames(jobspelldur_alt) <- str_replace_all(colnames(jobspelldur_alt), "\\s", "_")

#job spell duration (funccat_bn)
jobspelldur_bn <- managerlevel[,c("link", "black", "funccat_bn", "duration")]
jobspelldur_bn$jobid <- seq(1, dim(jobspelldur_bn)[1], 1)
jobspelldur_bn <- jobspelldur_bn %>%
  pivot_wider(
    names_from = funccat_bn,
    values_from = duration
  )
colnames(jobspelldur_bn)[4:dim(jobspelldur_bn)[2]] <- paste0("dur_bn_", colnames(jobspelldur_bn)[4:dim(jobspelldur_bn)[2]])
colnames(jobspelldur_bn) <- str_replace(colnames(jobspelldur_bn), "Youth_team_management/development", "Youth")
colnames(jobspelldur_bn) <- str_replace_all(colnames(jobspelldur_bn), "\\s", "_")


ms <- function(x) {
  paste0(formatC(mean(x, na.rm = TRUE), format = "f", digits = 2),
         " (", formatC(sd(x, na.rm = TRUE), format = "f", digits = 2), ")")
}
dt <- function(x, blk) {
  ok <- !is.na(x) & !is.na(blk)
  if (sum(ok) < 4 || length(unique(blk[ok])) < 2) return("--- (---)")
  tt <- t.test(x[ok] ~ blk[ok])
  paste0(formatC(round(diff(tt$estimate), 2), format = "f", digits = 2),
         " (", formatC(round(abs(tt$statistic), 2), format = "f", digits = 2), ")")
}
row_main <- function(label, x, blk, sp = "") {
  paste0(label, " & ", ms(x[blk == 0]), " & ", ms(x[blk == 1]), " & ", dt(x, blk), " \\\\ ", sp)
}
row_sub <- function(label, x, blk, sp = "") {
  paste0("  - ", label, " & ", ms(x[blk == 0]), " & ", ms(x[blk == 1]), " & ", dt(x, blk), " \\\\ ", sp)
}

# Panel A / D: player-level (playerlevelstats)
pA_blk <- playerlevelstats$black

# Panel B: one row per player who ever got a job
job_links_B  <- unique(managerlevel$link[as.character(managerlevel$funccat_alt) != "No Job"])
pnlB <- managerlevel[!duplicated(managerlevel$link) & managerlevel$link %in% job_links_B, ]
pB_blk <- pnlB$black

# Panel C: spell-level duration subsets
pC_mgr     <- managerlevel[managerlevel$funccat_alt == "Manager", ]
pC_mgr_eng <- managerlevel[managerlevel$funccat_bn  == "Manager England", ]
pC_mgr_row <- managerlevel[managerlevel$funccat_bn  == "Manager ROW", ]
pC_man     <- managerlevel[managerlevel$funccat_alt %in% c("Manager", "Assistant manager",
                                                           "Youth team management/development"), ]
pC_man_eng <- managerlevel[managerlevel$funccat_bn  %in% c("Manager England",
                                                           "Assistant manager England", "Youth England"), ]
pC_man_row <- managerlevel[managerlevel$funccat_bn  %in% c("Manager ROW",
                                                           "Assistant manager ROW", "Youth ROW"), ]
pC_all     <- managerlevel[managerlevel$funccat_alt != "No Job", ]

# DATA PREP 2: Propensity score matching → matchedlinks (used in Figures 2 & 3) -----------
regframe_tab2 <- playerlevelstats
regframe_tab2$yob <- substr(regframe_tab2$dob, nchar(regframe_tab2$dob)-3, nchar(regframe_tab2$dob))
regframe_tab2$uk <- ifelse(regframe_tab2$nationality %in% ukcountries | regframe_tab2$nationality2 %in% ukcountries | regframe_tab2$nationality3 %in% ukcountries, 1, 0)

glm1 <- glm(black ~ I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
              I(sum_minplayed_tier2/60) + sum_goals_tier2 + sum_assists_tier2 +
              I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
              attacker + midfielder + defender + ukonly + ukplus + foreign +
              natteam + natapps + as.factor(yob),
            family = binomial(link = "probit"), data = regframe_tab2)

set.seed(123)
rr1 <- Match(Y = regframe_tab2$coached, Tr = regframe_tab2$black,
             X = glm1$fitted, M = 1, ties = TRUE, replace = FALSE)

matchedlinks <- regframe_tab2$link[c(rr1$index.control, rr1$index.treated)]

# TABLE 1: Players summary statistics ----
header2 <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  " \\caption{Summary Statistics - Players}",
  "   \\label{tab:Summary Statistics Alt}",
  "\\footnotesize",
  "\\setlength\\tabcolsep{18pt}",
  "\\scalebox{0.9}{",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccc}",
  "\\hline",
  " & Non-Black & Black  & Difference \\\\",
  paste0(" & \\textit{(N = ", sum(pA_blk == 0, na.rm = TRUE),
         ")} & \\textit{(N = ", sum(pA_blk == 1, na.rm = TRUE),
         ")} & \\textit{(N = ", sum(!is.na(pA_blk)), ")}  \\\\"),
  " & Mean (std. dev.) & Mean (std. dev.) & Diff. (t-stat)\\\\",
  "\\hline",
  "\\\\[-1.8ex] ",
  "& \\multicolumn{3}{c}{Panel A: Playing Career Statistics}\\\\",
  "\\hline",
  "\\\\[-1.8ex] "
)

panelA2 <- c(
  row_main("Played in Tier 1",         playerlevelstats$played_tier1,              pA_blk),
  row_main("Played in Tier 2",         playerlevelstats$played_tier2,              pA_blk),
  row_main("Played in Lower Tiers",    playerlevelstats$played_lowtier,            pA_blk),
  row_main("Appearances",              playerlevelstats$sum_appearances_all,      pA_blk),
  row_main("Hours played",             playerlevelstats$sum_minplayed_all / 60,   pA_blk),
  row_main("Min. played per game",     playerlevelstats$avg_minplayedpergame_all, pA_blk),
  row_main("Goals",                    playerlevelstats$sum_goals_all,            pA_blk),
  row_main("Assists",                  playerlevelstats$sum_assists_all,          pA_blk),
  row_main("Goals per game",           playerlevelstats$avg_goalspergame_all,     pA_blk),
  row_main("Assists per game",         playerlevelstats$avg_assistspergame_all,   pA_blk)
)

panelB2_hdr <- c(
  "",
  "\\hline",
  "\\\\[-1.8ex] ",
  "& \\multicolumn{3}{c}{Panel B: Player Position Share}\\\\",
  "\\hline",
  "\\\\[-1.8ex] "
)

panelB2 <- c(
  row_main("  Attacker",    playerlevelstats$attacker,   pA_blk),
  row_main("  Midfielder",  playerlevelstats$midfielder,  pA_blk),
  row_main("  Defender",    playerlevelstats$defender,    pA_blk),
  row_main("  Goalkeeper",  playerlevelstats$goalkeeper,  pA_blk)
)

panelC2_hdr <- c(
  "\\hline",
  "\\\\[-1.8ex] ",
  "& \\multicolumn{3}{c}{Panel C: Nationality Share \\& National Team Stats}\\\\",
  "\\hline",
  "\\\\[-1.8ex] "
)

panelC2 <- c(
  row_main("UK Only",                    playerlevelstats$ukonly,  pA_blk),
  row_main("UK and Other",               playerlevelstats$ukplus,  pA_blk),
  row_main("Foreign",                    playerlevelstats$foreign, pA_blk),
  row_main("Played for National Team",   playerlevelstats$natteam, pA_blk),
  row_main("National Team Appearances",  playerlevelstats$natapps, pA_blk)
)

footer2 <- c(
  "\\hline",
  "\\\\[-1.8ex] ",
  "",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\t\t\\item[] \\scriptsize \\textit{Notes:}",
  paste0("This table presents players' career statistics and nationalities by race. Column 1 shows averages ",
         "for non-Black players, Column 2 for Black players, and Column 3 displays the difference between ",
         "the two along with the t-statistic to assess if the racial difference is significant. The nationality ",
         "shares include potential second nationalities. Statistics under Panel A and B refer to appearances in England."),
  "\t\t\\end{tablenotes}",
  "  \\end{threeparttable}",
  "  }",
  "\\end{table}"
)

tex2 <- c(header2, panelA2, panelB2_hdr, panelB2, panelC2_hdr, panelC2, footer2)
writeLines(tex2, con = file.path(outputwdtables, "table1_summary_stats_players.tex"))

# TABLE 2: Management summary statistics ----
header <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  " \\caption{Summary Statistics - Management}",
  "   \\label{tab:Summary Statistics man}",
  "\\footnotesize",
  "\\setlength\\tabcolsep{18pt}",
  "\\scalebox{0.8}{",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccc}",
  "\\hline",
  " & Non-Black & Black  & Difference \\\\",
  " & Mean (std. dev.) & Mean (std. dev.) & Diff. (t-stat)\\\\",
  "\\hline",
  "\\\\[-1.8ex] ",
  " & \\multicolumn{3}{c}{Panel A: Jobs ever obtained}\\\\",
  "\\hline",
  "\\\\[-1.8ex] "
)

panelA <- c(
  row_main("Manager",           playerlevelstats$Manager,                   pA_blk),
  row_sub("in England",         playerlevelstats$Manager_England,           pA_blk),
  row_sub("rest of world",      playerlevelstats$Manager_ROW,               pA_blk, "[1.8ex]"),
  row_main("Assistant manager", playerlevelstats$Assistant_manager,         pA_blk),
  row_sub("in England",         playerlevelstats$Assistant_manager_England, pA_blk),
  row_sub("rest of world",      playerlevelstats$Assistant_manager_ROW,     pA_blk, "[1.8ex]"),
  row_main("Youth team coach",  playerlevelstats$Youth_team_management,     pA_blk),
  row_sub("in England",         playerlevelstats$Youth_England,             pA_blk),
  row_sub("rest of world",      playerlevelstats$Youth_ROW,                 pA_blk, "[1.8ex]"),
  row_main("Other job",         playerlevelstats$Other_staff,               pA_blk),
  row_main("Any Job",           playerlevelstats$coached,                   pA_blk)
)

panelB_hdr <- c(
  "",
  "\\hline",
  "\\\\[-1.8ex] ",
  "& \\multicolumn{3}{c}{Panel B: First Job (if any job obtained)}\\\\",
  "\\hline",
  "\\\\[-1.8ex] "
)

panelB <- c(
  row_main("Manager",           pnlB$first_alt_Manager,                    pB_blk),
  row_sub("in England",         pnlB$first_bn_Manager_England,             pB_blk),
  row_sub("rest of world",      pnlB$first_bn_Manager_ROW,                 pB_blk, "[1.8ex]"),
  row_main("Assistant manager", pnlB$first_alt_Assistant_manager,          pB_blk),
  row_sub("in England",         pnlB$first_bn_Assistant_manager_England,   pB_blk),
  row_sub("rest of world",      pnlB$first_bn_Assistant_manager_ROW,       pB_blk, "[1.8ex]"),
  row_main("Youth team coach",  pnlB$first_alt_Youth,                      pB_blk),
  row_sub("in England",         pnlB$first_bn_Youth_England,               pB_blk),
  row_sub("rest of world",      pnlB$first_bn_Youth_ROW,                   pB_blk, "[1.8ex]"),
  row_main("Other job",         pnlB$first_alt_Other_staff,                pB_blk)
)

panelC_hdr <- c(
  "\\hline",
  "\\\\[-1.8ex] ",
  "& \\multicolumn{3}{c}{Panel C: Average Job Spell Duration in Days}\\\\",
  "\\hline",
  "\\\\[-1.8ex] "
)

panelC <- c(
  row_main("Manager",      pC_mgr$duration,     pC_mgr$black),
  row_sub("in England",    pC_mgr_eng$duration, pC_mgr_eng$black),
  row_sub("rest of world", pC_mgr_row$duration, pC_mgr_row$black, "[1.8ex]"),
  paste0("  Managerial job & ", ms(pC_man$duration[pC_man$black == 0]), " & ",
         ms(pC_man$duration[pC_man$black == 1]), " & ", dt(pC_man$duration, pC_man$black), " \\\\ "),
  row_sub("in England",    pC_man_eng$duration, pC_man_eng$black),
  row_sub("rest of world", pC_man_row$duration, pC_man_row$black, "[1.8ex]"),
  paste0("  Any Job & ", ms(pC_all$duration[pC_all$black == 0]), " & ",
         ms(pC_all$duration[pC_all$black == 1]), " & ", dt(pC_all$duration, pC_all$black), " \\\\ ")
)

panelD_hdr <- c(
  "\\hline",
  "\\\\[-1.8ex] ",
  "& \\multicolumn{3}{c}{Panel D: Coaching license obtained}\\\\",
  "\\hline",
  "\\\\[-1.8ex] "
)

panelD <- c(
  row_main("Any Coaching License", playerlevelstats$licencedum,       pA_blk),
  row_main("UEFA Pro License",     playerlevelstats$UEFA_Pro_Licence, pA_blk),
  row_main("UEFA A License",       playerlevelstats$UEFA_A_Licence,   pA_blk),
  row_main("Other License",        playerlevelstats$Licence_trainer,  pA_blk)
)

footer <- c(
  "\\hline",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\t\t\\item[] \\scriptsize \\textit{Notes:}",
  paste0("This table presents former players' managerial job roles (panel A), first jobs (panel B), ",
         "job duration (panel C), and coaching qualifications (panel D) by race. Column 1 shows averages ",
         "for non-Black players, Column 2 for Black players, and Column 3 displays the difference between ",
         "the two along with the t-statistic to assess if the racial difference is significant. Panel A and D ",
         "are based on observations for 4,770 former players, of which 3,933 are identified as non-Black and ",
         "837 as Black. Panel B and C are conditional on obtaining a job, such that the number of observations ",
         "decreases to 1,316 non-Black and 163 Black former players. The job category ``managerial job'' ",
         "combines managers, assistant managers and youth coaches."),
  "\t\t\\end{tablenotes}",
  "  \\end{threeparttable}",
  "  }",
  "\\end{table}"
)

tex <- c(header, panelA, panelB_hdr, panelB, panelC_hdr, panelC, panelD_hdr, panelD, footer)
writeLines(tex, con = file.path(outputwdtables, "table2_summary_stats_management.tex"))


# TABLE 3: chances of becoming executive (any job) -------------------------------------------

regframe_tab2$jobinuk <- ifelse(regframe_tab2$Manager_England == 1 |
                                  regframe_tab2$Assistant_manager_England == 1 |
                                  regframe_tab2$Youth_England == 1, 1, 0)
regframe_tab2$nother  <- ifelse(regframe_tab2$Manager == 1 | regframe_tab2$Assistant_manager == 1 | regframe_tab2$Youth_team_management == 1, 1, 0)

# any non-playing job located in England (pipeline in England + other staff in England).
# "Other staff" is not split by country in funccat_bn, so derive from raw coachdata.
alljobsinuk_links <- unique(coachdata$link[coachdata$country == "England" & coachdata$funccat != "No Job"])
regframe_tab2$alljobsinuk <- ifelse(regframe_tab2$link %in% alljobsinuk_links, 1, 0)

ctrl_ols <- "+ ukplus + foreign + attacker + midfielder + defender +
               I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
               I(sum_minplayed_tier2/60) + sum_goals_tier2  + sum_assists_tier2 +
               I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
               natteam + natapps | yob"

lm1 <- felm(I(coached)      ~ black, data = regframe_tab2)                                   # any job, few controls
lm2 <- felm(as.formula(paste("I(coached)    ~ black", ctrl_ols)), data = regframe_tab2)      # any job, all controls
lm3 <- felm(as.formula(paste("I(alljobsinuk) ~ black", ctrl_ols)), data = regframe_tab2)     # any job in England, all controls
lm4 <- felm(I(nother)       ~ black, data = regframe_tab2)                                   # pipeline anywhere, few controls
lm5 <- felm(as.formula(paste("I(nother)     ~ black", ctrl_ols)), data = regframe_tab2)      # pipeline anywhere, all controls
lm6 <- felm(as.formula(paste("I(jobinuk)    ~ black", ctrl_ols)), data = regframe_tab2)      # pipeline England, all controls

lms <- list(lm1, lm2, lm3, lm4, lm5, lm6)

.co <- function(v) round(table(v, regframe_tab2$black)[2, 1] / table(regframe_tab2$black)[1], 2)
corwhiteval <- c(.co(regframe_tab2$coached), .co(regframe_tab2$coached), .co(regframe_tab2$alljobsinuk),
                 .co(regframe_tab2$nother),  .co(regframe_tab2$nother),  .co(regframe_tab2$jobinuk))

.fe <- function(m, var) {
  s <- summary(m)$coefficients
  if (!var %in% rownames(s)) return(list(est = NA_real_, se = NA_real_, pv = NA_real_))
  list(est = s[var, 1], se = s[var, 2], pv = s[var, 4])
}
.fmte <- function(v) {
  if (is.na(v$est)) return("")
  st <- if (!is.na(v$pv) && v$pv < 0.01) "$^{***}$" else
    if (!is.na(v$pv) && v$pv < 0.05) "$^{**}$"  else
      if (!is.na(v$pv) && v$pv < 0.1)  "$^{*}$"   else ""
  e <- formatC(abs(v$est), format = "f", digits = 3)
  if (v$est < 0) paste0("$-$", e, st) else paste0(e, st)
}
.fmts <- function(v) {
  if (is.na(v$se)) return("")
  paste0("(", formatC(v$se, format = "f", digits = 3), ")")
}
.crow <- function(label, var, lms, blank = TRUE) {
  vs <- lapply(lms, .fe, var = var)
  el <- paste0(label, " & ", paste(sapply(vs, .fmte), collapse = " & "), " \\\\ ")
  sl <- paste0("  & ", paste(sapply(vs, .fmts), collapse = " & "), " \\\\ ")
  if (blank) c(el, sl, "  & & & & & \\\\ ") else c(el, sl)
}

n_obs <- sapply(lms, function(m) formatC(m$N, format = "d", big.mark = ","))
r2    <- sapply(lms, function(m) formatC(round(summary(m)$r.squared, 3), format = "f", digits = 3))

tex3 <- c(
  "\\begin{table}[!h]",
  "\\centering",
  " \\caption{Entry into Managerial Labour Market by Race}",
  "   \\label{tab:entryregsspecany}",
  "\\footnotesize",
  "\\scalebox{0.85}{",
  "\\begin{threeparttable}",
  "\\begin{tabular}{@{\\extracolsep{15pt}}lcccccc}",
  "\\\\[-1.8ex]\\hline",
  "\\\\[-1.8ex] Dependent variable: & \\multicolumn{3}{c}{Any Job in Football} & \\multicolumn{3}{c}{Pipeline Job} \\\\",
  "\\\\[-1.8ex] & All & All & in England & All & All & in England \\\\",
  "\\\\[-1.8ex] & (1) & (2) & (3) & (4) & (5) & (6)\\\\",
  "\\hline \\\\[-1.8ex]",
  .crow(" Black",        "black",      lms),
  .crow(" UK and Other", "ukplus",     lms),
  .crow(" Foreign",      "foreign",    lms),
  .crow(" Attacker",     "attacker",   lms),
  .crow(" Midfielder",   "midfielder", lms),
  .crow(" Defender",     "defender",   lms, blank = FALSE),
  "  & & & & & & \\\\ ",
  paste0("Average $Y_{i}$ Non-Black & ", paste(corwhiteval, collapse = " & "), " \\\\ "),
  "  & & & & & & \\\\ ",
  "\\hline \\\\[-1.8ex]",
  "Year of Birth FE & No & Yes & Yes & No & Yes & Yes \\\\",
  "Player performance & No & Yes & Yes & No & Yes & Yes \\\\",
  paste0("Observations & ", paste(n_obs, collapse = " & "), " \\\\"),
  paste0("R$^{2}$ & ", paste(r2, collapse = " & "), " \\\\"),
  "\\hline",
  "\\hline \\\\[-1.8ex]",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  paste0("\t\t\\item[] \\scriptsize \\textit{Notes:} This table presents linear probability regression results. ",
         "Columns (1)--(3) are for obtaining any non-playing job in football (anywhere, without controls in column 1 and with the full set of controls in column 2) and any non-playing job in England (column 3). ",
         "Columns (4)--(6) are for holding a pipeline (managerial) role of manager, assistant manager, or youth-team management: anywhere without controls (4), anywhere with controls (5), and in England (6). ",
         "Where included, controls are fixed effects for year of birth, nationality status (base category = UK-only citizen), playing-position controls (attacker, midfielder, defender; base = goalkeeper), national-team experience and ",
         "national-team appearances, and playing-performance statistics from Tier 1, Tier 2, and lower tiers ",
         "in England (hours played, goals, assists). Coefficients for the player career variables are omitted ",
         "from the output, but available on request. Significance is indicated as follows: ",
         "* p-value$<$ 0.1; ** p-value$<$ 0.05; *** p-value$<$ 0.01."),
  "\t\t\\end{tablenotes}",
  "  \\end{threeparttable}",
  "  }",
  "\\end{table}"
)

writeLines(tex3, con = file.path(outputwdtables, "table3_entry_regs.tex"))

# TABLE 3 (UK subsample): entry regressions ------------------------------------
# The paper runs each analysis either with nationality controls or on a subsample of UK
# nationals. This is the latter for the entry regressions. The subsample is defined as
# foreign == 0, i.e. holding a UK nationality, using the same natcat-based family as the
# controls themselves so the two are internally consistent. Within this subsample foreign is
# identically zero and is dropped from the specification; ukplus still varies, separating
# dual nationals holding a UK passport from UK-only nationals, and is retained.
#
# This section is self-contained and runs on its own once DATA PREP 2 has been executed. The
# outcome variables and the table-formatting helpers below are also built in the TABLE 3
# section above; they are repeated here, identically, so this section does not depend on
# having run that one first.

regframe_tab2$jobinuk <- ifelse(regframe_tab2$Manager_England == 1 |
                                  regframe_tab2$Assistant_manager_England == 1 |
                                  regframe_tab2$Youth_England == 1, 1, 0)
regframe_tab2$nother  <- ifelse(regframe_tab2$Manager == 1 | regframe_tab2$Assistant_manager == 1 | regframe_tab2$Youth_team_management == 1, 1, 0)

# any non-playing job located in England (pipeline in England + other staff in England).
# "Other staff" is not split by country in funccat_bn, so derive from raw coachdata.
alljobsinuk_links <- unique(coachdata$link[coachdata$country == "England" & coachdata$funccat != "No Job"])
regframe_tab2$alljobsinuk <- ifelse(regframe_tab2$link %in% alljobsinuk_links, 1, 0)

.fe <- function(m, var) {
  s <- summary(m)$coefficients
  if (!var %in% rownames(s)) return(list(est = NA_real_, se = NA_real_, pv = NA_real_))
  list(est = s[var, 1], se = s[var, 2], pv = s[var, 4])
}
.fmte <- function(v) {
  if (is.na(v$est)) return("")
  st <- if (!is.na(v$pv) && v$pv < 0.01) "$^{***}$" else
    if (!is.na(v$pv) && v$pv < 0.05) "$^{**}$"  else
      if (!is.na(v$pv) && v$pv < 0.1)  "$^{*}$"   else ""
  e <- formatC(abs(v$est), format = "f", digits = 3)
  if (v$est < 0) paste0("$-$", e, st) else paste0(e, st)
}
.fmts <- function(v) {
  if (is.na(v$se)) return("")
  paste0("(", formatC(v$se, format = "f", digits = 3), ")")
}
.crow <- function(label, var, lms, blank = TRUE) {
  vs <- lapply(lms, .fe, var = var)
  el <- paste0(label, " & ", paste(sapply(vs, .fmte), collapse = " & "), " \\\\ ")
  sl <- paste0("  & ", paste(sapply(vs, .fmts), collapse = " & "), " \\\\ ")
  if (blank) c(el, sl, "  & & & & & \\\\ ") else c(el, sl)
}

regframe_tab2_uk <- subset.data.frame(regframe_tab2, foreign == 0)

ctrl_ols_uk <- "+ ukplus + attacker + midfielder + defender +
                  I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
                  I(sum_minplayed_tier2/60) + sum_goals_tier2  + sum_assists_tier2 +
                  I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
                  natteam + natapps | yob"

lm1uk <- felm(I(coached)       ~ black, data = regframe_tab2_uk)
lm2uk <- felm(as.formula(paste("I(coached)     ~ black", ctrl_ols_uk)), data = regframe_tab2_uk)
lm3uk <- felm(as.formula(paste("I(alljobsinuk) ~ black", ctrl_ols_uk)), data = regframe_tab2_uk)
lm4uk <- felm(I(nother)        ~ black, data = regframe_tab2_uk)
lm5uk <- felm(as.formula(paste("I(nother)      ~ black", ctrl_ols_uk)), data = regframe_tab2_uk)
lm6uk <- felm(as.formula(paste("I(jobinuk)     ~ black", ctrl_ols_uk)), data = regframe_tab2_uk)

lms_uk <- list(lm1uk, lm2uk, lm3uk, lm4uk, lm5uk, lm6uk)

.co_uk <- function(v) round(table(v, regframe_tab2_uk$black)[2, 1] /
                              table(regframe_tab2_uk$black)[1], 2)
corwhiteval_uk <- c(.co_uk(regframe_tab2_uk$coached), .co_uk(regframe_tab2_uk$coached),
                    .co_uk(regframe_tab2_uk$alljobsinuk), .co_uk(regframe_tab2_uk$nother),
                    .co_uk(regframe_tab2_uk$nother), .co_uk(regframe_tab2_uk$jobinuk))

n_obs_uk <- sapply(lms_uk, function(m) formatC(m$N, format = "d", big.mark = ","))
r2_uk    <- sapply(lms_uk, function(m) formatC(round(summary(m)$r.squared, 3),
                                               format = "f", digits = 3))

tex3uk <- c(
  "\\begin{table}[!h]",
  "\\centering",
  " \\caption{Entry into Managerial Labour Market by Race, UK Nationals Only}",
  "   \\label{tab:entryregsspecanyuk}",
  "\\footnotesize",
  "\\scalebox{0.85}{",
  "\\begin{threeparttable}",
  "\\begin{tabular}{@{\\extracolsep{15pt}}lcccccc}",
  "\\\\[-1.8ex]\\hline",
  "\\\\[-1.8ex] Dependent variable: & \\multicolumn{3}{c}{Any Job in Football} & \\multicolumn{3}{c}{Pipeline Job} \\\\",
  "\\\\[-1.8ex] & All & All & in England & All & All & in England \\\\",
  "\\\\[-1.8ex] & (1) & (2) & (3) & (4) & (5) & (6)\\\\",
  "\\hline \\\\[-1.8ex]",
  .crow(" Black",        "black",      lms_uk),
  .crow(" UK and Other", "ukplus",     lms_uk),
  .crow(" Attacker",     "attacker",   lms_uk),
  .crow(" Midfielder",   "midfielder", lms_uk),
  .crow(" Defender",     "defender",   lms_uk, blank = FALSE),
  "  & & & & & & \\\\ ",
  paste0("Average $Y_{i}$ Non-Black & ", paste(corwhiteval_uk, collapse = " & "), " \\\\ "),
  "  & & & & & & \\\\ ",
  "\\hline \\\\[-1.8ex]",
  "Year of Birth FE & No & Yes & Yes & No & Yes & Yes \\\\",
  "Player performance & No & Yes & Yes & No & Yes & Yes \\\\",
  paste0("Observations & ", paste(n_obs_uk, collapse = " & "), " \\\\"),
  paste0("R$^{2}$ & ", paste(r2_uk, collapse = " & "), " \\\\"),
  "\\hline",
  "\\hline \\\\[-1.8ex]",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  paste0("\t\t\\item[] \\scriptsize \\textit{Notes:} This table repeats Table \\ref{tab:entryregsspecany} ",
         "on the subsample of former players holding a UK nationality, rather than controlling for nationality ",
         "across the full sample. Columns (1)--(3) are for obtaining any non-playing job in football (anywhere, without controls in column 1 and with the full set of controls in column 2) and any non-playing job in England (column 3). ",
         "Columns (4)--(6) are for holding a pipeline (managerial) role of manager, assistant manager, or youth-team management: anywhere without controls (4), anywhere with controls (5), and in England (6). ",
         "Because the subsample conditions on holding a UK nationality, the foreign-only indicator is identically zero and is omitted; ",
         "\\textit{UK and Other} continues to separate dual nationals from UK-only nationals, who are the base category. ",
         "Where included, further controls are fixed effects for year of birth, playing-position controls (attacker, midfielder, defender; base = goalkeeper), national-team experience and ",
         "national-team appearances, and playing-performance statistics from Tier 1, Tier 2, and lower tiers ",
         "in England (hours played, goals, assists). Coefficients for the player career variables are omitted ",
         "from the output, but available on request. Significance is indicated as follows: ",
         "* p-value$<$ 0.1; ** p-value$<$ 0.05; *** p-value$<$ 0.01."),
  "\t\t\\end{tablenotes}",
  "  \\end{threeparttable}",
  "  }",
  "\\end{table}"
)

writeLines(tex3uk, con = file.path(outputwdtables, "tableC1_entry_regs_uk.tex"))

# TABLE 3 (OLS, Matched subsample): entry regressions --------------------------

local({
  ms <- regframe_tab2[regframe_tab2$link %in% matchedlinks, ]
  
  ctrl <- "+ ukplus + foreign + attacker + midfielder + defender +
             I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
             I(sum_minplayed_tier2/60) + sum_goals_tier2  + sum_assists_tier2 +
             I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
             natteam + natapps | yob"
  
  lm1 <- felm(I(coached)     ~ black, data = ms)                                   # any job, few controls
  lm2 <- felm(as.formula(paste("I(coached)    ~ black", ctrl)), data = ms)         # any job, all controls
  lm3 <- felm(as.formula(paste("I(alljobsinuk) ~ black", ctrl)), data = ms)        # any job in England, all controls
  lm4 <- felm(I(nother)      ~ black, data = ms)                                   # pipeline anywhere, few controls
  lm5 <- felm(as.formula(paste("I(nother)     ~ black", ctrl)), data = ms)         # pipeline anywhere, all controls
  lm6 <- felm(as.formula(paste("I(jobinuk)    ~ black", ctrl)), data = ms)         # pipeline England, all controls
  
  lms <- list(lm1, lm2, lm3, lm4, lm5, lm6)
  
  co <- function(v) round(table(v, ms$black)[2, 1] / table(ms$black)[1], 2)
  corwhiteval <- c(co(ms$coached), co(ms$coached), co(ms$alljobsinuk),
                   co(ms$nother),  co(ms$nother),  co(ms$jobinuk))
  
  n_obs <- sapply(lms, function(m) formatC(m$N, format = "d", big.mark = ","))
  r2    <- sapply(lms, function(m) formatC(round(summary(m)$r.squared, 3), format = "f", digits = 3))
  
  tex <- c(
    "\\begin{table}[!h]",
    "\\centering",
    " \\caption{Entry into Managerial Labour Market by Race --- Matched Subsample}",
    "   \\label{tab:entryregsspecanymatched}",
    "\\footnotesize",
    "\\scalebox{0.85}{",
    "\\begin{threeparttable}",
    "\\begin{tabular}{@{\\extracolsep{15pt}}lcccccc}",
    "\\\\[-1.8ex]\\hline",
    "\\\\[-1.8ex] Dependent variable: & \\multicolumn{3}{c}{Any Job in Football} & \\multicolumn{3}{c}{Pipeline Job} \\\\",
    "\\\\[-1.8ex] & All & All & in England & All & All & in England \\\\",
    "\\\\[-1.8ex] & (1) & (2) & (3) & (4) & (5) & (6)\\\\",
    "\\hline \\\\[-1.8ex]",
    .crow(" Black",        "black",      lms),
    .crow(" UK and Other", "ukplus",     lms),
    .crow(" Foreign",      "foreign",    lms),
    .crow(" Attacker",     "attacker",   lms),
    .crow(" Midfielder",   "midfielder", lms),
    .crow(" Defender",     "defender",   lms, blank = FALSE),
    "  & & & & & & \\\\ ",
    paste0("Average $Y_{i}$ Non-Black & ", paste(corwhiteval, collapse = " & "), " \\\\ "),
    "  & & & & & & \\\\ ",
    "\\hline \\\\[-1.8ex]",
    "Year of Birth FE & No & Yes & Yes & No & Yes & Yes \\\\",
    "Player performance & No & Yes & Yes & No & Yes & Yes \\\\",
    paste0("Observations & ", paste(n_obs, collapse = " & "), " \\\\"),
    paste0("R$^{2}$ & ", paste(r2, collapse = " & "), " \\\\"),
    "\\hline",
    "\\hline \\\\[-1.8ex]",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    paste0("\t\t\\item[] \\scriptsize \\textit{Notes:} This table replicates Table 3 on the propensity-score matched subsample. ",
           "Columns (1)--(3) are linear probability models for obtaining any non-playing job in football (anywhere, without controls in column 1 and with the full set of controls in column 2) and any non-playing job in England (column 3). ",
           "Columns (4)--(6) are for holding a pipeline (managerial) job: anywhere without controls (4), anywhere with controls (5), and in England (6). ",
           "Where included, controls are year-of-birth fixed effects, nationality status (base = UK-only), playing position (base = goalkeeper), national-team experience and appearances, and playing performance by tier. ",
           "Significance is indicated as follows: ",
           "* p-value$<$ 0.1; ** p-value$<$ 0.05; *** p-value$<$ 0.01."),
    "\t\t\\end{tablenotes}",
    "  \\end{threeparttable}",
    "  }",
    "\\end{table}"
  )
  writeLines(tex, con = file.path(outputwdtables, "tableB2_entry_regs_matched.tex"))
})

# FIGURE 2 + 3: Time bar charts + cumulative days ------------------------------------------------------

# ── Shared helpers ────────────────────────────────────────────────────────────

funccat_bn_order  <- rev(c("No Job", "Other staff", "Youth ROW", "Youth England",
                           "Assistant manager ROW", "Assistant manager England",
                           "Manager ROW", "Manager England"))
funccat_alt_order <- rev(c("No Job", "Other staff", "Youth team management/development",
                           "Assistant manager", "Manager"))

# Collapse England/ROW pairs into combined categories
bn_to_alt <- function(x) {
  x <- ifelse(x %in% c("Manager England", "Manager ROW"), "Manager", x)
  x <- ifelse(x %in% c("Assistant manager England", "Assistant manager ROW"), "Assistant manager", x)
  ifelse(x %in% c("Youth England", "Youth ROW"), "Youth team management/development", x)
}

# Build a daily date_sequence (10 years) from a prepared plotframe
build_date_sequence <- function(pf) {
  ds <- pf %>%
    group_by(link) %>%
    summarize(min_date = min(startdatenum),
              max_date = as.numeric(as.Date("2024-01-01")), .groups = "drop") %>%
    mutate(date = purrr::map2(min_date, max_date, seq)) %>%
    unnest(date)
  jd1 <- unique(pf[, c("link", "startdatenum", "funccat_bn", "black")])
  colnames(jd1)[2] <- "date"
  jd2 <- unique(pf[, c("link", "enddatenum", "next_funccat_bn", "black")])
  jd2$enddatenum <- jd2$enddatenum + 1
  colnames(jd2)[c(2, 3)] <- c("date", "funccat_bn")
  ds <- plyr::join(ds, unique(rbind(jd1, jd2)))
  ds <- ds %>% fill(funccat_bn, black)
  ds$days_elapsed <- ds$date - ds$min_date
  subset.data.frame(ds, days_elapsed <= 3650)
}

# Stacked daily share barplot saved as PNG
make_barplot_bn <- function(ds, nblack_label, nwhite_label, savename,
                            palette = cbp7, fill_label = "",
                            height = 2.4 * 0.7 * 127, dpi = 600, bg = NULL) {
  ps <- ds %>%
    group_by(days_elapsed, black, funccat_bn) %>%
    count() %>%
    group_by(days_elapsed, black) %>%
    mutate(share = prop.table(n)) %>%
    subset(funccat_bn != "No Job")
  ps$funccat_bn <- factor(ps$funccat_bn, levels = funccat_bn_order)
  p <- ggplot(ps, aes(x = days_elapsed, y = share, fill = funccat_bn)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size = 12), legend.position = "bottom",
          axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize),
          axis.title.x = element_text(size = axislabsize), axis.title.y = element_text(size = axislabsize),
          plot.title = element_text(hjust = 0.5, size = textsize),
          axis.ticks.length = unit(0.2, "cm"),
          legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
    geom_bar(stat = "identity", position = "stack", width = 1, alpha = 1) +
    facet_grid(~ black, labeller = labeller(black = c(
      "1" = paste("Black Former Players, n =", nblack_label),
      "0" = paste("Non-Black Former Players, n =", nwhite_label)))) +
    scale_x_continuous(name = "Years since end of playing career",
                       breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = "Share of Sample",
                       breaks = seq(0, 0.25, 0.05), labels = seq(0, 0.25, 0.05)) +
    scale_fill_manual(values = palette, name = fill_label) +
    geom_vline(xintercept = c(1, seq(0, 3650, 365)),
               linetype = "dashed", color = "white", linewidth = 0.2) +
    geom_text(data = subset(ps, days_elapsed %in% seq(182, 3650 - 182, 365)),
              aes(label = scales::percent(share, accuracy = 1)),
              position = position_stack(vjust = 0.5), size = geomtextsize)
  args <- list(savename, device = "pdf", plot = p,
               width = 2.2 * 127, height = height, dpi = dpi, units = "mm")
  if (!is.null(bg)) args$bg <- bg
  do.call(ggsave, args)
  invisible(p)
}

# Year-by-year collapse loop (funccat_bn); returns list(set, resultsave)
build_set_bn <- function(ds) {
  out <- data.frame(); counter <- 0; resultsave <- NULL
  for (i in seq(365, 3650, 365)) {
    dater         <- subset.data.frame(ds, days_elapsed <= i)
    dater$counter <- 1
    pcollap       <- collap(dater, counter ~ link + funccat_bn + black, FUN = sum)
    all_comb      <- expand.grid(link = unique(pcollap$link),
                                 funccat_bn = unique(pcollap$funccat_bn))
    result        <- merge(all_comb, pcollap, by = c("link", "funccat_bn"), all.x = TRUE)
    result$counter[is.na(result$counter)] <- 0
    race          <- subset.data.frame(unique(result[, c("link", "black")]), !is.na(black))
    result        <- plyr::join(result[setdiff(colnames(result), "black")], race)
    if (i == 3650) resultsave <- result
    counter       <- counter + 1
    set           <- collap(result, counter ~ funccat_bn + black, FUN = c(mean, sd))
    set$year      <- counter
    out           <- rbind(out, set)
  }
  list(set = out, resultsave = resultsave)
}

# Collapse bn set to alt (England + ROW summed)
set_bn_to_alt <- function(set_bn) {
  set_bn$funccat_alt <- bn_to_alt(as.character(set_bn$funccat_bn))
  set_bn %>%
    group_by(funccat_alt, black, year) %>%
    summarize(mean.counter = sum(mean.counter),
              sd.counter   = sqrt(sum(sd.counter^2)), .groups = "drop")
}

# Cumulative days bar chart saved as PDF
days_bar <- function(data, y_axis, y_limit, savename, legend_pos = "none") {
  # savename may arrive as a bare filename; anchor it so output never lands in the data dir
  if (!grepl("/", savename)) savename <- file.path(outputwdfigures, savename)
  p <- ggplot(data, aes(x = year, y = mean.counter, fill = as.factor(black))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(legend.position = legend_pos, plot.title = element_blank()) +
    scale_fill_manual(name = "Race:", values = c("#56B4E9", "#CC79A7"),
                      labels = c("Non-Black", "Black")) +
    scale_x_continuous(name = "Years elapsed since end of playing career",
                       breaks = seq(0, 10, 1), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = y_axis, limits = c(0, y_limit)) +
    geom_text(aes(label = round(mean.counter, 0)),
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3)
  ggsave(savename, device = "pdf", p,
         width = 0.8 * 127, height = 1.2 * 0.7 * 127, dpi = 1200, units = "mm")
  invisible(p)
}

# ── All players ───────────────────────────────────────────────────────────────

plotframe <- jobevol
plotframe <- subset.data.frame(plotframe, !(link %in% rmlinks))

#plyr::join in race codings
plotframe <- plyr::join(plotframe, responsesperid)
plotframe$black <- ifelse(plotframe$blackmaj == 1, 1, 0)
plotframe <- subset.data.frame(plotframe, nomaj == 0)

#define race counts (also used in matched sample %)
races <- table(unique(plotframe[, c("link", "black")])$black)
sum(races)

plotframe <- subset.data.frame(plotframe,
                               max.season <= as.numeric(as.Date("2013-06-30")) & black %in% c(0, 1))
plotframe$next_funccat_bn <- ifelse(
  plotframe$next_startspell - plotframe$enddatenum > 1, "No Job", plotframe$next_funccat_bn)

ds_all <- build_date_sequence(plotframe)

make_barplot_bn(ds_all,
                nblack_label = length(unique(ds_all$link[ds_all$black == 1])),
                nwhite_label = length(unique(ds_all$link[ds_all$black == 0])),
                savename = file.path(outputwdfigures, "figure4_career_progression_all.pdf"))

r1          <- build_set_bn(ds_all)
set1        <- r1$set
resultsave1 <- r1$resultsave
set1_alt    <- set_bn_to_alt(set1)

# Regressions (bn)
rs1 <- resultsave1
rs1$funccat_bn <- factor(rs1$funccat_bn, levels = funccat_bn_order)
rs1 <- within(rs1, funccat_bn <- relevel(funccat_bn, ref = length(funccat_bn_order)))
lm1_bn <- lm(I(counter - 3650) ~ black * funccat_bn, data = rs1)
summary(lm1_bn)

# Regressions (alt)
rs1_alt <- resultsave1 %>%
  mutate(funccat_alt = bn_to_alt(as.character(funccat_bn))) %>%
  group_by(link, funccat_alt, black) %>%
  summarize(counter = sum(counter), .groups = "drop")
rs1_alt$funccat_alt <- factor(rs1_alt$funccat_alt, levels = funccat_alt_order)
rs1_alt <- within(rs1_alt, funccat_alt <- relevel(funccat_alt, ref = length(funccat_alt_order)))
lm1_alt <- lm(I(counter - 3650) ~ black * funccat_alt, data = rs1_alt)
summary(lm1_alt)

# ── Matched sample ────────────────────────────────────────────────────────────

plotframe <- subset.data.frame(jobevol, !(link %in% rmlinks) & link %in% matchedlinks)
plotframe <- plyr::join(plotframe, responsesperid)
plotframe$black <- ifelse(plotframe$blackmaj == 1, 1, 0)
plotframe <- subset.data.frame(plotframe,
                               nomaj == 0 & max.season <= as.numeric(as.Date("2013-06-30")) & black %in% c(0, 1))
plotframe$next_funccat_bn <- ifelse(
  plotframe$next_startspell - plotframe$enddatenum > 1, "No Job", plotframe$next_funccat_bn)

ds_ms <- build_date_sequence(plotframe)

nblack_ms <- paste0(length(unique(ds_ms$link[ds_ms$black == 1])), " (",
                    100 * round(length(unique(ds_ms$link[ds_ms$black == 1])) / races[2], 2),
                    "% of black players in data)")
nwhite_ms <- paste0(length(unique(ds_ms$link[ds_ms$black == 0])), " (",
                    100 * round(length(unique(ds_ms$link[ds_ms$black == 0])) / races[1], 2),
                    "% of non-black players in data)")

# Matched sample uses collapsed (alt-style) category labels
ps_ms <- ds_ms %>%
  group_by(days_elapsed, black, funccat_bn) %>%
  count() %>%
  group_by(days_elapsed, black) %>%
  mutate(share = prop.table(n)) %>%
  subset(funccat_bn != "No Job")
ps_ms$funccat_bn <- ifelse(ps_ms$funccat_bn == "Job outside UK",
                           "Job outside England", ps_ms$funccat_bn)
ps_ms$funccat_bn <- ifelse(ps_ms$funccat_bn == "Youth team management/development",
                           "Youth", ps_ms$funccat_bn)
ps_ms$funccat_bn <- factor(ps_ms$funccat_bn,
                           levels = rev(c("No Job", "Job outside England", "Other staff",
                                          "Youth", "Assistant manager", "Manager")))
barplot_ms <- ggplot(ps_ms, aes(x = days_elapsed, y = share, fill = funccat_bn)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 12), legend.position = "bottom",
        axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize),
        axis.title.x = element_text(size = axislabsize), axis.title.y = element_text(size = axislabsize),
        plot.title = element_text(hjust = 0.5, size = textsize),
        axis.ticks.length = unit(0.2, "cm"),
        legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
  geom_bar(stat = "identity", position = "stack", width = 1, alpha = 1) +
  facet_grid(~ black, labeller = labeller(black = c(
    "1" = paste("Black Former Players, n =", nblack_ms),
    "0" = paste("Non-Black Former Players, n =", nwhite_ms)))) +
  scale_x_continuous(name = "Years since end of playing career",
                     breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
  scale_y_continuous(name = "Share of Sample",
                     breaks = seq(0, 0.25, 0.05), labels = seq(0, 0.25, 0.05)) +
  scale_fill_manual(values = cbp5, name = "Job England:") +
  geom_vline(xintercept = c(1, seq(0, 3650, 365)),
             linetype = "dashed", color = "white", linewidth = 0.2) +
  geom_text(data = subset(ps_ms, days_elapsed %in% seq(182, 3650 - 182, 365)),
            aes(label = scales::percent(share, accuracy = 1)),
            position = position_stack(vjust = 0.5), size = geomtextsize)
ggsave(file.path(outputwdother, "barplotday_jobcat_nomods_ms.pdf"), device = "pdf",
       barplot_ms, width = 2.2 * 127, height = 2 * 0.7 * 127, units = "mm")

r2          <- build_set_bn(ds_ms)
set2        <- r2$set
resultsave_ms <- r2$resultsave
set2_alt    <- set_bn_to_alt(set2)

resultsave_ms$funccat_bn <- factor(resultsave_ms$funccat_bn, levels = funccat_bn_order)
resultsave_ms <- within(resultsave_ms, funccat_bn <- relevel(funccat_bn, ref = length(funccat_bn_order)))
lm1 <- lm(I(counter - 3650) ~ black * funccat_bn, data = resultsave_ms)
summary(lm1)

# ── UK players ────────────────────────────────────────────────────────────────

plotframe <- custom_join(jobevol, natframe)
plotframe <- subset.data.frame(plotframe,
                               (nationality %in% ukcountries | nationality2 %in% ukcountries |
                                  nationality3 %in% ukcountries) & !(link %in% rmlinks))
plotframe <- plyr::join(plotframe, responsesperid)
plotframe$black <- ifelse(plotframe$blackmaj == 1, 1, 0)
plotframe <- subset.data.frame(plotframe,
                               nomaj == 0 & max.season <= as.numeric(as.Date("2013-06-30")) & black %in% c(0, 1))
plotframe$next_funccat_bn <- ifelse(
  plotframe$next_startspell - plotframe$enddatenum > 1, "No Job", plotframe$next_funccat_bn)

ds_uk <- build_date_sequence(plotframe)

make_barplot_bn(ds_uk,
                nblack_label = length(unique(ds_uk$link[ds_uk$black == 1])),
                nwhite_label = length(unique(ds_uk$link[ds_uk$black == 0])),
                savename = file.path(outputwdfigures, "figure4_career_progression_ukonly.pdf"), bg = "transparent")

r3          <- build_set_bn(ds_uk)
set3        <- r3$set
resultsave3 <- r3$resultsave
set3_alt    <- set_bn_to_alt(set3)

# Regressions (bn)
rs3 <- resultsave3
rs3$funccat_bn <- factor(rs3$funccat_bn, levels = funccat_bn_order)
rs3 <- within(rs3, funccat_bn <- relevel(funccat_bn, ref = length(funccat_bn_order)))
lm3_bn <- lm(I(counter - 3650) ~ black * funccat_bn, data = rs3)
summary(lm3_bn)

# Regressions (alt)
rs3_alt <- resultsave3 %>%
  mutate(funccat_alt = bn_to_alt(as.character(funccat_bn))) %>%
  group_by(link, funccat_alt, black) %>%
  summarize(counter = sum(counter), .groups = "drop")
rs3_alt$funccat_alt <- factor(rs3_alt$funccat_alt, levels = funccat_alt_order)
rs3_alt <- within(rs3_alt, funccat_alt <- relevel(funccat_alt, ref = length(funccat_alt_order)))
lm3_alt <- lm(I(counter - 3650) ~ black * funccat_alt, data = rs3_alt)
summary(lm3_alt)

# ── Cumulative days charts ────────────────────────────────────────────────────

# Manager anywhere (all / UK)
days_bar(subset(set1, funccat_bn == "Manager ROW"),
         "Cumulative Average Days as Manager", 150, "figure3_days_as_manager_row_all.pdf")
days_bar(collap(subset(set2, funccat_bn == "Manager ROW"), mean.counter ~ black + year, FUN = sum),
         "Cumulative Average Days as Manager", 150, "figureB2_days_as_manager_row_ms.pdf")
days_bar(subset(set3, funccat_bn == "Manager ROW"),
         "Cumulative Average Days as Manager", 150, "figure3_days_as_manager_row_uk.pdf")

# Manager England only (all / UK)
days_bar(subset(set1, funccat_bn == "Manager England"),
         "Cumulative Average Days as Manager England", 150, "figure3_days_as_manager_eng_all.pdf",
         legend_pos = c(0.15, 0.85))
days_bar(collap(subset(set2, funccat_bn == "Manager England"), mean.counter ~ black + year, FUN = sum),
         "Cumulative Average Days as Manager", 150, "figureB2_days_as_manager_eng_ms.pdf",
         legend_pos = c(0.15, 0.85))
days_bar(subset(set3, funccat_bn == "Manager England"),
         "Cumulative Average Days as Manager England", 150, "figure3_days_as_manager_eng_uk.pdf")

# The "any job" version of this chart is not reported in the paper and is not written.

# FIGURE 4: Transition Probabilities --------------------------------------

bnframe <- jobevol
bnframe <- subset.data.frame(bnframe, !(link %in% rmlinks))

bnframe <- plyr::join(bnframe, responsesperid)
bnframe$black <- ifelse(bnframe$blackmaj == 1, 1, 0)
bnframe <- subset.data.frame(bnframe, nomaj == 0)

# Recode next state to No Job if gap > 365 days
bnframe$next_funccat_bn <- ifelse(bnframe$next_startspell - bnframe$enddatenum > 365,
                                  "No Job", bnframe$next_funccat_bn)

# Conditional on a transition occurring
bnframe <- subset.data.frame(bnframe, !is.na(next_funccat_bn))

# Join nationality controls from playerlevelstats
bnframe <- plyr::join(bnframe, playerlevelstats[, c("link", "ukplus", "foreign")], by = "link")

# Category definitions (high to low status for panel ordering)
bn_categories <- c("Manager England", "Manager ROW",
                   "Assistant manager England", "Assistant manager ROW",
                   "Youth England", "Youth ROW", "Other staff", "No Job")
bn_labels <- c("Mgr Eng", "Mgr ROW", "Ass. Eng", "Ass. ROW",
               "Youth Eng", "Youth ROW", "Other", "No Job")

# 5-tier hierarchy (England/ROW collapsed) for promotion
bn_tier <- c("Manager England" = 5, "Manager ROW" = 5,
             "Assistant manager England" = 4, "Assistant manager ROW" = 4,
             "Youth England" = 3, "Youth ROW" = 3,
             "Other staff" = 2, "No Job" = 1)

# 64 LPMs: for each current state (subsetted), one per destination
# Controls: nationality dummies. Inference uses conventional OLS standard errors.
bn_coef <- matrix(NA, nrow = length(bn_categories), ncol = length(bn_categories),
                  dimnames = list(bn_labels, bn_labels))
bn_se   <- bn_coef
bn_p    <- bn_coef
bn_nev  <- bn_coef   # number of transitions into each destination
bn_nblk <- bn_coef   # number of Black movers into each destination
bn_n    <- integer(length(bn_categories))

transition_black_stats <- function(data, destination) {
  df <- data.frame(
    y       = data$next_funccat_bn == destination,
    black   = data$black,
    ukplus  = data$ukplus,
    foreign = data$foreign,
    link    = data$link
  )
  df <- df[complete.cases(df), ]
  n       <- nrow(df)
  n_event <- if (n > 0) sum(df$y) else 0                 # transitions into this destination
  n_black <- if (n > 0) sum(df$y & df$black == 1) else 0 # Black movers into this destination
  na_ret  <- c(coef = NA_real_, se = NA_real_, p = NA_real_, n_event = n_event, n_black_event = n_black)
  if (n == 0) return(c(coef = NA_real_, se = NA_real_, p = NA_real_, n_event = 0, n_black_event = 0))
  
  fit <- tryCatch(
    lm(y ~ black + ukplus + foreign, data = df),
    error = function(e) NULL
  )
  if (is.null(fit) || !("black" %in% names(coef(fit)))) return(na_ret)
  
  ct <- tryCatch(
    summary(fit)$coefficients,
    error = function(e) NULL
  )
  if (is.null(ct) || !("black" %in% rownames(ct))) return(na_ret)
  
  c(coef = ct["black", "Estimate"],
    se   = ct["black", "Std. Error"],
    p    = ct["black", "Pr(>|t|)"],
    n_event = n_event,
    n_black_event = n_black)
}

for (k in seq_along(bn_categories)) {
  sub <- subset.data.frame(bnframe, funccat_bn == bn_categories[k])
  bn_n[k] <- nrow(sub)
  for (j in seq_along(bn_categories)) {
    stats <- transition_black_stats(sub, bn_categories[j])
    bn_coef[k, j] <- stats["coef"]
    bn_se[k, j]   <- stats["se"]
    bn_p[k, j]    <- stats["p"]
    bn_nev[k, j]  <- stats["n_event"]
    bn_nblk[k, j] <- stats["n_black_event"]
  }
}

# "Any" panel: unconditional coefficients (all current states pooled)
any_dest_coef <- setNames(rep(NA_real_, length(bn_categories)), bn_labels)
any_dest_se   <- any_dest_coef
any_dest_p    <- any_dest_coef
any_dest_nev  <- any_dest_coef
any_dest_nblk <- any_dest_coef
any_n_total   <- nrow(bnframe)

for (j in seq_along(bn_categories)) {
  stats <- transition_black_stats(bnframe, bn_categories[j])
  any_dest_coef[j] <- stats["coef"]
  any_dest_se[j]   <- stats["se"]
  any_dest_p[j]    <- stats["p"]
  any_dest_nev[j]  <- stats["n_event"]
  any_dest_nblk[j] <- stats["n_black_event"]
}

# Build plot data frame
bn_plot_df <- expand.grid(current = bn_labels, destination = bn_labels,
                          stringsAsFactors = FALSE)
bn_plot_df$coef    <- as.vector(bn_coef)
bn_plot_df$se      <- as.vector(bn_se)
bn_plot_df$p       <- as.vector(bn_p)
bn_plot_df$n_event       <- as.vector(bn_nev)
bn_plot_df$n_black_event <- as.vector(bn_nblk)
bn_plot_df$ci_lo <- bn_plot_df$coef - 1.96 * bn_plot_df$se
bn_plot_df$ci_hi <- bn_plot_df$coef + 1.96 * bn_plot_df$se
bn_plot_df       <- subset(bn_plot_df, !is.na(coef))

# Drop "Other" panel (replaced by "Any") and impossible No Job transitions
bn_plot_df <- subset(bn_plot_df, current != "Other")
bn_plot_df <- subset(bn_plot_df, !(current == "No Job" & destination == "No Job"))

# "Any" panel destination rows
any_df <- data.frame(
  current     = "Any",
  destination = bn_labels,
  coef        = any_dest_coef,
  se          = any_dest_se,
  p           = any_dest_p,
  n_event     = any_dest_nev,
  n_black_event = any_dest_nblk,
  stringsAsFactors = FALSE
)
any_df$ci_lo <- any_df$coef - 1.96 * any_df$se
any_df$ci_hi <- any_df$coef + 1.96 * any_df$se
any_df       <- subset(any_df, !is.na(coef))

# Pooled pipeline origins (origins only): England pipeline and ROW pipeline
pooled_origin_df <- function(origins, label) {
  sub   <- subset.data.frame(bnframe, funccat_bn %in% origins)
  stats <- lapply(bn_categories, function(dest) transition_black_stats(sub, dest))
  d <- data.frame(
    current     = label,
    destination = bn_labels,
    coef        = sapply(stats, function(s) s["coef"]),
    se          = sapply(stats, function(s) s["se"]),
    p           = sapply(stats, function(s) s["p"]),
    n_event     = sapply(stats, function(s) s["n_event"]),
    n_black_event = sapply(stats, function(s) s["n_black_event"]),
    stringsAsFactors = FALSE
  )
  d$ci_lo <- d$coef - 1.96 * d$se
  d$ci_hi <- d$coef + 1.96 * d$se
  subset(d, !is.na(coef))
}
pipe_eng_df <- pooled_origin_df(c("Manager England", "Assistant manager England", "Youth England"), "Pipeline Eng")
pipe_row_df <- pooled_origin_df(c("Manager ROW", "Assistant manager ROW", "Youth ROW"), "Pipeline ROW")
pipe_all_df <- pooled_origin_df(
  c("Manager England", "Manager ROW",
    "Assistant manager England", "Assistant manager ROW",
    "Youth England", "Youth ROW"),
  "Pipeline Jobs"
)

bn_plot_df_all <- rbind(bn_plot_df, pipe_all_df, pipe_eng_df, pipe_row_df, any_df)
bn_plot_df <- rbind(bn_plot_df, pipe_eng_df, pipe_row_df, any_df)

# # Drop cells identified off too few Black movers into the destination (noisy estimates)
min_black_movers <- 10
bn_plot_df <- subset(bn_plot_df, n_black_event >= min_black_movers)

# Drop cells identified off too few Black movers into the destination (noisy estimates)
# min_movers <- 50
# bn_plot_df <- subset(bn_plot_df, n_event >= min_movers)


# Highlight estimates with conventional p-values below 5 percent.
bn_plot_df$sig   <- with(bn_plot_df, !is.na(p) & p < 0.05)

origin_levels <- c(bn_labels[bn_labels != "Other"], "Pipeline Eng", "Pipeline ROW", "Any")  # x-axis
dest_levels   <- bn_labels                                   # facet panels

bn_plot_df$current     <- factor(bn_plot_df$current,     levels = origin_levels)
bn_plot_df$destination <- factor(bn_plot_df$destination, levels = dest_levels)
bn_plot_df$type        <- "origin"

p_bn_transition <- ggplot(bn_plot_df, aes(x = current, y = coef,
                                          ymin = ci_lo, ymax = ci_hi,
                                          colour = sig,
                                          shape  = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(width = 0.25, linewidth = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c("FALSE" = "grey65", "TRUE" = "#c0392b"),
    labels = c("n.s.", "p<0.05"),
    na.translate = FALSE, name = NULL) +
  scale_shape_manual(
    values = c("FALSE" = 1, "TRUE" = 16),
    labels = c("n.s.", "p<0.05"),
    na.translate = FALSE, name = NULL) +
  facet_wrap(~ destination, ncol = 2, scales = "free_y",
             labeller = labeller(destination = function(x) {
               full_names <- c(
                 "Mgr Eng"   = "Manager England",
                 "Mgr ROW"   = "Manager ROW",
                 "Ass. Eng"  = "Assistant Manager England",
                 "Ass. ROW"  = "Assistant Manager ROW",
                 "Youth Eng" = "Youth England",
                 "Youth ROW" = "Youth ROW",
                 "Other"     = "Other staff",
                 "No Job"    = "No Job"
               )
               full_names[x]
             })) +
  scale_y_continuous(labels = scales::label_number(scale = 100, accuracy = 1)) +
  coord_flip() +
  labs(x = NULL, y = "Black coefficient (percentage points)") +
  theme_bw(base_size = 13) +
  theme(legend.position    = "bottom",
        strip.background   = element_rect(fill = "grey92", colour = NA),
        strip.text         = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 10))

print(p_bn_transition)
# This chart is not reported in the paper: the same model is shown as origin panels in
# Figure 5. The plot is still built above, only the saving is switched off.
# ggsave(file.path(outputwdfigures, "transition_bn_conditional_coefs.pdf"), p_bn_transition,
#        width = 9, height = 10, device = "pdf")

bn_plot_df_all$sig <- with(bn_plot_df_all, !is.na(p) & p < 0.05)

dest_levels <- bn_labels
panel_labels <- c(
  "Any"           = "Any",
  "Pipeline Jobs" = "Pipeline Jobs",
  "Pipeline Eng"  = "Pipeline England",
  "Pipeline ROW"  = "Pipeline Rest of World",
  "Mgr Eng"       = "Manager England",
  "Mgr ROW"       = "Manager Rest of World",
  "Ass. Eng"      = "Assistant Manager England",
  "Ass. ROW"      = "Assistant Manager Rest of World",
  "Youth Eng"     = "Youth England",
  "Youth ROW"     = "Youth Rest of World"
)

make_origin_plot <- function(data, origin_levels, filename, width, height) {
  plot_df <- subset(data, current %in% origin_levels)
  plot_df$current <- factor(plot_df$current, levels = origin_levels)
  plot_df$destination <- factor(plot_df$destination, levels = dest_levels)
  panel_n <- tapply(plot_df$n_event, plot_df$current, sum, na.rm = TRUE)
  panel_labels_fig <- panel_labels[origin_levels]
  names(panel_labels_fig) <- origin_levels
  panel_labels_fig[is.na(panel_labels_fig)] <- origin_levels[is.na(panel_labels_fig)]
  for (origin in names(panel_labels_fig)) {
    n_origin <- if (origin %in% names(panel_n)) panel_n[[origin]] else 0
    panel_labels_fig[[origin]] <- paste0(
      panel_labels_fig[[origin]],
      " (n = ", format(n_origin, big.mark = ",", scientific = FALSE), ")"
    )
  }
  
  p <- ggplot(plot_df, aes(x = destination, y = coef,
                           ymin = ci_lo, ymax = ci_hi,
                           colour = sig,
                           shape  = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_errorbar(width = 0.25, linewidth = 0.5) +
    geom_point(size = 2) +
    scale_colour_manual(
      values = c("FALSE" = "grey65", "TRUE" = "#c0392b"),
      labels = c("n.s.", "p<0.05"),
      na.translate = FALSE, name = NULL) +
    scale_shape_manual(
      values = c("FALSE" = 1, "TRUE" = 16),
      labels = c("n.s.", "p<0.05"),
      na.translate = FALSE, name = NULL) +
    facet_wrap(~ current, ncol = 2, scales = "free_y",
               labeller = labeller(current = as_labeller(panel_labels_fig))) +
    scale_y_continuous(labels = scales::label_number(scale = 100, accuracy = 1)) +
    coord_flip() +
    labs(x = NULL, y = "Black coefficient (percentage points)") +
    theme_bw(base_size = 13) +
    theme(legend.position    = "bottom",
          strip.background   = element_rect(fill = "grey92", colour = NA),
          strip.text         = element_text(face = "bold", size = 11),
          panel.grid.major.y = element_blank(),
          axis.text.y        = element_text(size = 10))
  
  print(p)
  ggsave(file.path(outputwdfigures, paste0(filename, ".pdf")),
         p, width = width, height = height, device = "pdf")
  p
}

main_origin_levels <- c("Any", "Pipeline Jobs",
                        "Pipeline Eng", "Pipeline ROW",
                        "Mgr Eng", "Mgr ROW",
                        "Ass. Eng", "Ass. ROW")
appendix_origin_levels <- c(main_origin_levels, "Youth Eng", "Youth ROW")

p_bn_transition_origin <- make_origin_plot(
  subset(bn_plot_df_all, n_black_event >= min_black_movers),
  main_origin_levels,
  "figure5_job_transitions",
  width = 9,
  height = 10
)

p_bn_transition_origin_unrestricted <- make_origin_plot(
  bn_plot_df_all,
  appendix_origin_levels,
  "figureA1_job_transitions_unrestricted",
  width = 9,
  height = 12
)

# FIGURE 7 + TABLE 4 + APPENDIX TABLE A4: Duration analysis firing --------------------------------------

# Spells are split by whether the club is in England or the rest of the world. Cumulative
# surprise is no longer available: it needed the England-only betting odds, while the games now
# come from the TM scrape covering all countries.

cox_none <- function(d) {
  s <- Surv(time = d$cumulgames, time2 = d$cumulgames_next, event = d$leftjob, type = "counting")
  coxph(s ~ black, data = d)
}
cox_perf <- function(d) {
  s <- Surv(time = d$cumulgames, time2 = d$cumulgames_next, event = d$leftjob, type = "counting")
  coxph(s ~ black + I(-cumulpoints_last10), data = d)
}
# No divisional fixed effects. League dummies over-fitted the ROW cells, and once the fixture
# scrape is merged in div holds 423 competition names rather than a clean league code, so
# there is no well-defined division to absorb. See the note in DATA PREP 8.
# Nationality follows the definition used elsewhere in the paper: ukplus and foreign enter,
# with UK-only as the omitted reference. Hours played in the English tiers are not controlled
# for; they were dropped because they are missing for part of the sample and were largely
# proxying the same career background as the nationality dummies.
cox_full <- function(d) {
  s <- Surv(time = d$cumulgames, time2 = d$cumulgames_next, event = d$leftjob, type = "counting")
  coxph(s ~ black + I(-cumulpoints_last10) + cumulgamesstart + ukplus + foreign, data = d)
}

# spells behind a model, counted on the rows the model actually kept
spell_count <- function(m, d) {
  keep <- seq_len(nrow(d))
  if (!is.null(m$na.action)) keep <- keep[-as.integer(m$na.action)]
  length(unique(d$jobspellid[keep]))
}

# --- Manager spells ---
regframe_mgr_eng <- subset.data.frame(regframe, funccat %in% c("Manager") & engspell == 1)
regframe_mgr_row <- subset.data.frame(regframe, funccat %in% c("Manager") & engspell == 0)

mgr_eng_cox1 <- cox_none(regframe_mgr_eng)
mgr_row_cox1 <- cox_none(regframe_mgr_row)
mgr_eng_cox3 <- cox_full(regframe_mgr_eng)
mgr_row_cox3 <- cox_full(regframe_mgr_row)
summary(mgr_eng_cox1); summary(mgr_row_cox1); summary(mgr_eng_cox3); summary(mgr_row_cox3)

# --- Assistant manager spells ---
regframe_ass_eng <- subset.data.frame(regframe, funccat %in% c("Assistant manager") & engspell == 1)
regframe_ass_row <- subset.data.frame(regframe, funccat %in% c("Assistant manager") & engspell == 0)

ass_eng_cox3 <- cox_full(regframe_ass_eng)
ass_row_cox3 <- cox_full(regframe_ass_row)
summary(ass_eng_cox3); summary(ass_row_cox3)

# --- 2x2 survival panel ---
# Curves are pulled out of survfit and drawn with one ggplot so the four panels share axes,
# styling and a single legend, which arranging four ggsurvplots does not give. Conditioning is
# on cumulative points only, held at each cell's sample mean, not on the full control set.
# Colours match the cumulative-days charts produced by days_bar earlier in this script.
surv_cells <- list(
  mgr_eng = list(d = regframe_mgr_eng, lab = "Manager England"),
  mgr_row = list(d = regframe_mgr_row, lab = "Manager ROW"),
  ass_eng = list(d = regframe_ass_eng, lab = "Asst. Manager England"),
  ass_row = list(d = regframe_ass_row, lab = "Asst. Manager ROW"))

panel_labels <- sapply(surv_cells, function(x)
  sprintf("%s (N = %s)", x$lab,
          format(length(unique(x$d$jobspellid)), big.mark = ",")))

survdat <- do.call(rbind, lapply(names(surv_cells), function(nm) {
  dd <- surv_cells[[nm]]$d
  nd <- data.frame(black = c(0, 1),
                   cumulpoints_last10 = rep(mean(dd$cumulpoints_last10, na.rm = TRUE), 2))
  f  <- survfit(cox_perf(dd), newdata = nd)
  do.call(rbind, lapply(1:2, function(k) data.frame(
    time  = c(0, f$time),
    surv  = c(1, f$surv[, k]),
    black = factor(c("Non-Black", "Black")[k], levels = c("Non-Black", "Black")),
    panel = panel_labels[[nm]])))
}))
survdat$panel <- factor(survdat$panel, levels = unname(panel_labels))

survplot2x2 <- ggplot(survdat, aes(x = time, y = surv, colour = black)) +
  geom_step(linewidth = 1.1) +
  facet_wrap(~ panel, nrow = 2) +
  coord_cartesian(xlim = c(0, 350), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 350, 100)) +
  scale_colour_manual(values = c("Non-Black" = "#56B4E9", "Black" = "#CC79A7")) +
  labs(x = "Time (games)", y = "Survival probability", colour = "Race") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"))

ggsave(file.path(outputwdfigures, "figure6_survival_2x2.pdf"), device = "pdf", survplot2x2,
       width = 1.5*1.2*127, height = 1.5*1.2*0.85*127, units = "mm")

# --- Combined table (6 columns: mgr eng/row no controls, mgr eng/row full, ass eng/row full) ---
models_cox  <- list(mgr_eng_cox1, mgr_row_cox1, mgr_eng_cox3, mgr_row_cox3, ass_eng_cox3, ass_row_cox3)
models_data <- list(regframe_mgr_eng, regframe_mgr_row, regframe_mgr_eng, regframe_mgr_row,
                    regframe_ass_eng, regframe_ass_row)
nobs_vec   <- sapply(models_cox, function(m) format(m$n, big.mark = ","))
rsq_vec    <- sapply(models_cox, function(m) rsq_fmt(summary(m)$rsq["rsq"]))
spell_vec  <- trimws(format(mapply(spell_count, models_cox, models_data), big.mark = ","))
evblack_vec <- sapply(models_data, function(d) format(sum(d$leftjob[d$black == 1]), big.mark = ","))

# The earlier merged builder for this table is no longer written: table4_hazard_model.tex
# below is the version the manuscript uses, and two files sharing the caption
# "Cox hazard model for leaving a job" invited inputting the wrong one.

# --- Appendix coverage table: share of eligible spells reaching the analysis ---
# The denominator is every spell that could in principle enter: manager and assistant manager
# spells in coachdata, starting 1990 or later, held by a coach whose race is coded, excluding
# national teams. A spell counts as covered when it appears in regframe, i.e. when at least ten
# of its games were matched. jobspellid is the row index of coachdata, assigned the same way in
# DATA PREP 8, so the two line up without needing a join key.

cov_dat <- coachdata
cov_dat$jobspellid <- seq_len(nrow(cov_dat))
cov_dat <- plyr::join(cov_dat, unique(responsesperid[, c("link", "nomaj")]), by = "link")
cov_dat$yr <- as.integer(format(as.Date(cov_dat$startdatenum, origin = "1970-01-01"), "%Y"))
cov_dat <- subset.data.frame(cov_dat,
                             funccat_alt %in% c("Manager", "Assistant manager") &
                               !is.na(nomaj) & nomaj == 0 & !is.na(yr) & yr >= 1990 &
                               !is.na(startdatenum) & !is.na(enddatenum) &
                               (is.na(country) | country != "international"))

cov_dat$covered <- cov_dat$jobspellid %in% unique(regframe$jobspellid)
cov_dat$reg     <- ifelse(!is.na(cov_dat$country) & cov_dat$country == "England",
                          "England", "ROW")
cov_dat$ctry    <- ifelse(is.na(cov_dat$country), "Unknown", cov_dat$country)

#countries below the threshold are pooled so the table stays readable
.cov_min <- 25
.cov_tot <- table(cov_dat$ctry)
cov_dat$ctry_grp <- ifelse(cov_dat$ctry %in% names(.cov_tot)[.cov_tot >= .cov_min],
                           cov_dat$ctry, "Other countries")

.cvcell <- function(d) {
  m <- d[d$funccat_alt == "Manager", ]
  a <- d[d$funccat_alt == "Assistant manager", ]
  data.frame(mgr_n = nrow(m), mgr_c = sum(m$covered),
             mgr_p = if (nrow(m) > 0) 100 * mean(m$covered) else NA_real_,
             ass_n = nrow(a), ass_c = sum(a$covered),
             ass_p = if (nrow(a) > 0) 100 * mean(a$covered) else NA_real_)
}
.cvn <- function(x) formatC(x, format = "d", big.mark = ",")
.cvp <- function(x) ifelse(is.na(x), "--",
                           paste0(formatC(x, format = "f", digits = 1), "\\%"))
.cvrow <- function(lab, cl, bold = FALSE) {
  v <- c(.cvn(cl$mgr_n), .cvn(cl$mgr_c), .cvp(cl$mgr_p),
         .cvn(cl$ass_n), .cvn(cl$ass_c), .cvp(cl$ass_p))
  if (bold) { lab <- paste0("\\textbf{", lab, "}"); v <- paste0("\\textbf{", v, "}") }
  paste0(lab, " & ", paste(v, collapse = " & "), " \\\\ ")
}

cov_grps  <- unique(cov_dat$ctry_grp)
cov_cells <- lapply(cov_grps, function(g) .cvcell(cov_dat[cov_dat$ctry_grp == g, ]))
names(cov_cells) <- cov_grps
#England first, remaining countries by size, pooled residual last
.cv_ord <- sapply(cov_grps, function(g) {
  if (g == "England") -Inf
  else if (g == "Other countries") Inf
  else -(cov_cells[[g]]$mgr_n + cov_cells[[g]]$ass_n)
})
cov_grps <- cov_grps[order(.cv_ord)]

tex_cov <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Duration Analysis Spell Coverage by Country}",
  "\\label{tab:coverage_merged}",
  "\\footnotesize",
  "\\setlength\\tabcolsep{4pt}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lcccccc}",
  "\\\\[-1.8ex]\\hline",
  "\\hline \\\\[-1.8ex]",
  " & \\multicolumn{3}{c}{Manager} & \\multicolumn{3}{c}{Assistant Manager} \\\\",
  "\\cline{2-4} \\cline{5-7}",
  "Country & Spells & Covered & \\% & Spells & Covered & \\% \\\\",
  "\\hline \\\\[-1.8ex]",
  sapply(cov_grps, function(g) .cvrow(g, cov_cells[[g]])),
  "\\hline \\\\[-1.8ex]",
  .cvrow("England (total)",       .cvcell(cov_dat[cov_dat$reg == "England", ]), bold = TRUE),
  .cvrow("Rest of world (total)", .cvcell(cov_dat[cov_dat$reg == "ROW", ]),     bold = TRUE),
  .cvrow("All spells",            .cvcell(cov_dat),                             bold = TRUE),
  "\\hline",
  "\\hline \\\\[-1.8ex]",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  paste0("\\item[] \\scriptsize \\textit{Notes:} Eligible spells are manager and assistant manager ",
         "job spells in \\texttt{coachdata} that start in 1990 or later and belong to a coach whose ",
         "race is coded. A spell counts as covered when at least ten of its games are found in the ",
         "game data, the threshold required to construct the cumulative points control. Countries ",
         "with fewer than ", .cov_min, " eligible spells are pooled into \\textit{Other countries}. ",
         "Rows are ordered by number of eligible spells, with England first and the pooled group ",
         "last. National team spells are excluded throughout: they are not club jobs, and youth ",
         "national sides do carry fixtures on Transfermarkt, so leaving them in would admit a ",
         "handful of spells that do not belong in the sample."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

writeLines(tex_cov, con = file.path(outputwdtables, "tableA4_spell_coverage.tex"))

# --- Interaction models: performance x race -----------------------------------
# Tests whether a dip in team performance raises the exit hazard more steeply for Black
# than for non-Black coaches. Cumulative points is mean-centred within each cell and still
# enters with a negative sign, so the Black main effect is the racial gap at that cell's
# average performance rather than at the meaningless value of zero cumulative points, and
# the interaction is the additional effect of a one-point shortfall.
#
# These models are REQUIRED: they supply columns 3 to 6 of the main text table built below.
# Do not delete this block. The standalone appendix table that follows is optional.

.int_cells <- list(
  mgr_eng = regframe_mgr_eng, mgr_row = regframe_mgr_row,
  ass_eng = regframe_ass_eng, ass_row = regframe_ass_row)

models_int <- lapply(.int_cells, function(d) {
  d$points_c <- -(d$cumulpoints_last10 - mean(d$cumulpoints_last10, na.rm = TRUE))
  s <- Surv(time = d$cumulgames, time2 = d$cumulgames_next, event = d$leftjob,
            type = "counting")
  coxph(s ~ black * points_c + cumulgamesstart + ukplus + foreign, data = d)
})

# --- Main text table: no-controls and interaction specs side by side ----------
# This is the table the manuscript inputs (tab:hazard_final). Columns 1 and 2 are the
# uncontrolled manager models; columns 3 to 6 are the interaction models, so the Black
# coefficient there is the racial gap at average performance rather than an average across
# performance levels. It reuses models already fitted above rather than refitting them.

models_fin <- list(mgr_eng_cox1, mgr_row_cox1,
                   models_int$mgr_eng, models_int$mgr_row,
                   models_int$ass_eng, models_int$ass_row)
data_fin   <- list(regframe_mgr_eng, regframe_mgr_row,
                   regframe_mgr_eng, regframe_mgr_row,
                   regframe_ass_eng, regframe_ass_row)

nobs_fin  <- sapply(models_fin, function(m) format(m$n, big.mark = ","))
rsq_fin   <- sapply(models_fin, function(m) rsq_fmt(summary(m)$rsq["rsq"]))
spell_fin <- trimws(format(mapply(spell_count, models_fin, data_fin), big.mark = ","))
evbl_fin  <- sapply(data_fin, function(d) format(sum(d$leftjob[d$black == 1]),
                                                 big.mark = ","))

.finrow <- function(label, varname, blank = TRUE) {
  cs <- lapply(models_fin, get_cell, varname = varname)
  out <- c(paste0(label, " & ", paste(sapply(cs, `[[`, "coef"), collapse = " & "), " \\\\ "),
           paste0("  & ",      paste(sapply(cs, `[[`, "se"),   collapse = " & "), " \\\\ "))
  if (blank) c(out, " & & & & & & \\\\ ") else out
}

tex_fin <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Cox hazard model for leaving a job}",
  "\\label{tab:hazard_final}",
  "\\footnotesize",
  "\\setlength\\tabcolsep{3pt}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lcccccc}",
  "\\\\[-1.8ex]\\hline",
  "\\hline \\\\[-1.8ex]",
  " & \\multicolumn{6}{c}{\\textit{Dependent variable: log(Hazard) (leaving job)}} \\\\",
  "& \\multicolumn{4}{c}{Manager (1-4)} & \\multicolumn{2}{c}{Asst. Manager (5-6)}\\\\",
  "\\cline{2-5} \\cline{6-7}",
  "& England & ROW & England & ROW & England & ROW\\\\",
  "\\\\[-1.8ex] & (1) & (2) & (3) & (4) & (5) & (6)\\\\",
  "\\hline \\\\[-1.8ex]",
  .finrow("Black",                                    "black"),
  .finrow("UK and other nationality",                 "ukplus"),
  .finrow("Foreign",                                  "foreign"),
  .finrow("-Cumulative Points",                       "points_c"),
  .finrow("-Cumulative Points $\\times$ Black",       "black:points_c"),
  .finrow("Experience at start of job",               "cumulgamesstart", blank = FALSE),
  "\\hline \\\\[-1.8ex]",
  paste0("Observations & ",     paste(nobs_fin,  collapse = " & "), " \\\\ "),
  paste0("Job spells & ",       paste(spell_fin, collapse = " & "), " \\\\ "),
  paste0("Black job exits & ",  paste(evbl_fin,  collapse = " & "), " \\\\ "),
  paste0("Pseudo-R$^{2}$ & ",   paste(rsq_fin,   collapse = " & "), " \\\\ "),
  "\\hline",
  "\\hline \\\\[-1.8ex]",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  paste0("\\item[] \\scriptsize \\textit{Notes:} Cox proportional hazard models for the hazard ",
         "of leaving a job. Coefficients are log hazard ratios; exponentiating gives the hazard ",
         "ratio. Columns 1--4 cover manager spells and columns 5--6 assistant manager spells. ",
         "Within each block England columns are spells at English clubs and ROW columns spells ",
         "at clubs elsewhere. Columns 1 and 2 have no controls, columns 3--6 add cumulative ",
         "points, cumulative points interacted with the Black indicator, managerial experience ",
         "at the start of the spell and nationality dummies. Cumulative points over the last ten ",
         "games is mean-centred within each column and enters with a negative sign, so that a ",
         "higher value denotes worse performance. In columns 3--6 the Black coefficient is ",
         "therefore the racial gap in the exit hazard at average performance, and the ",
         "interaction is the additional gap opened by a one-point shortfall over ten games. ",
         "Nationality dummies have UK-only nationals as the omitted reference category. ",
         "An observation is a game within a job spell. Only spells with at least ten games ",
         "enter the model. Standard errors in parentheses. ",
         "* p$<$0.1; ** p$<$0.05; *** p$<$0.01."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

writeLines(tex_fin, con = file.path(outputwdtables, "table4_hazard_model.tex"))

# APPENDIX TABLE A1: Job descriptions -------------------------------------
funccattable_bn <- data.frame(table(coachdata$func, coachdata$funccat_bn))
funccattable_bn <- spread(funccattable_bn, Var2, Freq)
bn_cols <- c("Manager England", "Manager ROW", "Assistant manager England", "Assistant manager ROW",
             "Youth England", "Youth ROW", "Other staff")
funccattable_bn <- funccattable_bn[, c("Var1", bn_cols)]
funccattable_bn$Total <- rowSums(funccattable_bn[, bn_cols])
funccattable_bn <- funccattable_bn[order(-funccattable_bn$Total), ]
colnames(funccattable_bn)[1] <- "Job Description"

local({
  sink(file.path(outputwdtables, "tableA1_job_descriptions.tex"))
  cat("\\begin{table}[!h]\n")
  cat("\\centering\n")
  cat("\\caption{Job Spell Descriptions and Classifications}\n")
  cat("\\label{tab:jobdescriptions}\n")
  cat("\\footnotesize\n")
  cat("\\setlength\\tabcolsep{3pt}\n")
  cat("\\scalebox{0.9}{\n")
  cat("\\begin{threeparttable}\n")
  cat("\\begin{tabular}{lrrrrrrrrr}\n")
  cat("  \\hline\n")
  cat("Job Description & Mgr Eng & Mgr ROW & Asst Eng & Asst ROW & Youth Eng & Youth ROW & Other & Total \\\\ \n")
  cat("  \\hline\n")
  above20 <- funccattable_bn[funccattable_bn$Total >= 20, ]
  below20 <- funccattable_bn[funccattable_bn$Total <  20, ]
  for (i in seq_len(nrow(above20))) {
    row  <- above20[i, ]
    desc <- as.character(row[["Job Description"]])
    vals <- paste(formatC(as.integer(row[, bn_cols]), width = 4), collapse = " & ")
    cat(desc, " & ", vals, " & ", row$Total, " \\\\ \n")
  }
  cat(". & . & . & . & . & . & . & . & . \\\\ \n")
  cat(". & . & . & . & . & . & . & . & . \\\\ \n")
  cat(". & . & . & . & . & . & . & . & . \\\\ \n")
  last_row  <- below20[nrow(below20), ]
  last_desc <- as.character(last_row[["Job Description"]])
  last_vals <- paste(formatC(as.integer(last_row[, bn_cols]), width = 4), collapse = " & ")
  cat(last_desc, " & ", last_vals, " & ", last_row$Total, " \\\\ \n")
  cat("  \\hline\n")
  totals <- colSums(funccattable_bn[, bn_cols])
  cat("Total & ", paste(formatC(as.integer(totals), width = 4), collapse = " & "), " & ",
      sum(totals), " \\\\ \n")
  cat("  \\hline \n")
  cat("\\end{tabular}\n")
  cat("\\begin{tablenotes}[flushleft]\n")
  cat("\\item[] \\scriptsize \\textit{Notes:} This table categorizes job spells into eight groups based on job type and location (England vs.~rest of world). Each cell indicates the frequency with which a specific job description was classified into a given category. Job descriptions appearing fewer than twenty times are excluded from the table rows but included in the column totals.\n")
  cat("\\end{tablenotes}\n")
  cat("\\end{threeparttable}\n")
  cat("}\n")
  cat("\\end{table}\n")
  sink()
})


# APPENDIX TABLE A2: Probability of getting a coaching degree ---------------------------

regframe_tab2 <- playerlevelstats #playerlevel stats is made in above section for TABLE 1
regframe_tab2$yob <- substr(regframe_tab2$dob, nchar(regframe_tab2$dob)-3, nchar(regframe_tab2$dob))
regframe_tab2$uk <- ifelse(regframe_tab2$nationality %in% ukcountries | regframe_tab2$nationality2 %in% ukcountries | regframe_tab2$nationality3 %in% ukcountries, 1 ,0)

lm1 <- felm(I(licencedum) ~ black, data = regframe_tab2)
summary(lm1)

lm2 <- felm(I(licencedum) ~ black + ukplus + foreign | yob  , data = regframe_tab2)
summary(lm2)

lm3 <- felm(I(licencedum) ~ black + I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 + I(sum_minplayed_tier2/60) + sum_goals_tier2 + sum_assists_tier2 + I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier + attacker + midfielder + defender +  ukplus + foreign + natteam + natapps | yob  , data = regframe_tab2)
summary(lm3)

lms_deg <- list(lm1, lm2, lm3)

local({
  sink(file.path(outputwdtables, "tableA2_coaching_licence.tex"))
  nobs_deg <- sapply(lms_deg, function(m) format(m$N, big.mark = ","))
  r2_deg   <- sapply(lms_deg, function(m) formatC(round(summary(m)$r2, 3), format = "f", digits = 3))
  cat("\\begin{table}[!h]\n")
  cat("\\centering\n")
  cat(" \\caption{Obtaining Coaching Licence by Race}\n")
  cat("   \\label{tab:degree}\n")
  cat("\\footnotesize\n")
  cat("\\scalebox{0.9}{\n")
  cat("\\begin{threeparttable}\n")
  cat("\\begin{tabular}{@{\\extracolsep{5pt}}lccc} \n")
  cat("\\\\[-1.8ex]\\hline \n")
  cat("\\hline \\\\[-1.8ex] \n")
  cat(" & \\multicolumn{3}{c}{\\textit{Dependent variable:}} \\\\ \n")
  cat("\\cline{2-4} \n")
  cat("\\\\[-1.8ex] & \\multicolumn{3}{c}{Obtained any Coaching Licence} \\\\ \n")
  cat("\\\\[-1.8ex] & (1) & (2) & (3)\\\\ \n")
  cat("\\hline \\\\[-1.8ex] \n")
  var_row_lm("Black",                          "black")
  var_row_lm("UK and Other",                   "ukplus")
  var_row_lm("Foreign",                        "foreign")
  var_row_lm("Attacker",                       "attacker")
  var_row_lm("Midfielder",                     "midfielder")
  var_row_lm("Defender",                       "defender")
  var_row_lm("Hours Played Tier 1",            "I(sum_minplayed_tier1/60)")
  var_row_lm("Hours Played Tier 2",            "I(sum_minplayed_tier2/60)")
  var_row_lm("Hours Played Tier $>$ 2",        "I(sum_minplayed_lowtier/60)")
  var_row_lm("Played for National Team",       "natteam")
  var_row_lm("Games Played for National Team", "natapps")
  var_row_lm("Constant",                       "(Intercept)")
  cat("\\hline \\\\[-1.8ex] \n")
  cat("Year of Birth FE & No & Yes & Yes \\\\ \n")
  cat("Observations & ", paste(nobs_deg, collapse = " & "), " \\\\ \n")
  cat("R$^{2}$ & ",      paste(r2_deg,   collapse = " & "), " \\\\ \n")
  cat("\\hline \n")
  cat("\\hline \\\\[-1.8ex] \n")
  cat("\\end{tabular} \n")
  cat("\\begin{tablenotes}[flushleft]\n")
  cat("\\item[] \\scriptsize \\textit{Notes:} The table reports regressions of an indicator for having obtained any coaching licence. Column 1 includes only the Black indicator (raw racial difference). Column 2 adds year-of-birth fixed effects and nationality controls: indicators for UK citizens with a second non-UK citizenship (``UK and Other'') and for foreign (non-UK) citizenship (base category = UK-only citizen). Column 3 further adds playing-position dummies (attacker, midfielder, defender; base = goalkeeper), tiered playing-time measures (hours in Tier 1, Tier 2, and lower tiers), and national-team experience (Played for National Team and Games Played for National Team). Coefficients for goals and assists by tier are included in the regressions but omitted from the table for brevity. Standard errors in parentheses. * p-value $<$ 0.1; ** p-value $<$ 0.05; *** p-value $<$ 0.01.\n")
  cat("\\end{tablenotes}\n")
  cat("\\end{threeparttable}\n")
  cat("}\n")
  cat("\\end{table}\n")
  sink()
})




# APPENDIX TABLE A3 (Logit): entry regressions -------------------------------------------

# APPENDIX TABLE A2 above resets regframe_tab2 to playerlevelstats, which drops the three
# outcome variables built in the TABLE 3 section. They are rebuilt here, identically, so this
# section runs whether the script is executed top to bottom or a chunk at a time.
regframe_tab2$jobinuk <- ifelse(regframe_tab2$Manager_England == 1 |
                                regframe_tab2$Assistant_manager_England == 1 |
                                regframe_tab2$Youth_England == 1, 1, 0)
regframe_tab2$nother  <- ifelse(regframe_tab2$Manager == 1 |
                                regframe_tab2$Assistant_manager == 1 |
                                regframe_tab2$Youth_team_management == 1, 1, 0)
alljobsinuk_links <- unique(coachdata$link[coachdata$country == "England" &
                                           coachdata$funccat != "No Job"])
regframe_tab2$alljobsinuk <- ifelse(regframe_tab2$link %in% alljobsinuk_links, 1, 0)

ctrl_glm <- "+ ukplus + foreign + attacker + midfielder + defender +
               I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
               I(sum_minplayed_tier2/60) + sum_goals_tier2  + sum_assists_tier2 +
               I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
               natteam + natapps + factor(yob)"

.gfit <- function(f) glm(as.formula(f), family = binomial, data = regframe_tab2)
glm1 <- .gfit("I(coached)     ~ black")                       # any job, few controls
glm2 <- .gfit(paste("I(coached)    ~ black", ctrl_glm))       # any job, all controls
glm3 <- .gfit(paste("I(alljobsinuk) ~ black", ctrl_glm))      # any job in England, all controls
glm4 <- .gfit("I(nother)      ~ black")                       # pipeline anywhere, few controls
glm5 <- .gfit(paste("I(nother)     ~ black", ctrl_glm))       # pipeline anywhere, all controls
glm6 <- .gfit(paste("I(jobinuk)    ~ black", ctrl_glm))       # pipeline England, all controls

glms <- list(glm1, glm2, glm3, glm4, glm5, glm6)

n_obs_g <- sapply(glms, function(m) formatC(nobs(m), format = "d", big.mark = ","))
pr2_g   <- sapply(glms, function(m) {
  ll_f <- as.numeric(logLik(m))
  ll_0 <- as.numeric(logLik(update(m, . ~ 1)))
  formatC(round(1 - ll_f / ll_0, 3), format = "f", digits = 3)
})

tex3_logit <- c(
  "\\begin{table}[!h]",
  "\\centering",
  " \\caption{Entry into Executive Labour Market by Race --- Logit}",
  "   \\label{tab:entryregsspecanylogit}",
  "\\footnotesize",
  "\\scalebox{0.85}{",
  "\\begin{threeparttable}",
  "\\begin{tabular}{@{\\extracolsep{15pt}}lcccccc}",
  "\\\\[-1.8ex]\\hline",
  "\\\\[-1.8ex] Dependent variable: & \\multicolumn{3}{c}{Any Job in Football} & \\multicolumn{3}{c}{Pipeline Job} \\\\",
  "\\\\[-1.8ex] & All & All & in England & All & All & in England \\\\",
  "\\\\[-1.8ex] & (1) & (2) & (3) & (4) & (5) & (6)\\\\",
  "\\hline \\\\[-1.8ex]",
  .crow(" Black",        "black",      glms),
  .crow(" UK and Other", "ukplus",     glms),
  .crow(" Foreign",      "foreign",    glms),
  .crow(" Attacker",     "attacker",   glms),
  .crow(" Midfielder",   "midfielder", glms),
  .crow(" Defender",     "defender",   glms, blank = FALSE),
  "  & & & & & & \\\\ ",
  paste0("Average $Y_{i}$ Non-Black & ", paste(corwhiteval, collapse = " & "), " \\\\ "),
  "  & & & & & & \\\\ ",
  "\\hline \\\\[-1.8ex]",
  "Year of Birth FE & No & Yes & Yes & No & Yes & Yes \\\\",
  "Player performance & No & Yes & Yes & No & Yes & Yes \\\\",
  paste0("Observations & ", paste(n_obs_g, collapse = " & "), " \\\\"),
  paste0("McFadden R$^{2}$ & ", paste(pr2_g, collapse = " & "), " \\\\"),
  "\\hline",
  "\\hline \\\\[-1.8ex]",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  paste0("\t\t\\item[] \\scriptsize \\textit{Notes:} This table presents logit regression results. ",
         "Columns (1)--(3) are for obtaining any non-playing job in football (anywhere, without controls in column 1 and with the full set of controls in column 2) and any non-playing job in England (column 3). ",
         "Columns (4)--(6) are for holding a pipeline (managerial) role of manager, assistant manager, or youth-team management: anywhere without controls (4), anywhere with controls (5), and in England (6). ",
         "Where included, controls are fixed effects for year of birth, nationality status (base category = UK-only citizen), playing-position controls (attacker, midfielder, defender; base = goalkeeper), national-team experience and ",
         "national-team appearances, and playing-performance statistics from Tier 1, Tier 2, and lower tiers ",
         "in England (hours played, goals, assists). Coefficients for the player career variables are omitted ",
         "from the output, but available on request. Significance is indicated as follows: ",
         "* p-value$<$ 0.1; ** p-value$<$ 0.05; *** p-value$<$ 0.01."),
  "\t\t\\end{tablenotes}",
  "  \\end{threeparttable}",
  "  }",
  "\\end{table}"
)

writeLines(tex3_logit, con = file.path(outputwdtables, "tableA3_entry_regs_logit.tex"))

# APPENDIX TABLE B1: Balancing table matching -----------------------------

# APPENDIX TABLE A3 above overwrites glm1 with one of its own entry regressions, so the
# propensity model cannot be reused from DATA PREP 2 here: it has to be refitted. Matching on
# the wrong fitted values silently produces an "after matching" sample that is no better
# balanced than the raw one, which is exactly what this table is meant to rule out. Same
# specification as DATA PREP 2, and as APPENDIX FIGURE B1 below.
glm1 <- glm(black ~ I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
              I(sum_minplayed_tier2/60) + sum_goals_tier2 + sum_assists_tier2 +
              I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
              attacker + midfielder + defender + ukonly + ukplus + foreign +
              natteam + natapps + as.factor(yob),
            family = binomial(link = "probit"), data = regframe_tab2)

set.seed(123)
rr1 <- Match(Y = regframe_tab2$coached, Tr = regframe_tab2$black, X = glm1$fitted, M = 1, ties = TRUE, replace = FALSE)
#generating sample
matchedlinks <- unique(regframe_tab2[c(rr1$index.control, rr1$index.treated),])
table(matchedlinks$black)
matchedlinks <- unique(regframe_tab2$link[c(rr1$index.control, rr1$index.treated)])

# Extract the before and after balance data
e <- MatchBalance(black ~ I(sum_minplayed_tier1/60) + sum_goals_tier1 + sum_assists_tier1 + 
                    I(sum_minplayed_tier2/60) + sum_goals_tier2 + sum_assists_tier2 + 
                    I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier + 
                    attacker + midfielder + defender + ukonly + ukplus + foreign + natteam + natapps,
                  match.out = rr1, nboots = 0, data = regframe_tab2)

# Extract means and p-values for before and after matching
before_means_treated <- sapply(e$BeforeMatching, function(x) x$mean.Tr)
before_means_control <- sapply(e$BeforeMatching, function(x) x$mean.Co)
before_pvalues <- sapply(e$BeforeMatching, function(x) x$p.value)

after_means_treated <- sapply(e$AfterMatching, function(x) x$mean.Tr)
after_means_control <- sapply(e$AfterMatching, function(x) x$mean.Co)
after_pvalues <- sapply(e$AfterMatching, function(x) x$p.value)

# Function to add significance stars based on p-value
add_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("***")
  } else if (p < 0.05) {
    return("**")
  } else if (p < 0.1) {
    return("*")
  } else {
    return("")
  }
}

# Calculate the differences and apply significance stars
before_differences <- before_means_treated - before_means_control
before_differences_with_stars <- paste0(round(before_differences, 3), sapply(before_pvalues, add_stars))

after_differences <- after_means_treated - after_means_control
after_differences_with_stars <- paste0(round(after_differences, 3), sapply(after_pvalues, add_stars))

# Covariate names as found in your table
covariate_names <- c(
  "Hours Played Tier 1", "Goals Tier 1", "Assists Tier 1", 
  "Hours Played Tier 2", "Goals Tier 2", "Assists Tier 2", 
  "Hours Played Tier > 2", "Goals Tier > 2", "Assists Tier > 2", 
  "Attacker", "Midfielder", "Defender", 
  "UK Only", "UK and Other", "Foreign","Played for National Team","Games Played for National Team"
)

# Create a dataframe for covariate names, means, differences, and stars
mean_table <- data.frame(
  Covariate = covariate_names,
  Before_Treated = round(before_means_treated, 3),
  Before_Control = round(before_means_control, 3),
  Before_Difference = before_differences_with_stars,
  After_Treated = round(after_means_treated, 3),
  After_Control = round(after_means_control, 3),
  After_Difference = after_differences_with_stars
)

# Output LaTeX table
display_names <- gsub("> 2", "$>$ 2", covariate_names)

sink(file.path(outputwdtables, "tableB1_covariate_balance.tex"))
cat("\\begin{table}[!h]\n")
cat("\\centering\n")
cat("\\caption{Covariate Balance Before and After Matching with Differences} \n")
cat("\\label{tab:balance}\n")
cat("\\footnotesize\n")
cat("%\\hspace*{-0.5cm}\n\n")
cat("\\begin{threeparttable}\n")
cat("\\begin{tabular}{lrrlrrl}\n")
cat("  \\toprule\n")
cat("  & \\multicolumn{3}{c}{Before Matching} & \\multicolumn{3}{c}{After Matching}\\\\\n")
cat("  \\cline{2-4} \\cline{5-7}\n")
cat("Covariate & Black & Non-Black & Difference & Black & Non-Black & Difference \\\\ \n")
cat("  \\midrule\n")
for (i in seq_len(nrow(mean_table))) {
  cat(sprintf("  %s & %s & %s & %s & %s & %s & %s \\\\ \n",
              display_names[i],
              formatC(mean_table$Before_Treated[i], format = "f", digits = 2),
              formatC(mean_table$Before_Control[i], format = "f", digits = 2),
              mean_table$Before_Difference[i],
              formatC(mean_table$After_Treated[i], format = "f", digits = 2),
              formatC(mean_table$After_Control[i], format = "f", digits = 2),
              mean_table$After_Difference[i]
  ))
}
cat("   \\bottomrule\n")
cat("\\end{tabular}\n")
cat("\t\\begin{tablenotes}[flushleft]\n")
cat("\t\t\\item[] \\scriptsize \\textit{Notes:} This table shows the pre and post matching differences in the control variables used throughout this paper. Significance is indicated as follows: * p-value\\(<\\) 0.1; ** p-value\\(<\\) 0.05; *** p-value\\(<\\) 0.01.\n")
cat("\t\t\\end{tablenotes}\n")
cat("  \\end{threeparttable}\n")
cat("  \n")
cat("\\end{table}\n")
sink()


# APPENDIX FIGURE B1: Propensity scorres ----------------------------------

# Propensity score model 
glm1 <- glm(black~ I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 + I(sum_minplayed_tier2/60) + sum_goals_tier2 + sum_assists_tier2 + I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier + attacker + midfielder + defender + ukonly + ukplus + foreign +natteam + natapps + as.factor(yob), family=binomial(link = "probit"), data=regframe_tab2)
summary(glm1)

regframe_tab2$pblack <- glm1$fitted

#plot distribs
plot <- ggplot(data = regframe_tab2, aes(y = pblack, fill = as.factor(black))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_histogram(binwidth = 0.01) +
  # geom_density(alpha = 0.2) +
  coord_flip() +
  ylab("Propensity Score") +
  xlab("Count") +
  scale_fill_manual(name = "Race:", labels = c("Non-Black","Black"), values = cbp5[c(2,5)])
plot

ggsave(file.path(outputwdfigures, "figureB1_propensity_scores.pdf"), device = "pdf", plot, width = 1.5*1.2*127, height = 1.5*1.2*0.7*127, units = "mm")

# Leave the working directory where it started. Each script sets pkgwd from getwd(), so a
# script that ends inside data/ would make the next one in the sequence resolve every path
# one level down and create stray data/ and output/ folders there.
setwd(pkgwd)
