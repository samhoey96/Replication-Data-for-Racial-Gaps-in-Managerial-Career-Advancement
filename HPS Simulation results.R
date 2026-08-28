# Preamble ----------------------------------------------------------------
# HPS Simulation results.R
#
# Every result in the paper that is built from the simulations: the counterfactual prevalence
# figures, the headline decomposition in Table A5, and the multinomial transition model in
# Appendix Table C1.
#
# Run "HPS Simulation preparation.R" first, once for each uksim setting, so that all three
# simulation environments exist in data/. This script loads them; it does not build them.
#
# These sections were part of HPS Racial Gaps.R before being split out. They were always
# self-contained: nothing here is used by the rest of the results code, and nothing here
# depends on it. Keeping them separate also stops the simulation environments from
# overwriting the analysis objects of the same name partway through a results run.

rm(list = ls())

# packages
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
pkgwd  <- getwd()   # or: pkgwd <- "/path/to/Replication Package"
datawd <- file.path(pkgwd, "data")

# The analysis environment supplies the plotting constants and the helper functions. It is a
# save.image(), so loading it also restores the working directories of the machine that built
# it; every path is re-derived below from the current location.
load(file.path(datawd, "analysis_environment.RData"))

pkgwd           <- getwd()
datawd          <- file.path(pkgwd, "data")
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
setwd(datawd)

# In a non-interactive run, printing a plot opens R's default device and leaves an Rplots.pdf
# in the working directory. Every figure here is written explicitly with ggsave, so the default
# device is not needed.
if (!interactive()) pdf(NULL)

# Which simulation environments to read. HPS Simulation preparation.R stamps the replication
# count into the file name, so this picks the run to report on: 1000 is the full one that
# ships with the package, and a smaller number loads a test run if you have built one.
reps <- 10

simenv <- function(sample) file.path(datawd,
  sprintf("data_environment_%s_%d.RData", sample, reps))

# Everything this script writes carries the same replication count, so output from a test run
# sits alongside the full one rather than replacing it. The prefix is untouched, so fig_dir()
# still routes numbered figures to output/figures and the rest to output/other.
with_reps <- function(name) sub("\\.(pdf|tex)$", paste0("_", reps, ".\\1"), name)

for (.f in vapply(c("all", "uk", "ms"), simenv, character(1)))
  if (!file.exists(.f))
    stop("missing simulation environment: ", basename(.f),
         "\nRun HPS Simulation preparation.R with reps = ", reps, " for each uksim setting.")
rm(.f)


# SIMULATION OUTPUT --------------------------------------------

for(uksim in c(0,1,2)){
  
  if (uksim == 0) load(simenv("all"))
  if (uksim == 1) load(simenv("uk"))
  if (uksim == 2) load(simenv("ms"))
  
  print(paste("uksim = ", uksim))
  
  if(uksim == 1){
    ticks <- seq(0, 0.04, 0.005)
  }
  if(uksim == 0){
    ticks <- seq(0, 0.25, 0.005)
  }
  plotframe_share <- subset.data.frame(plotframe_share, days_elapsed %in% seq(15, 3650, 14))
  
  #labels
  types <- c("Actual (A)", "Simulated actual (SA)", "First entry as other race (C1)", "Progression as other race (C2)", "First entry and progression as other race (C3)", "Enter and progress as other race\nInitial state modelled")
  plotframe_share$type <- "1 Actual"
  plotframe_share_simact$type <- "2 Simulated actual"
  plotframe_share_simeao$type <- "3 Enter as other race"
  plotframe_share_simpao$type <- "4 Progress as other race"
  plotframe_share_simom$type <- "5 Enter and progress as other race"
  
  # Simulation plots prevalence simulated vs actual -------------------------------------------------------------
  
  
  errorbounds <- do.call(rbind, statelist_simact)
  result <- errorbounds %>%
    group_by(days_elapsed, black, funccat_bn) %>%
    summarize(
      p25 = quantile(share, probs = 0.025, na.rm = TRUE),
      p975 = quantile(share, probs = 0.975, na.rm = TRUE),
    ) %>%
    ungroup()
  
  plotframe_share_simact_alt <- join(plotframe_share_simact, result)
  plotframe_share$p25 <- plotframe_share$share
  plotframe_share$p975 <- plotframe_share$share
  
  plotframe_share_sim_plot <- rbind(plotframe_share, plotframe_share_simact_alt)
  
  #manager
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, funccat_bn == "Manager England")
  
  # Create the bar chart using ggplot
  barplot <- ggplot(plotframe_share_sim_plot, aes(x = days_elapsed, y = share, color = as.factor(type))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=12),
          legend.position = "bottom",
          axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize), axis.title.x = element_text(size=axislabsize),
          axis.title.y = element_text(size=axislabsize), plot.title = element_text(hjust = 0.5, size = textsize), axis.ticks.length =unit(0.2, "cm"),
          legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
    geom_line(linewidth = 1.5) +
    geom_errorbar(aes(ymin = p25, ymax = p975), width = 0, alpha = 0.5, color = "grey") +
    facet_grid(~ black, labeller = labeller(black = c("1" = paste("Black Former Players") , "0" = paste("Non-Black Former Players")))) +
    labs(x = "Days since end of playing career", y = "Manager Prevalence", color = "Simulation:") +
    scale_x_continuous(name = "Years since end of playing career", breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = "Manager Prevalence", breaks = ticks, labels = ticks) +
    scale_fill_manual(values = cbp5) +
    geom_vline(xintercept = c(1, seq(0, 3650, 365)), linetype = "dashed", color = "white", linewidth = 0.2) +
    guides(color = guide_legend(nrow = 1)) +
    scale_color_manual(name = "Simulation:", values = cbp6, labels = types)
  barplot
  
  if(uksim==0){
    savename <- "figure7_manager_sim_vs_actual_eng.pdf"
  }
  if(uksim==1){
    savename <- "figureC1_manager_sim_vs_actual_eng_uk.pdf"
  }
  if(uksim==2){
    savename <- "figureB3_manager_sim_vs_actual_eng_ms.pdf"
  }
  ggsave(file.path(fig_dir(savename), with_reps(savename)), device = "pdf", barplot,
         width = 1.5*1.2*127, height = 1.5*1.2*0.7*127, units = "mm")
  
  # manager ROW
  plotframe_share_sim_plot <- rbind(plotframe_share, plotframe_share_simact_alt)
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, funccat_bn == "Manager ROW")
  
  barplot <- ggplot(plotframe_share_sim_plot, aes(x = days_elapsed, y = share, color = as.factor(type))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=12),
          legend.position = "bottom",
          axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize), axis.title.x = element_text(size=axislabsize),
          axis.title.y = element_text(size=axislabsize), plot.title = element_text(hjust = 0.5, size = textsize), axis.ticks.length =unit(0.2, "cm"),
          legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
    geom_line(linewidth = 1.5) +
    geom_errorbar(aes(ymin = p25, ymax = p975), width = 0, alpha = 0.5, color = "grey") +
    facet_grid(~ black, labeller = labeller(black = c("1" = paste("Black Former Players") , "0" = paste("Non-Black Former Players")))) +
    labs(x = "Days since end of playing career", y = "Manager ROW Prevalence", color = "Simulation:") +
    scale_x_continuous(name = "Years since end of playing career", breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = "Manager ROW Prevalence", breaks = ticks, labels = ticks) +
    scale_fill_manual(values = cbp5) +
    geom_vline(xintercept = c(1, seq(0, 3650, 365)), linetype = "dashed", color = "white", linewidth = 0.2) +
    guides(color = guide_legend(nrow = 1)) +
    scale_color_manual(name = "Simulation:", values = cbp6, labels = types)
  barplot
  
  if(uksim==0){
    savename <- "figure7_manager_sim_vs_actual_row.pdf"
  }
  if(uksim==1){
    savename <- "figureC1_manager_sim_vs_actual_row_uk.pdf"
  }
  if(uksim==2){
    savename <- "figureB3_manager_sim_vs_actual_row_ms.pdf"
  }
  ggsave(file.path(fig_dir(savename), with_reps(savename)), device = "pdf", barplot,
         width = 1.5*1.2*127, height = 1.5*1.2*0.7*127, units = "mm")
  
  #any job
  plotframe_share_sim_plot <- rbind(plotframe_share, plotframe_share_simact_alt)
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, funccat_bn == "No Job")
  plotframe_share_sim_plot$share <- 1 - plotframe_share_sim_plot$share
  plotframe_share_sim_plot$p25 <- 1 - plotframe_share_sim_plot$p25
  plotframe_share_sim_plot$p975 <- 1 - plotframe_share_sim_plot$p975
  
  # Create the bar chart using ggplot
  barplot <- ggplot(plotframe_share_sim_plot, aes(x = days_elapsed, y = share, color = as.factor(type))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=12),
          legend.position = "bottom",
          axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize), axis.title.x = element_text(size=axislabsize),
          axis.title.y = element_text(size=axislabsize), plot.title = element_text(hjust = 0.5, size = textsize), axis.ticks.length =unit(0.2, "cm"),
          legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
    geom_line(linewidth = 1.5) +
    geom_errorbar(aes(ymin = p25, ymax = p975), width = 0, alpha = 0.5, color = "grey") +
    facet_grid(~ black, labeller = labeller(black = c("1" = paste("Black Former Players") , "0" = paste("Non-Black Former Players")))) +
    labs(x = "Days since end of playing career", y = "Manager Prevalence", color = "Simulation:") +
    scale_x_continuous(name = "Years since end of playing career", breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = "Employment Prevalence", breaks = seq(0, 0.25, 0.05), labels = seq(0, 0.25, 0.05), limits = c(0,0.25)) +
    scale_fill_manual(values = cbp5) +
    geom_vline(xintercept = c(1, seq(0, 3650, 365)), linetype = "dashed", color = "white", linewidth = 0.2) +
    guides(color = guide_legend(nrow = 1)) +
    scale_color_manual(name = "Simulation:", values = cbp6, labels = types)
  barplot 
  
  if(uksim==0){
    savename <- "figureA2_employment_sim_vs_actual.pdf"
  }
  if(uksim==1){
    savename <- "employment_sim_vs_actual_uk.pdf"
  }
  if(uksim==2){
    savename <- "employment_sim_vs_actual_ms.pdf"
  }
  ggsave(file.path(fig_dir(savename), with_reps(savename)), device = "pdf", barplot,
         width = 1.5*1.2*127, height = 1.5*1.2*0.7*127, units = "mm")
  
  # Simulation plots prevalence manager-------------------------------------------------------------
  
  
  
  plotframe_share <- subset.data.frame(plotframe_share, days_elapsed %in% seq(15, 3650, 14))
  
  
  #labels
  types <- c("Simulated actual Non-Black", "Simulated actual Black", "Black Enter as Non-Black (C1)", "Black Progress as Non-Black (C2)", "Black Enter and Progress as Non-Black (C3)")
  plotframe_share$type <- "1 Actual"
  plotframe_share_simact$type <- "2 Simulated actual Black"
  plotframe_share_simeao$type <- "3 Black Enter as Non-Black (C1)"
  plotframe_share_simpao$type <- "4 Black Progress as Non-Black (C2)"
  plotframe_share_simom$type <- "5 Black Enter and Progress as Non-Black (C3)"
  # plotframe_share_simoma$type <- "6 Enter and progress as other race\nInitial state modelled"
  
  plotframe_share_sim_plot <- rbind(plotframe_share, plotframe_share_simact)
  
  
  #manager
  plotframe_share_sim_plot <- rbind(plotframe_share_simact , plotframe_share_simeao)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simpao)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simom)
  #plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simoma)
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, funccat_bn %in% c("Manager England", "Manager ROW"))
  
  #adjustments
  whitesa <- subset.data.frame(plotframe_share_sim_plot, black == 0 & type == "2 Simulated actual Black")
  whitesa$type <- "1 Simulated Actual Non-Black"
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, black == 1)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot, whitesa)
  
  # Create the bar chart using ggplot
  barplot <- ggplot(plotframe_share_sim_plot, aes(x = days_elapsed, y = share, color = as.factor(type))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=12),
          legend.position = "bottom", 
          axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize), axis.title.x = element_text(size=axislabsize),
          axis.title.y = element_text(size=axislabsize), plot.title = element_text(hjust = 0.5, size = textsize), axis.ticks.length =unit(0.2, "cm"), 
          legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
    geom_line(linewidth = 1.5) +
    facet_grid(~ funccat_bn, labeller = labeller(black = c("1" = paste("Black Former Players") , "0" = paste("Non-Black Former Players")))) +
    labs(x = "Days since end of playing career", y = "Manager Prevalence", color = "Simulation:") +
    scale_x_continuous(name = "Years since end of playing career", breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = "Manager Prevalence", breaks = ticks, labels = ticks) +
    scale_fill_manual(values = cbp5) +
    geom_vline(xintercept = c(1, seq(0, 3650, 365)), linetype = "dashed", color = "white", linewidth = 0.2) +
    guides(color = guide_legend(nrow = 3)) +
    scale_color_manual(name = "", values = cbp6[c(6, 2,3, 4,5)], labels = types) +
    expand_limits(x = 0, y = 0)
  barplot 
  
  
  if(uksim==0){
    savename <- "figure8_manager_sim.pdf"
  }
  if(uksim==1){
    savename <- "figureC2_manager_sim_uk.pdf"
  }
  if(uksim==2){
    savename <- "figureB4_manager_sim_ms.pdf"
  }
  ggsave(file.path(fig_dir(savename), with_reps(savename)), device = "pdf", barplot,
         width = 1.5*1.2*127, height = 1.5*1.2*0.7*127, units = "mm")
  
  
  # Simulation plots prevalence assistant manager-------------------------------------------------------------
  
  #assistant manager
  plotframe_share_sim_plot <- rbind(plotframe_share_simact , plotframe_share_simeao)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simpao)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simom)
  #plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simoma)
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, funccat_bn %in% c("Assistant manager England", "Assistant manager ROW"))
  
  #adjustments
  whitesa <- subset.data.frame(plotframe_share_sim_plot, black == 0 & type == "2 Simulated actual Black")
  whitesa$type <- "1 Simulated Actual Non-Black"
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, black == 1)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot, whitesa)
  
  
  # Create the bar chart using ggplot
  barplot <- ggplot(plotframe_share_sim_plot, aes(x = days_elapsed, y = share, color = as.factor(type))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=12),
          legend.position = "bottom", 
          axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize), axis.title.x = element_text(size=axislabsize),
          axis.title.y = element_text(size=axislabsize), plot.title = element_text(hjust = 0.5, size = textsize), axis.ticks.length =unit(0.2, "cm"), 
          legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
    geom_line(linewidth = 1.5) +
    facet_grid(~ funccat_bn, labeller = labeller(black = c("1" = paste("Black Former Players") , "0" = paste("Non-Black Former Players")))) +
    labs(x = "Days since end of playing career", y = "Assistant Manager Prevalence", color = "Simulation:") +
    scale_x_continuous(name = "Years since end of playing career", breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = "Manager Prevalence", breaks = ticks, labels = ticks) +
    scale_fill_manual(values = cbp5) +
    geom_vline(xintercept = c(1, seq(0, 3650, 365)), linetype = "dashed", color = "white", linewidth = 0.2) +
    guides(color = guide_legend(nrow = 3)) +
    scale_color_manual(name = "", values = cbp6[c(6, 2,3, 4,5)], labels = types) +
    expand_limits(x = 0, y = 0)
  barplot 
  
  if(uksim==0){
    savename <- "figure9_assman_sim.pdf"
  }
  if(uksim==1){
    savename <- "figureC3_assman_sim_uk.pdf"  
  }
  if(uksim==2){
    savename <- "figureB5_assman_sim_ms.pdf"  
  }
  ggsave(file.path(fig_dir(savename), with_reps(savename)), device = "pdf", barplot,
         width = 1.5*1.2*127, height = 1.5*1.2*0.7*127, units = "mm")
  
  # Simulation plots prevalence youth -------------------------------------------------------------
  
  #youth
  plotframe_share_sim_plot <- rbind(plotframe_share_simact , plotframe_share_simeao)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simpao)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simom)
  #plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot , plotframe_share_simoma)
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, funccat_bn %in% c("Youth ROW", "Youth England"))
  
  #adjustments
  whitesa <- subset.data.frame(plotframe_share_sim_plot, black == 0 & type == "2 Simulated actual Black")
  whitesa$type <- "1 Simulated Actual Non-Black"
  plotframe_share_sim_plot <- subset.data.frame(plotframe_share_sim_plot, black == 1)
  plotframe_share_sim_plot <- rbind(plotframe_share_sim_plot, whitesa)
  
  
  # Create the bar chart using ggplot
  barplot <- ggplot(plotframe_share_sim_plot, aes(x = days_elapsed, y = share, color = as.factor(type))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=12),
          legend.position = "bottom", 
          axis.text.y = element_text(size = ticksize), axis.text.x = element_text(size = ticksize), axis.title.x = element_text(size=axislabsize),
          axis.title.y = element_text(size=axislabsize), plot.title = element_text(hjust = 0.5, size = textsize), axis.ticks.length =unit(0.2, "cm"), 
          legend.text = element_text(size = legendtextsize), strip.text = element_text(size = facetsize)) +
    geom_line(linewidth = 1.5) +
    facet_grid(~ funccat_bn, labeller = labeller(black = c("1" = paste("Black Former Players") , "0" = paste("Non-Black Former Players")))) +
    labs(x = "Days since end of playing career", y = "Youth Manager Prevalence", color = "Simulation:") +
    scale_x_continuous(name = "Years since end of playing career", breaks = seq(0, 3650, 365), labels = seq(0, 10, 1)) +
    scale_y_continuous(name = "Manager Prevalence", breaks = ticks, labels = ticks) +
    scale_fill_manual(values = cbp5) +
    geom_vline(xintercept = c(1, seq(0, 3650, 365)), linetype = "dashed", color = "white", linewidth = 0.2) +
    guides(color = guide_legend(nrow = 3)) +
    scale_color_manual(name = "", values = cbp6[c(6, 2,3, 4,5)], labels = types) +
    expand_limits(x = 0, y = 0)
  barplot 
  
  
  if(uksim==0){
    savename <- "figureA3_youth_sim.pdf"
  }
  if(uksim==1){
    savename <- "youth_sim_uk.pdf"  
  }
  if(uksim==2){
    savename <- "youth_sim_ms.pdf"  
  }
  ggsave(file.path(fig_dir(savename), with_reps(savename)), device = "pdf", barplot,
         width = 1.5*1.2*127, height = 1.5*1.2*0.7*127, units = "mm")
  
  
  # End of simulation loop --------------------------------------------------
}

# TABLE A5: Headline decomposition of manager-prevalence gap ----------------

load(simenv("all"))

decomp_years <- c(3, 6, 10)
decomp_days <- sapply(decomp_years * 365, function(x) {
  plotframe_share_simact$days_elapsed[which.min(abs(plotframe_share_simact$days_elapsed - x))]
})

decomp_sources <- list(
  black_sa    = plotframe_share_simact,
  nonblack_sa = plotframe_share_simact,
  c1          = plotframe_share_simeao,
  c2          = plotframe_share_simpao,
  c3          = plotframe_share_simom
)

decomp_get <- function(df, state, day, black_value) {
  out <- df$share[df$funccat_bn == state & df$days_elapsed == day & df$black == black_value]
  if (length(out) == 0) return(NA_real_)
  out[1]
}

decomp_fmt_pct <- function(x) {
  if (is.na(x)) return("")
  formatC(100 * x, format = "f", digits = 1)
}

decomp_fmt_closed <- function(x) {
  if (is.na(x)) return("")
  paste0(formatC(100 * x, format = "f", digits = 0), "\\%")
}

decomp_state_label <- function(x) {
  dplyr::recode(x, "Manager England" = "Manager England", "Manager ROW" = "Manager ROW")
}

decomp_rows <- do.call(rbind, lapply(c("Manager England", "Manager ROW"), function(state) {
  do.call(rbind, lapply(seq_along(decomp_years), function(i) {
    day <- decomp_days[i]
    black_sa <- decomp_get(decomp_sources$black_sa, state, day, 1)
    nonblack_sa <- decomp_get(decomp_sources$nonblack_sa, state, day, 0)
    c1 <- decomp_get(decomp_sources$c1, state, day, 1)
    c2 <- decomp_get(decomp_sources$c2, state, day, 1)
    c3 <- decomp_get(decomp_sources$c3, state, day, 1)
    gap <- nonblack_sa - black_sa
    data.frame(
      state = state,
      year = decomp_years[i],
      day = day,
      black_sa = black_sa,
      nonblack_sa = nonblack_sa,
      c1 = c1,
      c2 = c2,
      c3 = c3,
      c1_closed = (c1 - black_sa) / gap,
      c2_closed = (c2 - black_sa) / gap,
      c3_closed = (c3 - black_sa) / gap,
      stringsAsFactors = FALSE
    )
  }))
}))

sink(file.path(outputwdtables, with_reps("tableA5_headline_decomposition.tex")))
cat("\\begin{table}[!h]\n")
cat("\\centering\n")
cat("\\caption{Counterfactual Decomposition of Manager Prevalence}\n")
cat("\\label{tab:headline_decomposition}\n")
cat("\\footnotesize\n")
cat("\\begin{threeparttable}\n")
cat("\\begin{tabular}{llrrrrrrrr}\n")
cat("\\toprule\n")
cat("State & Year & Black SA & NB SA & C1 & C2 & C3 & Gap C1 & Gap C2 & Gap C3 \\\\ \n")
cat("\\midrule\n")
for (i in seq_len(nrow(decomp_rows))) {
  row <- decomp_rows[i, ]
  cat(sprintf(
    "%s & %d & %s & %s & %s & %s & %s & %s & %s & %s \\\\ \n",
    decomp_state_label(row$state),
    row$year,
    decomp_fmt_pct(row$black_sa),
    decomp_fmt_pct(row$nonblack_sa),
    decomp_fmt_pct(row$c1),
    decomp_fmt_pct(row$c2),
    decomp_fmt_pct(row$c3),
    decomp_fmt_closed(row$c1_closed),
    decomp_fmt_closed(row$c2_closed),
    decomp_fmt_closed(row$c3_closed)
  ))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}[flushleft]\n")
cat("\\item[] \\scriptsize \\textit{Notes:} Entries report simulated manager prevalence, expressed as percentages of former players, at the closest available 14-day interval to years 3, 6, and 10 after the end of the playing career. SA denotes simulated actual prevalence and NB denotes non-Black. C1 sets the entry process for Black former players to the non-Black process until first observed employment; C2 sets post-entry progression to the non-Black process; C3 sets both entry and progression to the non-Black process. Gap closed is calculated as $(\\mathrm{C}k - \\mathrm{Black\\ SA})/(\\mathrm{NB\\ SA} - \\mathrm{Black\\ SA})$ and may exceed 100\\% when a counterfactual overshoots non-Black simulated prevalence. The table reports point estimates from the existing simulation output and does not include parameter-uncertainty intervals.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()

# APPENDIX TABLE C1: Simplified multinomial logit (two-panel: non-Black baseline + Black effects) -----
load(simenv("all"))

jobstates_all <- c("Manager England", "Manager ROW", "Assistant manager England",
                   "Assistant manager ROW", "Youth England", "Youth ROW", "Other staff", "No Job")

dss <- date_sequence_sim
dss$funccat_bn     <- relevel(factor(dss$funccat_bn,     levels = jobstates_all), ref = "No Job")
dss$funccat_bn_lag <- relevel(factor(dss$funccat_bn_lag, levels = jobstates_all), ref = "No Job")
dss <- dss[!is.na(dss$funccat_bn) & !is.na(dss$funccat_bn_lag), ]

sim_simple <- multinom(
  funccat_bn ~
    as.factor(black) * funccat_bn_lag +
    days_elapsed +
    I(sum_minplayed_tier1/60) + sum_goals_tier1  + sum_assists_tier1 +
    I(sum_minplayed_tier2/60) + sum_goals_tier2  + sum_assists_tier2 +
    I(sum_minplayed_lowtier/60) + sum_goals_lowtier + sum_assists_lowtier +
    attacker + midfielder + defender + as.numeric(yob) + ukplus + foreign + natteam + natapps,
  data  = dss,
  decay = 1e-3,
  maxit = 1000
)

s  <- summary(sim_simple)
cf <- s$coefficients
se <- s$standard.errors
pv <- 2 * (1 - pnorm(abs(cf / se)))

# All 7 estimated outcomes (rows of cf); No Job is the reference outcome
outcomes <- rownames(cf)

# 7 estimated (non-reference) outcomes only; No Job is the reference and not shown
display_outcomes <- outcomes

out_labels <- c(
  "Manager England"             = "Mgr Eng",
  "Manager ROW"                 = "Mgr ROW",
  "Assistant manager England"   = "Asst Mgr Eng",
  "Assistant manager ROW"       = "Asst Mgr ROW",
  "Youth England"               = "Youth Eng",
  "Youth ROW"                   = "Youth ROW",
  "Other staff"                 = "Other Staff"
)

lag_non_ref <- levels(dss$funccat_bn_lag)[levels(dss$funccat_bn_lag) != "No Job"]

# Panel A: non-Black previous-state (funccat_bn_lag) coefficients
lag_coef_names <- paste0("funccat_bn_lag", lag_non_ref)
lag_row_labs   <- c(
  "Mgr England",
  "Mgr ROW",
  "Asst Mgr England",
  "Asst Mgr ROW",
  "Youth England",
  "Youth ROW",
  "Other Staff"
)

# Panel B: Black main effect and Black x previous-state interactions
black_main     <- "as.factor(black)1"
black_int      <- paste0("as.factor(black)1:funccat_bn_lag", lag_non_ref)
all_black      <- c(black_main, black_int)
black_row_labs <- c(
  "Black",
  "\\quad $\\times$ Mgr England",
  "\\quad $\\times$ Mgr ROW",
  "\\quad $\\times$ Asst Mgr Eng",
  "\\quad $\\times$ Asst Mgr ROW",
  "\\quad $\\times$ Youth Eng",
  "\\quad $\\times$ Youth ROW",
  "\\quad $\\times$ Other Staff"
)

stars_sim <- function(p) ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$", ifelse(p < 0.10, "$^{*}$", "")))
fmt_b_sim <- function(x) {
  s <- formatC(round(abs(x), 3), format = "f", digits = 3)
  if (x < 0) paste0("$-$", s) else s
}

# loops over display_outcomes (7 non-reference states; No Job reference column omitted)
make_rows_sim <- function(col_name, label) {
  coef_cells <- sapply(display_outcomes, function(o) {
    if (!col_name %in% colnames(cf)) return("---")
    paste0(fmt_b_sim(cf[o, col_name]), stars_sim(pv[o, col_name]))
  })
  se_cells <- sapply(display_outcomes, function(o) {
    if (!col_name %in% colnames(cf)) return("")
    paste0("(", formatC(round(se[o, col_name], 3), format = "f", digits = 3), ")")
  })
  paste0(
    label, " & ", paste(coef_cells, collapse = " & "), " \\\\\n",
    "  & ", paste(se_cells, collapse = " & "), " \\\\\n",
    " & & & & & & & \\\\\n"
  )
}

nobs_sim <- format(nrow(dss), big.mark = ",")

sink(file.path(outputwdtables, with_reps("tableD1_simplified_transition_model.tex")))
cat("\\begin{table}[!h]\n")
cat("\\centering\n")
cat("\\caption{Simplified Career Transition Model}\n")
cat("\\label{tab:sim_multinom}\n")
cat("\\footnotesize\n")
cat("\\scalebox{0.72}{\n")
cat("\\begin{threeparttable}\n")
cat("\\begin{tabular}{lrrrrrrr}\n")
cat("\\\\[-1.8ex]\\hline\n")
cat("\\hline \\\\[-1.8ex]\n")
col_hdr_sim <- paste(sapply(display_outcomes, function(o) out_labels[o]), collapse = " & ")
cat("Destination: &", col_hdr_sim, "\\\\\n")
cat("\\cline{2-8}\n")
cat("\\\\[-1.8ex]\n")
# Panel A
cat("\\multicolumn{8}{l}{\\textit{Panel A: Non-Black baseline transitions (ref.\\ previous state: No Job)}} \\\\\n")
cat("\\\\[-1.8ex]\n")
for (i in seq_along(lag_coef_names)) cat(make_rows_sim(lag_coef_names[i], lag_row_labs[i]))
cat("\\hline \\\\[-1.8ex]\n")
# Panel B
cat("\\multicolumn{8}{l}{\\textit{Panel B: Black effect and interactions (ref.\\ previous state: No Job)}} \\\\\n")
cat("\\\\[-1.8ex]\n")
for (i in seq_along(all_black)) cat(make_rows_sim(all_black[i], black_row_labs[i]))
cat("\\hline \\\\[-1.8ex]\n")
cat("Observations &", paste(rep(nobs_sim, 7), collapse = " & "), "\\\\\n")
cat("\\hline\n")
cat("\\hline \\\\[-1.8ex]\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}[flushleft]\n")
cat("\\item[] \\scriptsize \\textit{Notes:} The table reports coefficients from a simplified multinomial logit career transition model. The dependent variable is the job state at the next 14-day interval; the reference outcome is No Job (shown as ``---''). Panel A shows the previous job state coefficients for non-Black players (reference previous state: No Job). Panel B shows the Black indicator and Black $\\times$ previous-state interaction terms, capturing how Black players' transitions differ from non-Black players. The model additionally controls for days elapsed since end of playing career (linear), hours played by tier, goals, assists, playing position, year of birth, and nationality. The full simulation model additionally interacts Black and the previous job state with days elapsed. Standard errors in parentheses. $^{*}$ $p<0.10$, $^{**}$ $p<0.05$, $^{***}$ $p<0.01$.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("}\n")
cat("\\end{table}\n")
sink()

# Leave the working directory where it started. Each script sets pkgwd from getwd(), so a
# script that ends inside data/ would make the next one in the sequence resolve every path
# one level down and create stray data/ and output/ folders there.
setwd(pkgwd)
