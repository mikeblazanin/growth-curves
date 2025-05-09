#Setup and libraries ----
library(gcplyr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(forcats)

mysplit <- strsplit(getwd(), split = "/")[[1]]
if(mysplit[length(mysplit)] != "Empirical") {setwd("./Empirical/")}

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Read data ----
trav_resis <- read.csv("trav-phage_resis_data.csv")
trav_resis_sum <-
  summarize(group_by(trav_resis, Proj, Pop, Treat, Timepoint, Isol),
            .groups = "drop",
            mean_eop = mean(EOP),
            bd = any(bd),
            bacteria = paste0(Proj[1], Pop[1], Treat[1], Isol[1]))
#Fix ancestor for easy merging
trav_resis_sum <- trav_resis_sum[-which(trav_resis_sum$Pop == "Anc")[-1], ]
trav_resis_sum <- 
  mutate(trav_resis_sum,
         bacteria = ifelse(Pop == "Anc", "Anc", bacteria),
         Proj = fct_recode(Proj, "Weak" = "7x", "Strong" = "125"),
         Treat = 
           fct_recode(Treat, "Control" = "C", "Local" = "L", "Global" = "G"),
         plot_names = ifelse(Pop == "Anc", "Ancestor\nEOP = 1",
                             paste0(Proj, " ", Treat, "\nPop ", Pop, 
                                    " T", Timepoint, " Isol ", Isol,
                                    "\nEOP = ", signif(mean_eop, 1))))

#Main data experiments:
#2021-11-03 Testing different isolates with initial 1e5 bact, 1e4 phage
#2022-06-07 Testing different isolates with initial 1e5 bact, 1e4 phage
#2022-07-15 Testing different isolates with initial 1e5 bact, 1e4 phage
#    Note that the wells of PF with no phage and 7xCLD with no phage were 
#    flipped in the 2022-07-15 run

gcdata <-
  read_wides(
    files = c("2021-11-03_Emma_Growth_Curve.csv",
              "2022-06-07_William_Growth_Curve.csv",
              "2022-07-15_William_Growth_Curve.csv"),
    startrow = c(29, 41, 41), startcol = c("B", "B", "B"))

#Drop temperature columns
gcdata <- lapply(X = gcdata,
                 FUN = function(x) x[, -grep("^T.*600$", colnames(x))])
              
gcdata_lng <- trans_wide_to_tidy(gcdata,
                                id_cols = c("file", "Time"),
                                values_to = "OD600")

#Create designs ----
design_isols <- 
  make_design(output_format = "tidy",
    nrows = 8, ncols = 12,
    bacteria = make_designpattern(c("PF", "125ALE", "7xEGC", "7xACD",
                                    "125ALA", "125CGE"),
                                  rows = 2:7, cols = 2:4,
                                  pattern = "111222333444555666",
                                  byrow = TRUE),
    bacteria = make_designpattern(c("7xCLD", "7xEGD",
                                    "7xALA", "125CGA", "Blank"),
                                  rows = 2:6, cols = 5:7,
                                  pattern = "111222333444555",
                                  byrow = TRUE),
    bacteria = make_designpattern(c("PF", "125ALE", "7xEGC", "7xACD",
                                    "125ALA", "125CGE"),
                                  rows = 2:7, cols = 8:9,
                                  pattern = "112233445566",
                                  byrow = TRUE),
    bacteria = make_designpattern(c("7xCLD", "7xEGD",
                                    "7xALA", "125CGA", "Blank"),
                                  rows = 2:6, cols = 10:11,
                                  pattern = "1122334455",
                                  byrow = TRUE)
  )

design_isols_2022_07_15 <- 
  make_design(output_format = "tidy",
              nrows = 8, ncols = 12,
              bacteria = make_designpattern(c("PF", "125ALE", "7xEGC", "7xACD",
                                              "125ALA", "125CGE"),
                                            rows = 2:7, cols = 2:4,
                                            pattern = "111222333444555666",
                                            byrow = TRUE),
              bacteria = make_designpattern(c("7xCLD", "7xEGD",
                                              "7xALA", "125CGA", "Blank"),
                                            rows = 2:6, cols = 5:7,
                                            pattern = "111222333444555",
                                            byrow = TRUE),
              bacteria = make_designpattern(c("7xCLD", "125ALE", "7xEGC", "7xACD",
                                              "125ALA", "125CGE"),
                                            rows = 2:7, cols = 8:9,
                                            pattern = "112233445566",
                                            byrow = TRUE),
              bacteria = make_designpattern(c("PF", "7xEGD",
                                              "7xALA", "125CGA", "Blank"),
                                            rows = 2:6, cols = 10:11,
                                            pattern = "1122334455",
                                            byrow = TRUE)
  )

design_phage <- make_design(output_format = "tidy",
                            nrows = 8, ncol = 12,
                            phage_added = mdp("Yes", rows = 2:7, cols = 2:7),
                            phage_added = mdp("No", rows = 2:7, cols = 8:11))
                            
#Merge design, measures, EOP data ----
gcdata_lng[[1]] <- merge_dfs(gcdata_lng[[1]], design_isols)
gcdata_lng[[2]] <- merge_dfs(gcdata_lng[[2]], design_isols)
gcdata_lng[[3]] <- merge_dfs(gcdata_lng[[3]], design_isols_2022_07_15)

gcdata_lng <- lapply(X = gcdata_lng, 
                     FUN = function(x, y) merge_dfs(x = x, y = y, drop = TRUE),
                     y = design_phage)

gcdata_tidy <- merge_dfs(gcdata_lng, collapse = TRUE, drop = TRUE)

#Rename PF to Anc for eventual merging
gcdata_tidy <- mutate(gcdata_tidy,
                     bacteria = ifelse(bacteria == "PF", "Anc", bacteria))

#Convert Time, subtract blank, smooth, calc_deriv, summarize
gcdata_tidy <- mutate(group_by(gcdata_tidy, file),
                      Time = time_length(hms(Time), unit = "hour"),
                      OD600 = OD600 - min(OD600[bacteria == "Blank"]))
gcdata_tidy <- mutate(group_by(gcdata_tidy, file, Well),
                      smoothed = smooth_data(x = Time, y = OD600,
                                             sm_method = "moving-average",
                                             window_width = 1.3),
                      deriv = calc_deriv(y = smoothed, x = Time,
                                         window_width_n = 3))

#Join EOP data
gcdata_tidy <- left_join(gcdata_tidy, trav_resis_sum)
gcdata_tidy$plot_names <- as.factor(gcdata_tidy$plot_names)
gcdata_tidy$plot_names <-
  fct_reorder(gcdata_tidy$plot_names, .x = gcdata_tidy$mean_eop,
              .desc = TRUE, na.rm = FALSE)

gcdata_sum <- dplyr::summarize(
  group_by(gcdata_tidy, 
           file, Well, bacteria, phage_added, Proj, Pop, Treat, Timepoint, Isol),
  peak_dens = first_maxima(y = smoothed, x = Time, return = "y"),
  peak_time = first_maxima(y = smoothed, x = Time, return = "x"),
  auc = auc(x = Time, y = smoothed),
  death_slope = min_gc(deriv),
  eop = mean(mean_eop),
  bd = any(bd))

gcdata_sum_sum <- summarize(
  group_by(gcdata_sum, 
           file, bacteria, phage_added, Proj, Pop, Treat, Timepoint, Isol),
  across(c(peak_dens, peak_time, auc, death_slope, eop), 
         list(mean = mean, sd = sd)))

gcdata_sum_sum_sum <- summarize(
  group_by(gcdata_sum_sum, 
           bacteria, phage_added, Proj, Pop, Treat, Timepoint, Isol),
  across(c(peak_dens_mean, peak_time_mean, auc_mean, death_slope_mean, eop_mean), 
         list(mean = mean, min = min, max = max)))

#Plots ----
#Example plot
png("./example_plot.png",
    width = 6, height = 4, units = "in", res = 150)
ggplot(data = filter(gcdata_tidy, file == "2021-11-03_Emma_Growth_Curve",
                     Well %in% c("B2", "B8"), Time <= 12),
       aes(x = Time, y = smoothed, color = phage_added)) +
  geom_point(size = 2) +
  scale_color_manual(values = my_cols[c(8, 3)], name = "Phage Added?") +
  scale_x_continuous(breaks = c(0, 6, 12)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.spacing.y = unit(.1, "in"),
        plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "in")) +
  labs(y = "Density", color = "Population", x = "Time (hr)")
dev.off()

png("figS21_empiricalcurves.png", width = 8, height = 5,
    units = "in", res = 300)
ggplot(data = filter(gcdata_tidy, bacteria != "Blank"), 
       aes(x = Time, y = smoothed, color = file)) +
  geom_line(aes(group = paste(file, Well), lty = phage_added),
            alpha = 0.7) +
  facet_wrap(~plot_names, nrow = 2) +
  scale_color_manual(values = my_cols[1:3],
                     name = "Batch", labels = 1:3) +
  scale_linetype_manual(name = "Phage added?", values = 2:1) +
  scale_x_continuous(breaks = c(0, 12, 24)) +
  labs(x = "Time (hr)", y = "OD600", subtitle = "Bacterial strain") +
  theme_bw()
dev.off()

fs13a <- 
  ggplot(data = filter(gcdata_sum_sum_sum, phage_added == "Yes"),
         aes(x = eop_mean_mean, y = peak_dens_mean_mean)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_linerange(aes(ymin = peak_dens_mean_min, ymax = peak_dens_mean_max)) +
  scale_x_log10() +
  labs(x = "Efficiency of plaquing", y = "Peak density (OD600)") +
  theme_bw()

fs13b <- 
  ggplot(data = filter(gcdata_sum_sum_sum, phage_added == "Yes"),
         aes(x = eop_mean_mean, y = peak_time_mean_mean)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_linerange(aes(ymin = peak_time_mean_min, ymax = peak_time_mean_max)) +
  scale_x_log10() +
  scale_y_continuous(breaks = c(0, 6, 12, 18, 24)) +
  labs(x = "Efficiency of plaquing", y = "Peak time (hr)") +
  theme_bw()

fs13c <- 
  ggplot(data = filter(gcdata_sum_sum_sum, phage_added == "Yes"),
         aes(x = eop_mean_mean, y = auc_mean_mean)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_linerange(aes(ymin = auc_mean_min, ymax = auc_mean_max)) +
  scale_x_log10() +
  labs(x = "Efficiency of plaquing", y = "Area under the curve\n(OD600 hrs)") +
  theme_bw()

fs13d <- 
  ggplot(data = filter(gcdata_sum_sum_sum, phage_added == "Yes"),
         aes(x = eop_mean_mean, y = death_slope_mean_mean)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_linerange(aes(ymin = death_slope_mean_min, ymax = death_slope_mean_max)) +
  scale_x_log10() +
  labs(x = "Efficiency of plaquing", y = "Minimum death slope\n(OD600/hr)") +
  theme_bw()

png("figS22_EOPvmetrics.png", width = 6, height = 5,
    units = "in", res = 300)
cowplot::plot_grid(fs13a, fs13b, fs13c, fs13d,
                   nrow = 2, labels = "AUTO", align = "hv", axis = "tblr")
dev.off()