#Import functions ----
library(growth.curves.pkg)
library(ggplot2)
library(dplyr)

setwd("./Empirical/")

#Read data ----
trav_resis <- read.csv("trav-phage_resis_data.csv")

#Growth curve runs:
#2021-10-15 Dift inoc dens
#2021-10-25 Dift isols (starting at 5e4 bact 5e3 phage)
#2021-10-27 Dift isols (starting at 1e5 bact 1e4 phage, 7xEGC not run)
#2021-11-03 Dift isols (starting at 1e5 bact 1e4 phage)
#2021-11-08 Dift isols (starting at 1e5 bact 1e4 phage, titered post-peak)
#    TODO: get these data stitched together across "plates"


gcdata <-
  import_widemeasures(
    files = c("2021-10-15_Emma_Growth_Curve.csv",
              "2021-10-25_Emma_Growth_Curve.csv",
              "2021-10-27_Emma_Growth_Curve.csv",
              "2021-11-03_Emma_Growth_Curve.csv"),
    startrow = 30, startcol = 2,
    endrow = 126, endcol = 99)
              
gcdata_lng <- pivot_wide_longer(gcdata,
                                id_cols = c("file", "Time", "T° 600"),
                                values_to = "OD600")
#Create designs ----
design_diftconcs <- 
  make_tidydesign(
    nrows = 8, ncols = 12,
    block_row_names = LETTERS[1:8],
    block_col_names = 1:12,
    wellnames_sep = "",
    init_bact = make_designpattern(c(10**5, 5*10**4, 10**4, 0),
                                   rows = 2:7, cols = 2:4,
                                   pattern = "444111222333222222",
                                   byrow = TRUE),
    init_bact = make_designpattern(c(10**5, 10**4),
                                   rows = 2:6, cols = 5:7,
                                   pattern = "111111111222222",
                                   byrow = TRUE),
    init_moi = make_designpattern(c(0, 0.1, 0.01),
                                  rows = 2:7, cols = 2:4,
                                  pattern = "111111111111222333",
                                  byrow = TRUE),
    init_moi = make_designpattern(c(0.1, 0.01, 0.001),
                                  rows = 2:6, cols = 5:7,
                                  pattern = "111222333111222",
                                  byrow = TRUE),
    #Row 7 cols 5:7 is actually empty but df rows will be dropped anyway
    bacteria = make_designpattern("PF",
                                  rows = 2:7, cols = 2:7,
                                  pattern = "1")
  )
design_diftconcs$init_phage <- design_diftconcs$init_bact * design_diftconcs$init_moi

design_isols <- 
  make_tidydesign(
    nrows = 8, ncols = 12,
    block_row_names = LETTERS[1:8],
    block_col_names = 1:12,
    wellnames_sep = "",
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

design_5e4_5e3 <-
  make_tidydesign(
    nrows = 8, ncols = 12,
    block_row_names = LETTERS[1:8],
    block_col_names = 1:12,
    wellnames_sep = "",
    init_bact = make_designpattern(c(5*10**4),
                                   rows = 2:7, cols = 2:11,
                                   pattern = "1"),
    init_phage = make_designpattern(c(5*10**3),
                                    rows = 2:7, cols = 2:7,
                                    pattern = "1"),
    init_phage = make_designpattern(c(0),
                                    rows = 2:7, cols = 8:11,
                                    pattern = "1")
  )

design_1e5_1e4 <-
  make_tidydesign(
    nrows = 8, ncols = 12,
    block_row_names = LETTERS[1:8],
    block_col_names = 1:12,
    wellnames_sep = "",
    init_bact = make_designpattern(c(10**5),
                                   rows = 2:7, cols = 2:11,
                                   pattern = "1"),
    init_phage = make_designpattern(c(10**4),
                                    rows = 2:7, cols = 2:7,
                                    pattern = "1"),
    init_phage = make_designpattern(c(0),
                                    rows = 2:7, cols = 8:11,
                                    pattern = "1")
  )

#Merge design and measures ----
gcdata_lng[["2021-10-15_Emma_Growth_Curve"]] <-
  merge_dfs(design_diftconcs,
            gcdata_lng[["2021-10-15_Emma_Growth_Curve"]],
            drop = TRUE)
gcdata_lng[["2021-10-25_Emma_Growth_Curve"]] <-
  merge_dfs(
    merge_dfs(design_isols, design_5e4_5e3),
    gcdata_lng[["2021-10-25_Emma_Growth_Curve"]],
    drop = TRUE)
gcdata_lng[["2021-10-27_Emma_Growth_Curve"]] <-
  merge_dfs(
    merge_dfs(design_isols, design_1e5_1e4),
    gcdata_lng[["2021-10-27_Emma_Growth_Curve"]],
    drop = TRUE)
gcdata_lng[["2021-11-03_Emma_Growth_Curve"]] <-
  merge_dfs(
    merge_dfs(design_isols, design_1e5_1e4),
    gcdata_lng[["2021-11-03_Emma_Growth_Curve"]],
    drop = TRUE)

gcdata_lng <- merge_dfs(gcdata_lng, collapse = TRUE)

gcdata_lng$Time <- lubridate::hms(gcdata_lng$Time)

#Smooth and summarize ----
gcdata_lng <- smooth_data(OD600 ~ Time,
                          data = gcdata_lng,
                          algorithm = "moving-average",
                          subset_by = paste(gcdata_lng$file, gcdata_lng$Well),
                          window_width = 5)

gcdata_lng <- group_by(gcdata_lng,
                             file, init_bact, init_phage, init_moi, Well, bacteria)
gcdata_sum <- dplyr::summarize(gcdata_lng,
                               peak_index = find_local_extrema(fitted,
                                                               return_minima = FALSE,
                                                               width_limit = 11,
                                                               na.rm = TRUE,
                                                               remove_endpoints = FALSE)[1],
                               peak_time = Time[peak_index],
                               peak_dens = fitted[peak_index])

gcdata_sum_sum <- dplyr::summarise(
  group_by(gcdata_sum, file, init_bact, init_phage, init_moi, bacteria),
  peak_time_avg = mean(as.numeric(peak_time)),
  peak_dens_avg = mean(peak_dens))

#Join in EOP data
trav_resis_sum <-
  dplyr::summarize(group_by(trav_resis, Proj, Pop, Treat, Timepoint, Isol),
                   mean_eop = mean(EOP),
                   bd = any(bd))
gcdata_sum_sum <- cbind(gcdata_sum_sum, 
                    data.frame(Proj = NA, Pop = NA, 
                               Treat = NA, Timepoint = NA, Isol = NA))
for (i in 1:nrow(gcdata_sum_sum)) {
  if(gcdata_sum_sum$bacteria[i] == "PF") {
    gcdata_sum_sum$Proj[i] <- "125"
    gcdata_sum_sum$Pop[i] <- "Anc"
    gcdata_sum_sum$Treat[i] <- "Anc"
    gcdata_sum_sum$Timepoint[i] <- 0
    gcdata_sum_sum$Isol[i] <- "Anc"
  } else if (substr(gcdata_sum_sum$bacteria[i], 1, 2) == "7x") {
    gcdata_sum_sum$Proj[i] <- "7x"
    gcdata_sum_sum$Pop[i] <- substr(gcdata_sum_sum$bacteria[i], 3, 3)
    gcdata_sum_sum$Treat[i] <- substr(gcdata_sum_sum$bacteria[i], 4, 4)
    gcdata_sum_sum$Timepoint[i] <- 14
    gcdata_sum_sum$Isol[i] <- substr(gcdata_sum_sum$bacteria[i], 5, 5)
  } else if (substr(gcdata_sum_sum$bacteria[i], 1, 3) == "125") {
    gcdata_sum_sum$Proj[i] <- "125"
    gcdata_sum_sum$Pop[i] <- substr(gcdata_sum_sum$bacteria[i], 4, 4)
    gcdata_sum_sum$Treat[i] <- substr(gcdata_sum_sum$bacteria[i], 5, 5)
    gcdata_sum_sum$Timepoint[i] <- 14
    gcdata_sum_sum$Isol[i] <- substr(gcdata_sum_sum$bacteria[i], 6, 6)
  }
}

gcdata_sum_sum <- left_join(gcdata_sum_sum, trav_resis_sum)


#Plot all data sloppily ----
for (filenm in unique(gcdata_lng$file)) {
  temp <- gcdata_lng[gcdata_lng$file == filenm, ]
  print(ggplot(data = temp, 
               aes(x = as.numeric(Time), y = fitted, 
                   color = bacteria, group = Well)) +
          geom_line() +
          geom_point(data = gcdata_sum[gcdata_sum$file == filenm, ],
                     aes(x = peak_time, y = peak_dens)) +
          NULL)
}

#Plots of run varying init dens & moi ----
temp <- gcdata_lng[gcdata_lng$file == "2021-10-15_Emma_Growth_Curve", ]
temp_sum <- gcdata_sum[gcdata_sum$file == "2021-10-15_Emma_Growth_Curve", ]
temp_sum_sum <- group_by(temp_sum[temp_sum$init_bact > 0 &
                                    temp_sum$peak_dens < 0.75, ], 
                         init_bact, init_moi) %>%
  dplyr::summarise(peak_dens = mean(peak_dens),
                   peak_time = mean(as.numeric(peak_time)))

print(ggplot(data = temp, 
             aes(x = as.numeric(Time)/3600, y = fitted+1, 
                 color = as.factor(init_bact), group = Well)) +
        geom_line() +
        scale_y_continuous(trans = "log10", name = "Smoothed OD600") +
        labs(color = "Initial Bacteria\n(cfu/mL)", x = "Time (h)") +
        NULL)

print(ggplot(data = temp[temp$init_bact > 0, ], 
             aes(x = as.numeric(Time)/3600, y = fitted+1, 
                 color = as.factor(init_moi), group = Well)) +
        geom_line() +
        facet_wrap(~init_bact) +
        scale_y_continuous(trans = "log10", name = "Smoothed OD600") +
        labs(color = "Initial MOI\n(pfu/cfu)", x = "Time (h)",
             subtitle = "Initial Bacteria (cfu/mL)") +
        NULL)

print(ggplot(data = temp_sum[temp_sum$init_bact > 0, ],
             aes(x = as.factor(init_moi), y = as.numeric(peak_time)/3600)) +
        geom_point(alpha = 0.5, position = position_jitter(0.2)) +
        facet_grid(~init_bact) +
        scale_y_continuous(trans = "log10", name = "Peak time (h)") +
        labs(subtitle = "Initial Bacteria (cfu/mL)",
             x = "Initial MOI (pfu/cfu)") +
        NULL)

print(ggplot(data = temp_sum[temp_sum$init_bact > 0, ],
             aes(x = as.factor(init_moi), y = as.numeric(peak_dens))) +
        geom_point(alpha = 0.5, position = position_jitter(0.2)) +
        facet_grid(~init_bact) +
        # geom_point(data = temp_sum_sum, size = 3, alpha = 0.5,
        #            aes(x = as.factor(init_moi), y = as.numeric(peak_dens))) +
        scale_y_continuous(trans = "log10", name = "Peak density (OD600)") +
        labs(subtitle = "Initial Bacteria (cfu/mL)",
             x = "Initial MOI (pfu/cfu)") +
        NULL)

#Plots of varying isolates ----
temp <- gcdata_lng[gcdata_lng$file %in% 
                     c("2021-10-25_Emma_Growth_Curve", 
                       "2021-10-27_Emma_Growth_Curve",
                       "2021-11-03_Emma_Growth_Curve"), ]
temp_sum <- gcdata_sum[gcdata_sum$file %in% 
                         c("2021-10-25_Emma_Growth_Curve", 
                           "2021-10-27_Emma_Growth_Curve",
                           "2021-11-03_Emma_Growth_Curve"), ]
temp_sum_sum <- gcdata_sum_sum[gcdata_sum_sum$file %in% 
                                 c("2021-10-25_Emma_Growth_Curve", 
                                   "2021-10-27_Emma_Growth_Curve",
                                   "2021-11-03_Emma_Growth_Curve"), ]

ggplot(data = temp[temp$bacteria != "Blank", ],
       aes(x = as.numeric(Time)/3600, y = fitted, lty = init_phage > 0,
           group = Well)) +
  geom_line() +
  scale_x_continuous(breaks = c(0, 12)) +
  facet_grid(file~bacteria,
             labeller = labeller(file = 
                                   as_labeller(c("2021-10-25_Emma_Growth_Curve" = "1", 
                                      "2021-10-27_Emma_Growth_Curve" = "2",
                                      "2021-11-03_Emma_Growth_Curve" = "3")))) +
  labs(x = "Time (h)", y = "OD600") +
  guides(lty = guide_legend(title = "Phage added?"))

ggplot(data = temp_sum[temp_sum$bacteria != "Blank", ],
       aes(x = file, y = peak_dens, 
           shape = as.factor(init_phage))) +
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~bacteria, nrow = 2) +
  scale_x_discrete(labels = 1:3) +
  labs(x = "Batch", y = "Peak density (OD600)") +
  guides(shape = guide_legend(title = "Initial Phage\n(pfu/mL)"))

ggplot(data = temp_sum[temp_sum$bacteria != "Blank", ],
       aes(x = file, y = as.numeric(peak_time)/3600, 
           shape = as.factor(init_phage))) +
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~bacteria, nrow = 2) +
  scale_x_discrete(labels = 1:3) +
  labs(x = "Batch", y = "Peak time (h)") +
  guides(shape = guide_legend(title = "Initial Phage\n(pfu/mL)"))

ggplot(data = temp_sum_sum[temp_sum_sum$bacteria != "Blank" &
                             temp_sum_sum$init_phage > 0 &
                             temp_sum_sum$peak_time_avg > 10000, ],
       aes(x = log10(mean_eop), y = peak_time_avg/3600, color = bacteria)) +
  geom_point(size = 2, alpha = 0.5) +
  facet_grid(~file, 
             labeller = as_labeller(c("2021-10-25_Emma_Growth_Curve" = "1", 
                                      "2021-10-27_Emma_Growth_Curve" = "2",
                                      "2021-11-03_Emma_Growth_Curve" = "3"))) +
  labs(subtitle = "Batch", y = "Peak time (h)", x = "log10(EOP)")

ggplot(data = temp_sum_sum[temp_sum_sum$bacteria != "Blank" &
                             temp_sum_sum$init_phage > 0 &
                             temp_sum_sum$peak_time_avg > 10000, ],
       aes(x = log10(mean_eop), y = peak_dens_avg, color = bacteria)) +
  geom_point(size = 2, alpha = 0.5) +
  facet_grid(~file, 
             labeller = as_labeller(c("2021-10-25_Emma_Growth_Curve" = "1", 
                                      "2021-10-27_Emma_Growth_Curve" = "2",
                                      "2021-11-03_Emma_Growth_Curve" = "3"))) +
  labs(subtitle = "Batch", y = "Peak density (OD600)", x = "log10(EOP)")

