#Import functions ----
library(growth.curves.pkg)
library(ggplot2)
library(dplyr)

#Read data ----
trav_resis <- read.csv("./Empirical/trav-phage_resis_data.csv")

gcdata_10_25 <- 
  import_widemeasures("./Empirical/2021-10-25_Emma_Growth_Curve_Assay.csv",
                    startrow = 30, startcol = 2,
                    endrow = 126, endcol = 99)
colnames(gcdata_10_25)[2] <- "Temp"

gcdata_10_25_lng <- pivot_wide_longer(gcdata_10_25,
                                   id_cols = c("Time", "Temp"),
                                   values_to = "OD600")

design_10_25 <- 
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
                                  byrow = TRUE),
    phage = make_designpattern(c("Phage Added", "No Phage"),
                               rows = 2:7, cols = 2:11,
                               pattern = "111111222211111122221111112222111111222211100022001110002200",
                               byrow = TRUE),
  )
                                  
gcdata_10_25_lng <- merge_tidydesign_tidymeasures(design_10_25, gcdata_10_25_lng,
                                                  by = "Well", drop = TRUE)

gcdata_10_25_lng$Time <- lubridate::hms(gcdata_10_25_lng$Time)
gcdata_10_25_lng <- smooth_data(OD600 ~ Time,
                                data = gcdata_10_25_lng,
                                algorithm = "moving-average",
                                subset_by = gcdata_10_25_lng$Well,
                                window_width = 5)

ggplot(data = gcdata_10_25_lng, 
       aes(x = as.numeric(Time), y = as.numeric(OD600), 
           group = Well, color = phage)) +
  geom_point(size = 0.5) +
  geom_line(aes(y = fitted)) +
  facet_wrap(~bacteria) +
  NULL

gcdata_10_25_lng <- group_by(gcdata_10_25_lng,
                             Well, bacteria, phage)
gcdata_10_25_sum <- dplyr::summarize(gcdata_10_25_lng,
                              peak_index = find_local_extrema(fitted,
                                                 return_minima = FALSE,
                                                 width_limit = 11,
                                                 na.rm = TRUE)[1],
                              peak_time = Time[peak_index],
                              peak_dens = fitted[peak_index])

ggplot(data = gcdata_10_25_lng[gcdata_10_25_lng$bacteria == "PF", ], 
       aes(x = as.numeric(Time), y = as.numeric(fitted), 
           group = Well, color = phage)) +
  geom_line() +
  geom_point(data = gcdata_10_25_sum[gcdata_10_25_sum$bacteria == "PF", ],
             aes(x = peak_time, y = peak_dens),
             size = 2, color = "black") +
  facet_wrap(~bacteria) +
  NULL

ggplot(data = gcdata_10_25_sum,
       aes(x = bacteria, color = phage, y = as.numeric(peak_time))) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge(width = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
