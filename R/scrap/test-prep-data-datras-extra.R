library(DATRASextra)
library(dplyr)
library(ggplot2)


# Main tutorial; https://github.com/tokami/DATRASextra/blob/main/vignettes/datrasextra-tutorial.Rmd
#list_surveys()
#plot_surveys()

# Download data
surveys <- c("BTS", "FR-CGFS", "IE-IGFS", "NS-IBTS", "SCOWCGFS", "SWC-IBTS")
  
# Create temporary directory
#tmp <- tempdir()
dir <- paste0(here::here(), "/data-dir")

# Download
#download_datras(surveys = surveys, years = 1992:2024, dir = dir)

# Load data
surv0 <- read_datras(file.path(dir, surveys))
surv0[["HH"]]
#cod <- subset(surv0, )

# Clean data (only valid hauls etc), impute missing depths
surv <- clean_datras(surv0)



# From here, species specific!
cod <- subset(surv, Valid_Aphia == 126436)

# Check for outliers
cod <- check_outliers(cod, pct = TRUE)

head(attr(cod, "outlier_report"))
 
# Plot data (returns error)
# FIXME error
plot(cod)
 
# # To calculate total numbers and weight by haul, the information in the `HL` table
# # first needs to be raised to the haul level. A convenient first step is to add
# # numbers-at-length:
# 
# surv <- add_numbers_at_length(surv)
# 
# FIXME: Throws error when multiple classes!
# class(surv)
# class(surv) <- "DATRASraw"
# surv <- add_numbers_at_length(surv)


# Custom length tutorial: https://github.com/tokami/DATRASextra/blob/main/vignettes/articles/custom-length-classes.Rmd
lpars <- check_lengths(cod)

class(cod) <- "DATRASraw"

# Add numbers at length before converting to biomass. Then we add the custom size classes
cod <- add_numbers_at_length(cod)

wpars <- check_weights(cod)

cod <- add_weight_at_length(cod)
cod[["HH"]]

length_cuts <- c(0, 35, Inf)

cod <- add_total_numbers_by_haul(cod, length_cuts = length_cuts)
cod <- add_total_weight_by_haul(cod, length_cuts = length_cuts)

head(cod[["HH"]][["HaulN"]])
head(cod[["HH"]][["HaulWgt"]])

colnames(cod[["HH"]][["HaulN"]]) <- c("small", "large")
colnames(cod[["HH"]][["HaulWgt"]]) <- c("small", "large")

# ncols <- ncol(cod[["HH"]][["HaulWgt"]])
# 
# par(mfrow = n2mfrow(ncols, asp = 2),
#     mar = c(3, 3, 2, 2),
#     oma = c(2, 2, 0, 0))
# 
# for (i in seq_len(ncols)) {
#   plot(cod[["HH"]]$lon, cod[["HH"]]$lat,
#        type = "n",
#        xlab = "", ylab = "",
#        main = colnames(cod[["HH"]]$HaulN)[i])
#   
#   ind <- which(cod[["HH"]]$HaulN[, i] > 0)
#   
#   points(cod[["HH"]]$lon[ind], cod[["HH"]]$lat[ind],
#          cex = cod[["HH"]]$HaulN[, i][ind] /
#            max(cod[["HH"]]$HaulN[, i][ind]) * 3,
#          col = i)
# }
# 
# mtext("Longitude", 1, outer = TRUE)
# mtext("Latitude", 2, outer = TRUE)

cod

t <- cod[["HH"]]
# t |> 
#   dplyr::summarise(n = dplyr::n(), .by = haul.id) |> 
#   dplyr::distinct(n)

cod[["HH"]][["HaulWgt"]]

d <- cod[["HH"]]

names(t)

str(t)

ggplot(t, aes(lon, lat ))

# clean data:
haul_dur





