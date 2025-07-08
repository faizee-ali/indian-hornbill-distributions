library(auk)
library(tidyverse)
setwd("D:/NCF/")

##making the eBird data from single species into a big data frame of all 5 species
ebd_brn<- read_tsv("ebd/ebd_IN_brnhor1_smp_relFeb-2025/ebd_IN_brnhor1_smp_relFeb-2025.txt")
ebd_gre<- read_tsv("ebd/ebd_IN_grehor1_smp_relFeb-2025/ebd_IN_grehor1_smp_relFeb-2025.txt")
ebd_orp<- read_tsv("ebd/ebd_IN_orphor1_relFeb-2025/ebd_IN_orphor1_relFeb-2025.txt")
ebd_run<- read_tsv("ebd/ebd_IN_runhor1_smp_relFeb-2025/ebd_IN_runhor1_smp_relFeb-2025.txt")
ebd_wre<- read_tsv("ebd/ebd_IN_wrehor1_relFeb-2025/ebd_IN_wrehor1_relFeb-2025.txt")

all <- rbind(ebd_brn, ebd_gre, ebd_orp, ebd_run, ebd_wre)
all_ap <- filter(all, all$STATE == "Arunachal Pradesh")

write_tsv(all, file = "all-hornbills.txt")
write_tsv(all_ap, file = "all-hornbills-ap.txt")
##keep data inside NE only
library(rgdal)
library(sf)
library(readr)
library(dplyr)

ooty <- st_read("data/shapefiles/NE_smaller.kml")
ebd_sampling = read.delim("data/ebd_sampling_NE.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F, na.strings = c(""," ",NA))
ebd_sf <- st_as_sf(ebd_sampling, coords = c("longitude", "latitude"), crs = 4326)


ooty_sf <- st_as_sf(ooty)
ooty_sf <- st_transform(ooty_sf, crs = st_crs(ebd_sf))


points_in_ooty <- st_within(ebd_sf, ooty_sf, sparse = FALSE)
# subset data frame
ebd_in_ooty <- ebd_sampling[points_in_ooty, ]
#saving the WG subsetted sampling file  
write_tsv(
    ebd_in_ooty,
    "ebd_sampling_NE.txt",
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    escape = c("double", "backslash", "none"),
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress(),
)


# doing the same for gbif data
#first just convert to character when needed filter then convert to posix 

gbif_occ1<- read_tsv("data/gbif/0001933-250402121839773/occurrence.txt")
gbif_occ1$eventDate <- as.character(gbif_occ1$eventDate)
gbif_occ1$eventTime <- as.character(gbif_occ1$eventTime)

gbif_occ2<- read_tsv("data/gbif/0001939-250402121839773/occurrence.txt")
gbif_occ2$eventDate <- as.character(gbif_occ2$eventDate)
gbif_occ2$eventTime <- as.character(gbif_occ2$eventTime)

gbif_occ3<- read_tsv("data/gbif/0001941-250402121839773/occurrence.txt")
gbif_occ3$eventDate <- as.character(gbif_occ3$eventDate)
gbif_occ3$eventTime <- as.character(gbif_occ3$eventTime)

gbif_occ4<- read_tsv("data/gbif/0001947-250402121839773/occurrence.txt")
gbif_occ4$eventDate <- as.character(gbif_occ4$eventDate)
gbif_occ4$eventTime <- as.character(gbif_occ4$eventTime)

gbif_occ5<- read_tsv("data/gbif/0001948-250402121839773/occurrence.txt")
gbif_occ5$eventDate <- as.character(gbif_occ5$eventDate)
gbif_occ5$eventTime <- as.character(gbif_occ5$eventTime)

gbif_all <- rbind(gbif_occ1, gbif_occ2, gbif_occ3, gbif_occ4, gbif_occ5)
gbif_all_ap <- gbif_all %>% filter(stateProvince == "Arunachal Pradesh")


library(ggplot2)
state_counts <- as.data.frame(table(gbif_all$stateProvince))  # Creates a frequency table
ggplot(state_counts, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Category", y = "Count", title = "Number of Occurrences in Each Category") +
    scale_fill_viridis_d() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

write_tsv(gbif_all, file = "data/all-hornbills-gbif.txt")
write_tsv(gbif_all_ap, file = "data/all-hornbills-ap-gbif.txt")

ebd_sampling = read.delim("data/gbif-all-hornbills-NE.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F, na.strings = c(""," ",NA))

ebd_sf <- st_as_sf(ebd_sampling, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)


ooty_sf <- st_as_sf(ooty)
ooty_sf <- st_transform(ooty_sf, crs = st_crs(ebd_sf))


points_in_ooty <- st_within(ebd_sf, ooty_sf, sparse = FALSE)
# subset data frame
ebd_in_ooty <- ebd_sampling[points_in_ooty, ]
#saving the WG subsetted sampling file  
write_tsv(
    ebd_in_ooty,
    "gbif-all-hornbills-NE.txt",
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    escape = c("double", "backslash", "none"),
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress(),
)
# for Hornbill watch data
library(readxl)

hw <- read_excel("data/Hornbill Watch_31 March 2025.xlsx")
hw_ap <- filter(hw, State == "Arunachal Pradesh")
state_counts <- as.data.frame(table(hw$State))  # Creates a frequency table
write_tsv(hw_ap, file = "data/all-hornbills-ap-HW.txt")

### extracting North-East hornbill data for 5 species
## ebd
all_ebd <- read_tsv("data/ebd/all-hornbills.txt")
NE_ebd <- filter(all_ebd, STATE %in% c("West Bengal", "Mizoram", "Assam", "Arunachal Pradesh", "Nagaland", "Manipur", "Tripura", "Meghalaya"))
write_tsv(NE_ebd, file = "data/NE_ebd.txt")
##HW
hw <- read_excel("data/HW/Hornbill Watch_31 March 2025.xlsx")
hw_NE <- filter(hw, State %in% c("West Bengal", "Mizoram", "Assam", "Arunachal Pradesh", "Nagaland", "Manipur", "Tripura", "Meghalaya"))
state_counts <- as.data.frame(table(hw_NE$State))  # Creates a frequency table
write_tsv(hw_NE, file = "data/HW_NE.txt")

ebd_sampling = read.delim("data/HW_NE.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F, na.strings = c(""," ",NA))
#verify on mapview
ebd_sampling <- ebd_sampling[!is.na(ebd_sampling$Latitude),]
ebd_sf <- st_as_sf(ebd_sampling, coords = c("Longitude", "Latitude"), crs = 4326)


ooty_sf <- st_as_sf(ooty)
ooty_sf <- st_transform(ooty_sf, crs = st_crs(ebd_sf))


points_in_ooty <- st_within(ebd_sf, ooty_sf, sparse = FALSE)
# subset data frame
ebd_in_ooty <- ebd_sampling[points_in_ooty, ]
#saving the WG subsetted sampling file  
write_tsv(
    ebd_in_ooty,
    "HW-all-hornbills-NE.txt",
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    escape = c("double", "backslash", "none"),
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress(),
)


##gbif
gbif_all <- read_tsv("data/gbif/all-hornbills-gbif.txt")
gbif_NE <- filter(gbif_all, stateProvince %in% c("West Bengal", "Mizoram", "Assam", "Arunachal Pradesh", "Nagaland", "Manipur", "Tripura", "Meghalaya"))
write_tsv(gbif_NE, file = "data/gbif_NE.txt")

