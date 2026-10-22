#A cleaned up R file that does the age Fisher test analysis

# Produces: 
# 1.  Fisher tests between years for proportion of adult vs child cases
# 2.  Plot of proportion of adult cases by year, with number of sequences with age data indicated (Supp Fig 3)



#What's slightly confusing about this file is that it uses the Nextstrain output files to include all sequences included in the phylogeny (those downloaded online)
#with clade and etc - but this excludes some 2025 sequences that were excluded from the phylogeny for various reasons (too short, duplicates, etc).  #
#So these are manually added back in, and get their subgenotype and age data from other sources.

# date of last phylogenetic run (used as one method to identify the '2025' sequences)
date_phylo_run <- "2026-10-10"

#Read in the metadata file of the Nextstrain analysis, as it contains all
#data used, with subgenotype (from the phylogeny) and age.
# However this will *not* contain any sequences excluded, e.g., for length or 
setwd("C:/Users/Emma/wsl/enteroviruses/ev_githubs/enterovirus_d68-oct2025/vp1/results")
meta <- read.csv("metadata_subgeno_seqlen_dupcol_extramet.tsv", sep="\t", header=T)

#convert dates, and also catch month and year in case of dates with missing days (-XX)
meta$dateraw <- meta$date
meta$year <- substr(meta$dateraw,1,4)
meta$month <- substr(meta$dateraw,6,7)
meta$date <- as.Date(meta$date, format="%Y-%m-%d")

# read in the original 2025 metadata file as well, as it contains all sequences from 2025
# from this we can get the ones excluded from the phylogeny, and manually add their subgenotype
# running in Nextcldae, and age 

# original metadata file (no subgenotype, no age)
enpenMeta <- read.csv("../../data/ENPEN_2026-07-19.csv", header=T, stringsAsFactors=F)
# read in age data
agedata <- read.csv("../../data/age-data.tsv", header=T, stringsAsFactors=F, sep="\t")

# Find out wihch sequences are from ENPEN in 2025 - these were 'added' in most recent
# phylogenetic run 2026-07-19
# Find out which sequences from 2025 were excluded from the phylogeny
meta2025 <- meta[which(meta$year == "2025" & meta$date_added == date_phylo_run),]
# 230 sequences
# There should be 242 - so 12 are missing, as expected - but exclude the bad one.

# get the age data for those that are ENPEN 2025 sequences
ageDataEnpen <- agedata[which(agedata$accession %in% enpenMeta$strain),] #238

# how many are in the age-data file, but not in the meta file?
ageDataEnpen[which(ageDataEnpen$age != "" & !(ageDataEnpen$accession %in% meta$strain) ),]
#12 sequences - as expected.
# These are 11 that are <300 bp and one that was excluded manually

#remove the one we're now excluding overall: HUVH-111881340
ageDataEnpen <- ageDataEnpen[which(ageDataEnpen$accession != "HUVH-111881340"),]

## These need to be added to meta
# get the accessions
missingAccessions <- ageDataEnpen[which(ageDataEnpen$age != "" & !(ageDataEnpen$accession %in% meta$strain) ),"accession"]

# get the data
missEnpenData <- enpenMeta[which(enpenMeta$strain %in% missingAccessions),]
# add the age ranges from the age-data file
missEnpenData$age <- ageDataEnpen[match(missEnpenData$strain, ageDataEnpen$accession), "age"]

# convert into age_range3 (<18y or >=18y)
for (i in 1:nrow(missEnpenData)) {
  #if isn't a range, create two parts of the range that are that age
  if (!grepl("-", missEnpenData$age[i])) {
	agerng <- c(missEnpenData$age[i], missEnpenData$age[i])
  } else{
    #break into parts of the range if a range
    agerng <- strsplit(missEnpenData$age[i],"-")[[1]]
  }
  if (missEnpenData$age[i] == "") {
	missEnpenData$age_range3[i] <- ""
  } else if (as.numeric(agerng[1]) < 18 & as.numeric(agerng[2]) < 18) {
	missEnpenData$age_range3[i] <- "<18y"
  } else if (as.numeric(agerng[1]) >= 18 & as.numeric(agerng[2]) >= 18) {
	missEnpenData$age_range3[i] <- ">=18y"
  } else { #if age range spans both sides of 18, can't categorize
	missEnpenData$age_range3[i] <- ""
  }
}

#ran these 11 through Nextclade and all are firmly A2 
missEnpenData$subgenogroup <- "A2"

# add data info 
missEnpenData$date_added <- date_phylo_run
missEnpenData$dateraw <- missEnpenData$date
missEnpenData$year <- substr(missEnpenData$dateraw,1,4)
missEnpenData$month <- substr(missEnpenData$dateraw,6,7)
missEnpenData$date <- as.Date(missEnpenData$date, format="%Y-%m-%d")


# add to main metadata
# Select only the columns we need
cols_needed <- c("strain", "subgenogroup", "age", "date", "date_added", 
                 "age_range3", "dateraw", "year", "month", "country")

# Select from each data frame
meta_subset <- meta[, cols_needed]
missEnpen_subset <- missEnpenData[, cols_needed]

# Combine them
allmeta <- rbind(meta_subset, missEnpen_subset)

################################
# Check - do we have all the data we are expecting?
# we expect 241 sequences from ENPEN 2025 (with 1 excluded)
# We expect 114 of these to be clinical, and have age range data
# Therefore a final dataset of 114 ENPEN 2025 sequences with age_range3 data

meta2025 <- allmeta[which(allmeta$year == "2025" & allmeta$date_added == date_phylo_run),]
# 241 sequences
length(which(meta2025$age_range3 != "")) # 114 sequences with age range data


#################################
# Do the Fishers tests - first divide into years we're interested in,
# and age and subgenotype groups

# Put counts into a table we can draw from
yrs <- c("2016", "2018", "2021", "2025")

yrCts <- matrix(nrow=4, ncol=2)
colnames(yrCts) <- c("adult", "child")
rownames(yrCts) <- yrs

for(i in 1:4) {
    m <- allmeta[which(allmeta$year == yrs[i]),]
    yrCts[i,1] <- length(which(m$age_range3 == ">=18y"))
    yrCts[i,2] <- length(which(m$age_range3 == "<18y"))
}

# Output the counts and results for each year, in an easy to read format

# Calculate percentages and totals for each year
calc_stats <- function(year) {
  adult <- as.numeric(yrCts[year, "adult"])
  child <- as.numeric(yrCts[year, "child"])
  total <- adult + child
  pct <- round(adult / total * 100,2)
  return(list(pct = pct, n = total, adult=adult))
}

stats_2018 <- calc_stats("2018")
stats_2016 <- calc_stats("2016")
stats_2021 <- calc_stats("2021")
stats_2025 <- calc_stats("2025")

# Build the sentence
result <- paste0(
  "In 2018, the previous year when A2 was more dominant, ",
  "of those cases with ages, ", stats_2018$pct, "% were adults (n=", stats_2018$n, "). ",
  "In a more 'typical' year with B3 dominance, the percent of adult cases is much lower ",
  "(2016, ", stats_2016$pct, "% with n=", stats_2016$n, " and 2021, ", stats_2021$pct, "% with n=", stats_2021$n, "). ",
  "In 2025, is ", stats_2025$pct, "%, though number of sequences with ages is fairly small (n=", stats_2025$n, ")."
)

print(result)

###### Now do the Fishers tests
# Comparison 1: 2025 vs. 2021

table_2025_vs_2021 <- matrix(c(yrCts["2025","adult"], yrCts["2025","child"],
                               yrCts["2021","adult"], yrCts["2021","child"]),
                               nrow = 2,
                               byrow = FALSE,
                               dimnames = list(
                               Age = c("Adults", "Children"),
                               Year = c("2025", "2021")
                              ))
# view the table
print(table_2025_vs_2021)

# run Fishers exact test
fisher_2025_vs_2021 <- fisher.test(table_2025_vs_2021)
print(fisher_2025_vs_2021)

# comparison 2: 2025 vs. 2018
table_2025_vs_2018 <- matrix(c(yrCts["2025","adult"], yrCts["2025","child"],
                               yrCts["2018","adult"], yrCts["2018","child"]),
                               nrow = 2,
                               byrow = FALSE,
                               dimnames = list(
                               Age = c("Adults", "Children"),
                               Year = c("2025", "2018")
                              ))

# view the table
print(table_2025_vs_2018)

# run Fishers exact test
fisher_2025_vs_2018 <- fisher.test(table_2025_vs_2018)
print(fisher_2025_vs_2018)


#############################################################################
#############################################################################
#############################################################################
# Create plot of ages < & > 18 by year, where available (not available for all)

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
# Prepare data for age range by year
# Filter out rows without valid year or age data
valid_age_data <- allmeta[allmeta$year != "" & allmeta$age_range3 != "", ]

# Count sequences per year and age range
year_age <- as.data.frame(table(year = valid_age_data$year, 
                                age_range = valid_age_data$age_range3))
colnames(year_age) <- c("year", "age_range", "n")

# Get all unique years from 2013-2025, sorted
all_years <- sort(as.character(2013:2025))

# Create data frame with all combinations of years and age ranges
all_combos_age <- expand.grid(year = all_years, 
                              age_range = c("<18y", ">=18y"), 
                              stringsAsFactors = FALSE)

# Merge to include all combinations, fill missing with 0
year_age <- merge(all_combos_age, year_age, 
                  by = c("year", "age_range"), all.x = TRUE)
year_age$n[is.na(year_age$n)] <- 0

# Convert year to ordered factor
year_age$year <- factor(year_age$year, levels = all_years, ordered = TRUE)

# Calculate total per year and proportion
year_totals <- tapply(year_age$n, year_age$year, sum)
year_age$total <- year_totals[as.character(year_age$year)]
year_age$prop <- ifelse(year_age$total == 0, 0, 
                        year_age$n / year_age$total)


# Create stacked bar plot by year and age range
p_age_by_year <- ggplot(year_age, aes(x = year, y = prop, fill = age_range)) +
  geom_col(position = "stack") +
geom_text(data = year_age[year_age$age_range == ">=18y", ], 
          aes(x = year, y = 1.02, label = paste0("n=", total)),
          vjust = 0, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("<18y" = "#1f77b4", ">=18y" = "#ff7f0e")) +
  scale_y_continuous(labels = percent, limits = c(0, 1.1)) +
  labs(x = "Year", y = "Proportion of Samples", fill = "Age Range", 
       title = "Proportion of Samples by Year and Age Range in Europe") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_age_by_year)

ggsave("../../figs/age_by_year.pdf", p_age_by_year, width = 10, height = 6)



