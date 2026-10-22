# R script to read in a metadata file with few columns and add missing columns with NA values or
# default values, then write the updated metadata to a new file.

### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE THAT AGES MAY BE REMOVED HERE, AND SHOULD BE ADDED IN EXTRA_META FILE!!

library(seqinr)

setwd("C:\\Users\\Emma\\wsl\\enteroviruses\\ev_githubs\\enterovirus_d68-oct2025\\data")
#orig_meta_file <- "ENPEN_2025-11-19-raw.csv"
#orig_seqs_file <- "ENPEN_2025-11-19.fasta"
#new_meta_file <- "ENPEN_2025-11-19.csv"

#read in age file to check ages are there too
age_file <- "age-data.tsv"
age_data <- read.csv(age_file, sep="\t", stringsAsFactors = FALSE)

orig_meta_file <- "ENPEN_2025-11-24-raw.csv"
orig_seqs_file <- "ENPEN_2025-11-24-raw.fasta"
#new_meta_file <- "ENPEN_2025-11-24.csv" #note this name doesn't reflect the true date when created, now
#new_seqs_file <- "ENPEN_2025-11-24.fasta" #note this name doesn't reflect the true date when created, now
#new_meta_file <- "ENPEN_2026-07-19.csv"  #after addition new seqs & meta 
#new_seqs_file <- "ENPEN_2026-07-19.fasta" #after addition new seqs & meta 
new_meta_file <- "ENPEN_2026-10-10.csv"  #after addition new seqs & meta
new_seqs_file <- "ENPEN_2026-10-10.fasta" #after addition new seqs & meta

#sequences from rotterdam 
nl_seq_file <- "NL_Rotterdam_2025-11-24.fasta"
nl_meta_file <- "NL_Rotterdam_2025-11-24.csv"

#sequences from Sweden - some replace those in orig_meta_file
# and some are new!
swe_seq_file <- "swe_replace_2025-11-25.fasta"
swe_meta_file <- "swe_replace_2025-11-25.csv"

# add some sewage sequences from the UK - want to be clear on origin here
# RECEIVED NEW FILES - use the latest files instead
uk_seq_file <- "D68_UK_sewage_2026-01-17.fasta" #"D68_UK_sewage_2025-12-19.fasta"
uk_meta_file <- "D68_UK_sewage_2026-01-17.csv" #"D68_UK_sewage_2025-12-19.csv"

#add new sequence from Italy 
ita_seq_file <- "EV-D68-Italy-LOM-MI-26-00164-2025.fasta"
ita_meta_file <- "EV-D68-Italy-LOM-MI-26-00164-2025.csv"

#add new sequences from Iceland, UK, belgium (note these have been combined by Emma)
ukicebel_seq_file <- "Uk-Iceland-Belg_2026-05-10.fasta"
ukicebel_meta_file <- "Uk-Iceland-Belg_2026-05-10.csv"

#add new sequences from Netherlands
nl1b_seq_file <- "NL1b_2026-06-01.fasta"
nl1b_meta_file <- "NL1b_2026-06-01.csv"

#add final set of new sequences from Esp, UK, Den, Fra, NL, Iceland - 19 Jul
jul19_seq_file <- "esp-uk-den-fra-nl-ice_2026-Jul.fasta"
jul19_meta_file <- "esp-uk-den-fra-nl-ice_2026-Jul_raw.csv"

#add final set of sequences from netherlands - 19 Sept
nl2_seq_file <- "NL_2026-09-19.fasta"
nl2_meta_file <- "NL_2026-09-19.csv"


##### Start processing ######

# read in original metadata
orig_meta <- read.csv(orig_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- orig_meta

# read in rotterdam metadata and append to updated_meta
nl_meta <- read.csv(nl_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- rbind(updated_meta, nl_meta)

#read in sweden metadata and update or append to updated_meta
swe_meta <- read.csv(swe_meta_file, sep=";", stringsAsFactors = FALSE)
# for each row in swe_meta, check if strain exists in updated_meta  
for (i in 1:nrow(swe_meta)) {
  strain <- swe_meta$seq_name[i]
  if (strain %in% updated_meta$seq_name) {
  # if exists, update the row
  updated_meta[updated_meta$seq_name == strain, ] <- swe_meta[i, ]
  } else {
  # if not exists, append the row
  updated_meta <- rbind(updated_meta, swe_meta[i, ])
  }
}

#read in UK sewage metadata and append to updated_meta
uk_meta <- read.csv(uk_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- rbind(updated_meta, uk_meta)

#read in Italy metadata and append to updated_meta
ita_meta <- read.csv(ita_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- rbind(updated_meta, ita_meta)

#read in UK/Iceland/Belgium metadata and append to updated_meta
ukicebel_meta <- read.csv(ukicebel_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- rbind(updated_meta, ukicebel_meta)

#read in Netherlands metadata and append to updated_meta
nl1b_meta <- read.csv(nl1b_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- rbind(updated_meta, nl1b_meta)

#read in july 19 metadata and append to updated_meta
jul19_meta <- read.csv(jul19_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- rbind(updated_meta, jul19_meta)

#read in sept data from NL
nl2_meta <- read.csv(nl2_meta_file, sep=";", stringsAsFactors = FALSE)
updated_meta <- rbind(updated_meta, nl2_meta)

# list of necessary columns
necessary_columns <- c("Lab.ID","strain","virus","subgenogroup","accession",
				"seq.len","age","sex","date","date_added","host","genbank.host",
				"country","region","collab_country","moltype","authors","title","url",
				"paper_url")

# rename existing columns if needed
# seq_name should be Lab.ID & strain
# collection_date should be date
if ("seq_name" %in% colnames(updated_meta)) {
  updated_meta$Lab.ID <- updated_meta$seq_name
  updated_meta$strain <- updated_meta$seq_name
  updated_meta$seq_name <- NULL
}
if ("collection_date" %in% colnames(updated_meta)) {
  updated_meta$date <- updated_meta$collection_date
  updated_meta$collection_date <- NULL
}

# add missing columns with NA values or default values
for (col in necessary_columns) {
  if (!(col %in% colnames(updated_meta))) {
	if (col == "virus") {
	  updated_meta[[col]] <- "Enterovirus D68"
	} else if (col == "moltype") {
	  updated_meta[[col]] <- "genomic RNA"
	} else if (col == "host" || col == "genbank.host") {
	  updated_meta[[col]] <- "Human"
	} else {
	  updated_meta[[col]] <- NA	
	}
  }
}

#change date_added to current date
updated_meta$date_added <- Sys.Date()

# add the regions column based on country and the regions file
#regions_df contains columns country and region to match country to region
regions_df <- read.csv("../config/geo_regions.tsv", sep="\t", stringsAsFactors = FALSE)
# update region column based on country
for (i in 1:nrow(updated_meta)) {
  country <- updated_meta$country[i]
  region <- regions_df$region[regions_df$country == country]
  if (length(region) > 0) {
	updated_meta$region[i] <- region
  } else {
	updated_meta$region[i] <- NA
  }
}

### TRIM MESSY ENDS
# Let's do something about sequences that have suspected 'messy ends' - there may be a better way to
# automate detection, but for now we'll do this manually

#create table of sequences to trim: name, number bases from start, number of bases from end (to trim)
# First, trim up to 130 in front & after 490 at end, as that seems where clusters of bad happen
trim_seqs <- data.frame(
  strain = c("HUVH-111881340",
             "HUVH-102693808",
             "HUAV-198335645",
             "HUAV-198335026",
             "NL_RTM/EV-D68/4/2025",
             "HUVH-102698999",
             "Portugal/1/2025",
             "Portugal/2/2025",
             "HUVH-102508169",
             "HUVH-102510181",
             "HUAV-198320175",
             "HUJO-560002275",
             "HUJT-198343418",
             #custom trimming
             "NL_RTM/EV-D68/3/2025" #trim ambig bases at end
             ),
  trim_start = c(27, 75, 16, 19, 25, 16, 25, 25, 3, 12, 12, 16, 8, 0),
  trim_end =   c(42, 0,   0,  0,  0,  0,  7,  7, 5,  2,  0,  0, 0, 7)
)

# read in sequences

# main sequence file
seqs <- read.fasta(orig_seqs_file, as.string = TRUE, forceDNAtolower = FALSE)
# rotterdam sequences
nl_seqs <- read.fasta(nl_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# uk sequences
uk_seqs <- read.fasta(uk_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# Italy sequence
ita_seqs <- read.fasta(ita_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# UK/Iceland/Belgium sequences
ukicebel_seqs <- read.fasta(ukicebel_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# Netherlands new sequences
nl1b_seqs <- read.fasta(nl1b_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# July 19 new sequences
jul19_seqs <- read.fasta(jul19_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# September 19 new sequences from netherlands
nl2_seqs <- read.fasta(nl2_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# combine sequences
seqs <- c(seqs, nl_seqs, uk_seqs, ita_seqs, ukicebel_seqs, nl1b_seqs, jul19_seqs,nl2_seqs)

#for sweden sequences, need to replace some sequences and add some new ones
swe_seqs <- read.fasta(swe_seq_file, as.string = TRUE, forceDNAtolower = FALSE)
# create a named list for easy replacement
swe_seqs_named <- setNames(swe_seqs, sapply(swe_seqs, function(x) attr(x, "name")))
# replace or add sequences in seqs
for (i in 1:length(swe_seqs_named)) {
  strain <- names(swe_seqs_named)[i]
  if (strain %in% sapply(seqs, function(x) attr(x, "name"))) {
  # if exists, replace the sequence
  seqs[[which(sapply(seqs, function(x) attr(x, "name")) == strain)]] <- swe_seqs_named[[i]]
  } else {
  # if not exists, append the sequence
  seqs <- c(seqs, swe_seqs_named[i])
  }
}
orig_seqs <- seqs

# trim sequences as per trim_seqs table
for (i in 1:nrow(trim_seqs)) {
  strain <- trim_seqs$strain[i]
  trim_start <- trim_seqs$trim_start[i]
  trim_end <- trim_seqs$trim_end[i]
  # find the sequence in seqs
  seq_index <- which(sapply(seqs, function(x) attr(x, "name")) == strain)
  if (length(seq_index) > 0) {
  seq <- seqs[[seq_index]]
  seq_length <- nchar(seq[[1]])
  # trim the sequence
  new_start <- trim_start + 1
  new_end <- seq_length - trim_end
  trimmed_seq <- substr(seq[[1]], new_start, new_end)
  # update the sequence in seqs
  seqs[[seq_index]][[1]] <- trimmed_seq
  }
}




#### CHECKS!!! #####

#do a quick check to ensure we can find the sequences for all strains in the fasta file

# get sequence names
seq_names <- sapply(seqs, function(x) attr(x, "name"))
# check if all strains in metadata have corresponding sequences
missing_seqs <- setdiff(updated_meta$strain, seq_names)
if (length(missing_seqs) > 0) {
  cat("Warning: The following metadata strains are missing sequences:\n")
  print(missing_seqs)
} else {
  cat("All metadata strains have corresponding sequences.\n")
}
# check if all sequences have corresponding metadata entries
missing_meta <- setdiff(seq_names, updated_meta$strain)
if (length(missing_meta) > 0) {
  cat("Warning: The following sequences are missing metadata entries:\n")
  print(missing_meta)
} else {
  cat("All sequences have corresponding metadata entries.\n")
}

#specific instructions to fix any names per file
# Italy ones not needed anymore as fixed by Robert 24 Nov
# for ENPEN 2025 data from Nov 2025:
# Italy-LOM-MI-25-05200/2025, should be Italy-LOM-MI-25-05200
# Italy-LOM-MI-25-05201/2025, should be Italy-LOM-MI-25-05201
# Italy-LOM-MI-25-05403/2025, should be Italy-LOM-MI-25-05403
# Italy-LOM-MI-25-05418/2025, should be Italy-LOM-MI-25-05418
# These need to be fixed one by one as other sequences hvae /2025 as part of their name
#updated_meta$strain[updated_meta$strain == "Italy-LOM-MI-25-05200/2025"] <- "Italy-LOM-MI-25-05200"
#updated_meta$strain[updated_meta$strain == "Italy-LOM-MI-25-05201/2025"] <- "Italy-LOM-MI-25-05201"
#updated_meta$strain[updated_meta$strain == "Italy-LOM-MI-25-05403/2025"] <- "Italy-LOM-MI-25-05403"
#updated_meta$strain[updated_meta$strain == "Italy-LOM-MI-25-05418/2025"] <- "EV-D68-Italy-LOM-MI-25-05418"


#print out the number of entries in metadata and the number of sequences
cat("Number of metadata entries:", nrow(updated_meta), "\n")
cat("Number of sequences:", length(seqs), "\n")

#how many have age data?
cat("This many sequences appear in the age file:", length(which(updated_meta$strain %in% age_data$accession)), "\n")
enpenEntries <- age_data[which(age_data$accession %in% updated_meta$strain),]
#which ones actually have entries in the age column
#count how many have "" in the enpenEntries$age column
cat("This many sequences have age data in the age column:", length(which(enpenEntries$age != "")), "\n")
#count how many are marked as sewage in sampleType column?
cat("This many sequences are marked as sewage in the sampletype column:", length(which(enpenEntries$sampletype == "sewage")), "\n")
# display how many are clinical then
cat("This many sequences are therefore clinical in the sampletype column:", length(which(enpenEntries$sampletype != "sewage")), "\n")
# how many are <300 bases?
cat("This many sequences are <300 bases:", length(which(nchar(sapply(seqs, function(x) x[[1]])) < 300)), "\n")

#check any ages missing
#updated_meta$strain[which(!(updated_meta$strain %in% age_data$accession))]

#rearrange columns to match necessary_columns order
updated_meta <- updated_meta[, necessary_columns]
# write updated metadata to new file
write.csv(updated_meta, new_meta_file, row.names = FALSE, na = "")

#write updated sequences to new fasta file
write.fasta(sequences = seqs, names = sapply(seqs, function(x) attr(x, "name")), 
            file.out = new_seqs_file)
