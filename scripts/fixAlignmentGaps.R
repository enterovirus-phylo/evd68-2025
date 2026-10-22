args = commandArgs(trailingOnly=TRUE)

orig_seqs <- args[1]
fixed_seqs <- args[2]


####

library(seqinr)


##### Check for two types of gap in 300 and vp1 seqs
# It'll process this for full-genome too but it won't fix it currently...
# This should be replaced by a codon-aware alignment method, which will solve this problem.

##### Now extended to also fix gaps in nuc 418 area (AA 140 area)

seqs <- read.fasta(orig_seqs, as.string=T, forceDNAtolower=F)
seq_names <- attr(seqs,"name")

#As an example, locate a specific sequence in the file: EV68-B3-SEL-DIJ-21-E0-298-RSE2-21-1004-FRA21
#specific_seq <- "EV68-B3-SEL-DIJ-21-E0-298-RSE2-21-1004-FRA21"
#locate_seq <- grep(specific_seq, seq_names)
#seq <- seqs[[locate_seq]][1]

#First, 300 area of vp1

#Gap type number 1
withGap <- grep("G---TAA", seqs)

countr <- 0
for(wG in withGap){
    seq <- seqs[[wG]][1]
    findr <- regexpr("G---TAA", seq) #should be at 425 -- or can be at 428
    if(findr == 425) {
        newSeq <- gsub("G---TAA","GT---AA", seq)
        seqs[[wG]][1] <- newSeq
        countr = countr+1
    } else if(findr == 428) {
        newSeq <- gsub("G---TAA","GT---AA", seq)
        seqs[[wG]][1] <- newSeq
        countr = countr+1
    }
}
print(paste("replaced",countr,"type 1 gaps"))

#Gap type number 2

withGap <- grep("A---A", seqs)

countr <- 0
for(wG in withGap){
    seq <- seqs[[wG]][1]
    findr <- regexpr("A---A", seq) #should be at 427
    if(427 %in% findr[[1]]){
        newSeq <- gsub("A---A","---AA", seq)
        seqs[[wG]][1] <- newSeq
        countr = countr+1
    }
}
print(paste("replaced",countr,"type 2 gaps"))


#Gap type number 3
withGap <- grep("---G", seqs)

countr <- 0
for(wG in withGap){
    seq <- seqs[[wG]][1]
    findr <- gregexpr("---G", seq) #should be at 419
    if(419 %in% findr[[1]]){
        newSeq <- sub("A---G","---AG", seq)
        newSeq <- sub("G---G","---GG", newSeq)
        seqs[[wG]][1] <- newSeq
        countr = countr+1
    }
}
print(paste("replaced",countr,"type 3 gaps"))

#Gap type number 4
# This is more complex as need to move the gap futher!

withGap <- grep("[^-]---T", seqs) #only 3 gaps - not longer

countr <- 0
for(wG in withGap){
    seq <- seqs[[wG]][1]
    findr <- gregexpr("---T", seq) #should be at 420
    if(420 %in% findr[[1]]){
        #If the gap starts at 420, skip 3 bases (the gap),
        # and then take the next 6 bases and move these before the gap
        preGap <- substr(seq, 1, 419)
        gap <- substr(seq, 420, 422)
        postGap <- substr(seq, 423, nchar(seq))
        moveBases <- substr(postGap, 1, 7)
        restBases <- substr(postGap, 8, nchar(postGap))
        newSeq <- paste0(preGap, moveBases, gap, restBases)

        seqs[[wG]][1] <- newSeq
        countr = countr+1
    }
}
print(paste("replaced",countr,"type 4 gaps"))

write.fasta(seqs, seq_names, file.out=fixed_seqs, nbchar=10000)
