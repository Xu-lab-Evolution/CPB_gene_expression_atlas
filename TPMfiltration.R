#Script filtration by TPM

Counts.newanno <- read.table("D:/path/to/the/file/<table of raw counts obtained from featureCounts>.txt", sep="\t", stringsAsFactors=FALSE,
                             header=TRUE)

#### Removing lowly expressed genes ----

ccols <- c(6:66) # sample1_column:sampleN_column
lcol <- 6 # Column containing length of sum of exons per gene (Please place it in last column to avoid confusion)

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpms <- apply(Counts.newanno[,7:67], 2, function(x) tpm(x, Counts.newanno[,lcol]))
row.names(tpms) <- Counts.newanno[,1]

#filtration
Counts.tpm.flt <- as.data.frame(tpms[rowSums(tpms>1)>=2,])

Geneid <- rownames(Counts.tpm.flt)
rownames(Counts.tpm.flt) <- NULL
Counts.tpm.flt <- cbind(Geneid,Counts.tpm.flt)

Counts.flt <- semi_join(Counts.newanno, Counts.tpm.flt, by="Geneid")
write.table(Counts.flt, file = "D:/path/to/the/new/table/<table of filtered counts>.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# After filtration we keep 15578 genes out of 34350
