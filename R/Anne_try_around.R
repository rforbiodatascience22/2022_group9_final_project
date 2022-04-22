
# checks data
dim(proteomes_raw)
view(proteomes_raw)
view(PAM50_raw)
view(patients_raw)

proteomes_raw$RefSeq_accession_number

length(unique(PAM50_raw$RefSeqProteinID))
length(unique(proteomes_raw$RefSeq_accession_number))
length(PAM50_raw$RefSeqProteinID)
