files <- list.files("~/Downloads/GSE110372_RAW", full.names = T)

for (file in files) {
  if(file == files[1]){
    name <- strsplit(file, "_")[[1]][3]
    print(name)
    counts <- read.table(file, col.names = name, row.names = 1)
  } else {
    name <- strsplit(file, "_")[[1]][3]
    print(name)
    tmp <- read.table(file, col.names = name, row.names = 1)
    counts <- cbind(counts, tmp)
  }
}

colnames(counts) <- gsub("counts24C", "24hrV1_", colnames(counts))
colnames(counts) <- gsub("countsC", "6hrV1_", colnames(counts))
colnames(counts) <- gsub("counts24T3", "24T3_", colnames(counts))
colnames(counts) <- gsub("countsT3", "6hrT3_", colnames(counts))
colnames(counts) <- gsub("countsCHX[A-D]", "6hrCHXV1_", colnames(counts))
colnames(counts) <- gsub("countsCHX", "6hrCHXT3_", colnames(counts))
colnames(counts)

write.csv(counts,"~/Downloads/GSE110372_RAW/Total_counts.csv", row.names = TRUE)
