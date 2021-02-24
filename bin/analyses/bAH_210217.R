islets <- read.table('~/data/21013-01/nf_ATACseq_210219/bwa/mergedReplicate/macs/broadPeak/islet.mRp.clN_peaks.annotatePeaks.txt', header = T, sep = '\t')
names(islets)[1] <- 'PeakID'
pancreas <- read.table('~/data/21013-01/nf_ATACseq_210219/bwa/mergedReplicate/macs/broadPeak/pancreas.mRp.clN_peaks.annotatePeaks.txt', header = T, sep = '\t')
names(pancreas)[1] <- 'PeakID'

islets$feature <- unlist(sapply(islets$Annotation, function(x) {unlist(strsplit(x, ' '))[1]}))
pancreas$feature <- unlist(sapply(pancreas$Annotation, function(x) {unlist(strsplit(x, ' '))[1]}))

islets <- islets[islets$feature == 'promoter-TSS',]
islets <- islets[!is.na(islets$PeakID),]
pancreas <- pancreas[pancreas$feature == 'promoter-TSS',]
pancreas <- pancreas[!is.na(pancreas$PeakID),]

int <- intersect(pancreas$Annotation, islets$Annotation)
islets_only <- setdiff(islets$Annotation, pancreas$Annotation)

islets_only <- islets[islets$Annotation %in% islets_only,]

score_diff <- as.data.frame(int)
names(score_diff) <- 'Annotation'
score_diff$score_diff <- NA
score_diff$chr <- NA
score_diff$start <- NA
score_diff$end <- NA
score_diff$name <- NA
for (i in 1:length(score_diff$Annotation)){
  ann <- score_diff$Annotation[i]
  tmp_islet <- islets[islets$Annotation == ann,]
  tmp_pancreas <- pancreas[pancreas$Annotation == ann, ]
  
  if (length(tmp_islet[,1]) == 1 & length(tmp_pancreas[,1]) == 1){
    score_diff$score_diff[i] <- tmp_islet$Peak.Score - tmp_pancreas$Peak.Score
  }
  else{
    i_avg <- mean(tmp_islet$Peak.Score)
    p_avg <- mean(tmp_islet$Peak.Score)
    score_diff$score_diff[i] <- i_avg - p_avg
  }
  
  score_diff$chr[i] <- tmp_islet$Chr[1]
  score_diff$name[i] <- tmp_islet$Gene.Name[1]
  score_diff$start[i] <- min(c(tmp_islet$Start, tmp_pancreas$Start))
  score_diff$end[i] <- max(c(tmp_islet$Start, tmp_pancreas$Start))

}

score_diff <- score_diff[order(score_diff$score_diff, decreasing = TRUE),]


vs <- read.table('~/data/21013-01/nf_ATACseq_210219/bwa/mergedLibrary/macs/broadPeak/consensus/deseq2/isletvspancreas/isletvspancreas.mLb.clN.deseq2.results.txt', header = T)
vs_filt <- vs[vs$padj < 0.05,]
vs_filt <- drop_na(vs_filt)
vs_filt <- vs[vs$log2FoldChange > 0,]
