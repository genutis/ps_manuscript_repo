# setting up R env with custom package functions
install.packages("workflow/scripts/pspkg_1.0.tar.gz")
library(pspkg)

args <- commandArgs(trailingOnly = TRUE)

df_name=args[1]
cf_col=args[2]
features=args[3]
bins=args[4]
string="LTR"
te_vector=c("Ty1", "Ty2", "Ty3", "Ty4", "Ty5")

df_Y12 <- chrsplit(df, "Y12")
df_DBVPG6044 <- chrsplit(df, "DBVPG6044")


Y12_means_medians<-compute_means_medians(tes_bed, df_Y12)
DBVPG6044_means_medians<-compute_means_medians(tes_bed, df_DBVPG6044)


write.table(DBVPG6044_means_medians, file=args[5], quote = FALSE, row.names = FALSE, sep = "\t")
write.table(Y12_means_medians, file=args[6], quote = FALSE, row.names = FALSE, sep = "\t")

tes_granges<-makeGRangesFromDataFrame(tes_bed,  keep.extra.columns=TRUE)
DBVPG6044_granges<-makeGRangesFromDataFrame(df_DBVPG6044,  keep.extra.columns=TRUE)
Y12_granges<-makeGRangesFromDataFrame(df_Y12,  keep.extra.columns=TRUE)
test_gr <-DBVPG6044_granges
overlaps <- findOverlaps(query = tes_granges, subject = test_gr)
mcols(test_gr)$te_family <- "No TE"
mcols(test_gr)$te_family[subjectHits(overlaps)] <- mcols(tes_granges)$te_family[queryHits(overlaps)]

# get the medians now
test_gr_medians <- sapply(sort(unique(test_gr$te_family)), function (x) {median(test_gr[test_gr$"te_family"==x]$pairing_score, na.rm=TRUE)})
test_gr_means <- sapply(sort(unique(test_gr$te_family)), function (x) {mean(test_gr[test_gr$"te_family"==x]$pairing_score, na.rm=TRUE)})
