# setting up R env with custom package functions
install.packages("workflow/scripts/pspkg_1.0.tar.gz")
library(pspkg)

args <- commandArgs(trailingOnly = TRUE)


df_name_du_mat=args[1]
df_name_du_pat=args[2]
cf_col_du=args[3]
bins_du_mat=args[4]
bins_du_pat=args[5]
tes_du=args[6]


print(paste(Sys.time(), "Step 2/19: reading in Du et. al 2017 maternal mouse data"))
if(!exists("dat_du_mat")){   dat_du_mat <- read_and_avg_df(df_name_du_mat, cf_col_du)   }
print(paste(Sys.time(), "Step 3/19: reading in Du et al. 2017 paternal mouse data"))
if(!exists("dat_du_pat")){    dat_du_pat <- read_and_avg_df(df_name_du_pat, cf_col_du)   }
print(paste(Sys.time(), "Step 5/19: combining bins with chromosomal coordinates for mouse maternal data"))
if(!exists("df_du_mat")){
    df_du_mat <- combine_bins_and_df(dat_du_mat, bins_du_mat)   
    df_du_mat[[1]]<-gsub("^.{0,3}", "", df_du_mat[[1]]) # this fixes the chr1 to 1 in the first column
}
print(paste(Sys.time(), "Step 6/19: combining bins with chromosomal coordinates for mouse paternal data"))
if(!exists("df_du_pat")){
    df_du_pat <- combine_bins_and_df(dat_du_pat, bins_du_pat)
    df_du_pat[[1]]<-gsub("^.{0,3}", "", df_du_pat[[1]]) # this fixes the chr1 to 1 in the first column
}
print(paste(Sys.time(), "Step 10/19: parsing fasta for transposon coordinates for mm9 mouse genome"))
if(!exists("tes_bed_du")){    tes_bed_du <-lookup_tes_mm9(tes_du)   }
print(paste(Sys.time(), "Step 11/19: counting TEs for mouse maternal data"))
if(!exists("te_counts_du_mat")){    te_counts_du_mat <- te_counter(tes_bed_du, df_du_mat)    }  
print(paste(Sys.time(), "Step 12/19: counting TEs for mouse paternal data"))
if(!exists("te_counts_du_pat")){    te_counts_du_pat <- te_counter(tes_bed_du, df_du_pat)    }  
print(paste(Sys.time(), "Step 17/19: computing mean and median pairing scores for TEs and total genomic regions for mouse maternal data"))
if(!exists("du_mat_means_medians")){    du_mat_means_medians<-compute_means_medians(tes_bed_du, df_du_mat)    }   
print(paste(Sys.time(), "Step 18/19: computing mean and median pairing scores for TEs and total genomic regions for mouse paternal data"))
if(!exists("du_pat_means_medians")){    du_pat_means_medians<-compute_means_medians(tes_bed_du, df_du_pat)    }  

write.table(du_mat_means_medians, file=args[7], quote = FALSE, row.names = FALSE, sep = "\t")
write.table(du_pat_means_medians, file=args[8], quote = FALSE, row.names = FALSE, sep = "\t")
write.table(te_counts_du_mat, file=args[9], quote = FALSE, row.names = FALSE, sep = "\t")
write.table(te_counts_du_pat, file=args[10], quote = FALSE, row.names = FALSE, sep = "\t")

