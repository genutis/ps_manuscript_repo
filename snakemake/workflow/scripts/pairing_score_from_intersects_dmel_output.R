# setting up R env with custom package functions
install.packages("workflow/scripts/pspkg_1.0.tar.gz")
library(pspkg)



#args <- commandArgs(trailingOnly = TRUE) 
args <- c("/scratch1/genutis/hic/tmp/beds/multicov.dm3.bed", "/scratch1/genutis/hic/figs", "dm3")
# path to bed input file with intersects
bedfile=args[1]
# output figure path and filename
figpath=args[2]
genome=args[3]
dir.create(file.path(figpath), showWarnings = FALSE)
figname="pairing_score_plots_"
lmname="pairing_score_mreg_lm_"


# true or false to run older RNA/DNA TE plotting analysis this run (for debugging and adding new functions.
oldplots <- TRUE
# make summary table or not
stable <- FALSE
#oldplots <- TRUE

# cutoff for chip read depth
chip_cutoff=100

# true or false to log transform certain columns to agree with log transformed ps numbers
#logtform <- TRUE
logtform <- FALSE
if(logtform == TRUE){
    figname <- paste0(figname, "log2_")
    lmname <- paste0(lmname, "log2_")
    chip_cutoff <- log2(chip_cutoff)
}

# vector of columns to apply transformation to
# 4:8 is dna and rna1-4 sequencing coverage
# 28 is number of snps
# 43, 51, 59, 67 are rna 1-4 te reads
# 80:120 are the chip binding coverage for different proteins from public data
# 97-136 is genes per bin, avg rna te coverage, dna and rna coverage, avg chip coverage, and dna te coverage.
tformcols <- c(4:8, 28, 59, 67, 75, 83, 97:136)

# string for pulling TEs
# this one only pulls new insertions events (in_repeat refers to events known in reference, pass is new ones)
#fstring="FILTER=PASS"
# for this string, it omits strand biased detected TEs (Reverse<forward^(0.8) or Reverse>forward^(1.25))
# and also TEs with weak evidence (less than 10 reads supporting)
fstring="FILTER=PASS|FILTER=in_repeat"

# importing data from bedtools pipelines and some wrangling:

data <- read.table(bedfile, sep = "\t")

# add column names
# note for chip names, mrg15 are rep 1 and rep2, fs1h-L are rep 1 and rep2. Last two GAFs are from same experiment, other GAFs are from different experiment
# last two drefs are from same experiment, others different experiment.
# beaf are all different experiments
# see /scratch/genutis/hic/data/chip/sra_gsm_table_chip.tsv table for more explanation
# patched name file
# this has been fixed again as of 2/14/25 to remove single ended chip seq samples
names(data) <- c("chrom", "bin_start", "bin_end", "chrom", "bin_start", "bin_end", "chrom", "bin_start", "bin_end", "dna_coverage", "rna1_coverage", "rna2_coverage", "rna3_coverage", "rna4_coverage", "ps_chrom", "ps_start", "ps_end", "pairing_score", "cis_score", "bp_overlap", "gtf_seqname", "gtf_source", "gtf_feature", "gtf_start", "gtf_stop", "gtf_score", "gtf_strand", "gtf_frame", "gtf_attribute", "gtf_bp_overlap", "pi_chrom", "pi_start", "pi_end", "snps", "pi", "pi_bp_overlap", "dna_te_chrom", "dna_te_start", "dna_te_end", "dna_te_family", "dna_te_reads", "dna_te_dot", "dna_te_quality", "dna_te_bp_overlap", "mat_te_chrom", "mat_te_start", "mat_te_end", "mat_te_family", "mat_te_reads", "mat_te_dot", "mat_te_quality", "mat_te_bp_overlap", "pat_te_chrom", "pat_te_start", "pat_te_end", "pat_te_family", "pat_te_reads", "pat_te_dot", "pat_te_quality", "pat_te_bp_overlap", "rna1_te_chrom", "rna1_te_start", "rna1_te_end", "rna1_te_family", "rna1_te_reads", "rna1_te_dot", "rna1_te_quality", "rna1_te_bp_overlap", "rna2_te_chrom", "rna2_te_start", "rna2_te_end", "rna2_te_family", "rna2_te_reads", "rna2_te_dot", "rna2_te_quality", "rna2_te_bp_overlap", "rna3_te_chrom", "rna3_te_start", "rna3_te_end", "rna3_te_family", "rna3_te_reads", "rna3_te_dot", "rna3_te_quality", "rna3_te_bp_overlap", "rna4_te_chrom", "rna4_te_start", "rna4_te_end", "rna4_te_family", "rna4_te_reads", "rna4_te_dot", "rna4_te_quality", "rna4_te_bp_overlap", "pct_at", "pct_gc", "num_A", "num_C", "num_G", "num_T", "num_N", "num_oth", "seq_len", "Mod(mdg4)2.2_SRR442320", "TFIIIC_SRR1636805", "CTCF_SRR1636769", "GAF_SRR1151106", "ZIPIC_SRR3452736", "CapH2_SRR1636752", "Rad21_SRR1636795", "GAF_SRR3452729", "DREF_SRR1636771", "CapH2_SRR1636753", "Nup98_SRR3452737", "Chromator_SRR1636761", "CP190_SRR1636766", "DREF_SRR1636770", "Fs1h-L_SRR1636773", "Fs1h-L_SRR1636774", "Rad21_SRR1636797", "Ibf2_SRR3452734", "CBP_SRR1636757", "Ibf1_SRR3452733", "RNAPII_SRR1636800", "Chromator_SRR1636762", "Rad21_SRR1636796", "BEAF_SRR1636749", "TFIIIC_SRR1636806", "Chromator_SRR1151105", "Pita_SRR3452735", "CBP_SRR1636756")


# removing positions that don't have pairing score information (shows up as ".")
data <- data[data$pairing_score != ".",]

# removing regions with zero dna reads
data <- data[data$dna_coverage != 0,]

# count genes per bin and put it in a new column
# pulls geneID from gtf column and counts uniques per rowrna3_te_quality
data$genes_per_bin <- unlist(lapply(data$gtf_attribute, function(x) {length(unique(regmatches(x, gregexpr("\\gene_id (.*?)\\;", x))[[1]]))}))

# average out the number of TE reads per bin with this line
data[grepl("rna._te_reads", names(data))] <- sapply(as.data.frame(sapply(data[grepl("rna._te_reads", names(data))], function(x) { gsub("\\.", 0, x) })), function (x) { sapply(strsplit((x),  split=","), function(y) mean(as.numeric(y))) } )
# this line removes the . values and splits strings and takes averages
#data$avg_rna_te_coverage <- rowMeans(sapply(as.data.frame(sapply(data[grepl("rna._te_reads", names(data))], function(x) { gsub("\\.", 0, x) })), function (x) { sapply(strsplit((x), split=","), function(y) mean(as.numeric(y))) } ))
# this line uses precomputed averages from above averaging of te reads 
data$avg_rna_te_coverage <- rowMeans(data[grepl("rna._te_reads", names(data))])

# collect average gene expression and chip coverage across four RNA replicates per bin and put it in a new column
data$avg_rna_coverage <- rowMeans(data[grepl("rna._coverage", names(data))])
data$avg_chip_coverage <- rowMeans(data[grepl("SRR", names(data))])
data$avg_dna_te_coverage <- rowMeans(sapply(as.data.frame(sapply(data[grepl("dna_te_reads", names(data))], function(x) { gsub("\\.", 0, x) })), function (x) { sapply(strsplit((x), split=","), function(y) mean(as.numeric(y))) } ))
data$avg_mat_te_coverage <- rowMeans(sapply(as.data.frame(sapply(data[grepl("mat_te_reads", names(data))], function(x) { gsub("\\.", 0, x) })), function (x) { sapply(strsplit((x), split=","), function(y) mean(as.numeric(y))) } ))
data$avg_pat_te_coverage <- rowMeans(sapply(as.data.frame(sapply(data[grepl("pat_te_reads", names(data))], function(x) { gsub("\\.", 0, x) })), function (x) { sapply(strsplit((x), split=","), function(y) mean(as.numeric(y))) } ))

# end data input section
# making table of summary statistics for pairing score and TEs
# note needs to be reworked with new format for data subsetting for plotting.
if(stable == TRUE){
    summarytable <- round(as.data.frame(list(mean=sapply(ls(pattern="te_"), function(x){ c(n_bins=nrow(get(x)), apply(na.omit(as.data.frame(lapply(get(x)[c(4,5,6,7,8,12,13,71)],as.numeric))),2,mean)) }), median=sapply(ls(pattern="te_"), function(x){ c(n_bins=nrow(get(x)), apply(na.omit(as.data.frame(lapply(get(x)[c(4,5,6,7,8,12,13,71)],as.numeric))),2,median)) })))[unlist(lapply(c(1,7:11,2:6), function(x) { seq(x, x+33, 11) }))], digits = 3)
    names(summarytable) <- gsub("\\.", "_", names(summarytable))

    # write to file
    write.csv(t(summarytable), file=paste0(stablepath,"/","summarytable_te_analysis.", genome, ".csv"))
}

# log transform colums specified in vector tformcols above in the script.
if(logtform == TRUE){
   print("Log2 transforming the following colums from the data:")
   print(names(data[tformcols]))
   data[, tformcols] <- log2(sapply(data[, tformcols], as.numeric)+1)
}

    df<-data.frame(V1 = rep(ps, sapply(s, length)), V2 = unlist(s))
    # convert ps to numeric values
    df2<-data.frame(as.numeric(df$V1), gsub("\\.", "No TE", df$V2))
    return(df2)
}

pdf(file = paste0(figpath, "/", figname, ".", genome, ".pdf"),


    width = 5, height = 5)


par(mar=c(7,4,4,2))


if(oldplots == TRUE){
    # old plots below here

# making boxplots of specific terms like pairing score
    for (string in c("pairing_score", "pct_gc")){

        # ps 5'utr containing 4kb bin vs not

        bp(
            data[grepl("CDS", data$gtf_feature) ,],
            data[!grepl("CDS", data$gtf_feature) ,],
            string,
            "PnM CDS-sequence windows",
            "With CDS", "Without CDS"
        ) 

        # boxplots for expressed in all replicate tes

        bp(
            data[grepl("CDS", data$gtf_feature)& grepl(fstring, data$rna1_te_quality) & grepl(fstring, data$rna2_te_quality) & grepl(fstring, data$rna3_te_quality) & grepl(fstring, data$rna4_te_quality) ,],
            data[!grepl(fstring, data$rna1_te_quality) & !grepl(fstring, data$rna2_te_quality) & !grepl(fstring, data$rna3_te_quality) & !grepl(fstring, data$rna4_te_quality) ,],
             string,
             "PnM all RNA replicates CDS",
             "With TEs", "Without TEs"
        )   
        # boxplots for expressed in all replicate tes and dna

        bp(
            data[grepl("CDS", data$gtf_feature) & grepl(fstring, data$dna_te_quality) & grepl(fstring, data$rna1_te_quality) & grepl(fstring, data$rna2_te_quality) & grepl(fstring, data$rna3_te_quality) & grepl(fstring, data$rna4_te_quality) ,],
            data[!grepl(fstring, data$dna_te_quality) & !grepl(fstring, data$rna1_te_quality) & !grepl(fstring, data$rna2_te_quality) & !grepl(fstring, data$rna3_te_quality) & !grepl(fstring, data$rna4_te_quality) ,],
            string,
            "PnM DNA and all RNA replicates CDS",
            "With TEs", "Without TEs"
        )     
        # ps 5'utr regions with and without tes
        sapply(
#            names(data[grepl("_te_quality", names(data))])[-c(1:3)],
            names(data[grepl("_te_quality", names(data))]),
            function(name){
                bp(
                    data[grepl("CDS", data$gtf_feature) & grepl(fstring, data[[name]]) ,],
                    data[grepl("CDS", data$gtf_feature) & !grepl(fstring, data[[name]]) ,],
                    string,
                    paste("PnM", toupper(gsub("_te_quality", "", name)), "CDS-sequence windows"),
                    "With TEs", "Without TEs"
                )   
            }
        )

    # ps 5'utr containing 4kb bin vs not

        bp(
            data[grepl("5UTR", data$gtf_feature) ,],
            data[!grepl("5UTR", data$gtf_feature) ,],
            string,
            "PnM 5'-UTR windows",
            "With 5'-UTR", "Without 5'-UTR"
        ) 

        # boxplots for expressed in all replicate tes

        bp(
            data[grepl("5UTR", data$gtf_feature)& grepl(fstring, data$rna1_te_quality) & grepl(fstring, data$rna2_te_quality) & grepl(fstring, data$rna3_te_quality) & grepl(fstring, data$rna4_te_quality) ,],
            data[!grepl(fstring, data$rna1_te_quality) & !grepl(fstring, data$rna2_te_quality) & !grepl(fstring, data$rna3_te_quality) & !grepl(fstring, data$rna4_te_quality) ,],
             string,
             "PnM all RNA replicates 5'-UTR",
             "With TEs", "Without TEs"
        )   
        # boxplots for expressed in all replicate tes and dna

        bp(
            data[grepl("5UTR", data$gtf_feature) & grepl(fstring, data$dna_te_quality) & grepl(fstring, data$rna1_te_quality) & grepl(fstring, data$rna2_te_quality) & grepl(fstring, data$rna3_te_quality) & grepl(fstring, data$rna4_te_quality) ,],
            data[!grepl(fstring, data$dna_te_quality) & !grepl(fstring, data$rna1_te_quality) & !grepl(fstring, data$rna2_te_quality) & !grepl(fstring, data$rna3_te_quality) & !grepl(fstring, data$rna4_te_quality) ,],
            string,
            "PnM DNA and all RNA replicates 5'-UTR",
            "With TEs", "Without TEs"
        )     
        # ps 5'utr regions with and without tes
        sapply(
#            names(data[grepl("_te_quality", names(data))])[-1],
            names(data[grepl("_te_quality", names(data))]),
            function(name){
                bp(
                    data[grepl("5UTR", data$gtf_feature) & grepl(fstring, data[[name]]) ,],
                    data[grepl("5UTR", data$gtf_feature) & !grepl(fstring, data[[name]]) ,],
                    string,
                    paste("PnM", toupper(gsub("_te_quality", "", name)), "5'-UTR windows"),
                    "With TEs", "Without TEs"
                )   
            }
        )


        # boxplot for expressed and dna tes
        sapply(
#            names(data[grepl("_te_quality", names(data))])[-1],
            names(data[grepl("_te_quality", names(data))])[-c(1:3)],
            function(name){
                bp(
                    data[grepl(fstring, data$dna_te_quality) & grepl(fstring, data[[name]]) ,],
                    data[!grepl(fstring, data$dna_te_quality) & !grepl(fstring, data[[name]]) ,],
                    string,
                    paste("PnM DNA and", toupper(gsub("_te_quality", "", name)), "Intersect"),
                    "With TEs", "Without TEs"
                )   
            }
        )
        # boxplot for expressed tes
        sapply(
            names(data[grepl("_te_quality", names(data))])[-c(1:3)],
            function(name){
                bp(
                    data[grepl(fstring, data[[name]]) ,],
                    data[!grepl(fstring, data[[name]]) ,],
                    string,
                    paste("PnM", toupper(gsub("_te_quality", "", name))),
                    "With TEs", "Without TEs"
                )   
            }
        )
        # boxplots for expressed in all replicate tes

        bp(
            data[grepl(fstring, data$rna1_te_quality) & grepl(fstring, data$rna2_te_quality) & grepl(fstring, data$rna3_te_quality) & grepl(fstring, data$rna4_te_quality) ,],
            data[!grepl(fstring, data$rna1_te_quality) & !grepl(fstring, data$rna2_te_quality) & !grepl(fstring, data$rna3_te_quality) & !grepl(fstring, data$rna4_te_quality) ,],
             string,
             "PnM all RNA replicates",
             "With TEs", "Without TEs"
        )   
        # boxplot of detecte tes in dna
        bp(
            data[grepl(fstring, data$dna_te_quality) ,],
            data[!grepl(fstring, data$dna_te_quality) ,],
            string,
            "PnM DNA",
            "With TEs", "Without TEs"
        )   
        # boxplot of detecte tes in mat dna
        bp(
            data[grepl(fstring, data$mat_te_quality) ,],
            data[!grepl(fstring, data$mat_te_quality) ,],
            string,
            "PnM maternal DNA",
            "With TEs", "Without TEs"
        )   
        # boxplot of detecte tes in dna
        bp(
            data[grepl(fstring, data$pat_te_quality) ,],
            data[!grepl(fstring, data$pat_te_quality) ,],
            string,
            "PnM paternal DNA",
            "With TEs", "Without TEs"
        )   

        # boxplots for expressed in all replicate tes and dna

        bp(
            data[grepl(fstring, data$dna_te_quality) & grepl(fstring, data$rna1_te_quality) & grepl(fstring, data$rna2_te_quality) & grepl(fstring, data$rna3_te_quality) & grepl(fstring, data$rna4_te_quality) ,],
            data[!grepl(fstring, data$dna_te_quality) & !grepl(fstring, data$rna1_te_quality) & !grepl(fstring, data$rna2_te_quality) & !grepl(fstring, data$rna3_te_quality) & !grepl(fstring, data$rna4_te_quality) ,],
            string,
            "PnM DNA and all RNA replicates",
            "With TEs", "Without TEs"
        )   

        
        # boxplot for pairing score for chip binding regions of individual proteins using >99 reads per 4kb window just as a cutoff ( i was selecting a lot of windows with this at zero)
            sapply(
                names(data[grepl("SRR", names(data))]),
                function(name){
                    bp(
                        data[data[[name]] >=chip_cutoff,],
                        data[data[[name]] <chip_cutoff,],
                        string,
                        paste("ChIP-seq", name),
                        paste("With", gsub("_SRR.*", "", name)), paste("Without", gsub("_SRR.*", "", name))
                    )
                }
            )

    }
    # pairing score plotted against avg rna coverage per bin
    sp(data$avg_rna_coverage, data$pairing_score, "Pairing Score and RNA Coverage", "Number of RNA reads per 4kb window", "Pairing Score per 4kb window",10)

    # pairing score plotted against genes per bin
    sp(data$genes_per_bin, data$pairing_score, "Pairing Score and number of Genes", "Number of Genes within each 4kb window", "Pairing Score per 4kb window",3)


    # pairing score plotted against snps
    sp(data$snps, data$pairing_score, "Pairing Score and number of SNPs", "Number of SNPs per 4kb window", "Pairing Score per 4kb window",0)

    # pairing score plotted against gc content
    sp(data$pct_gc, data$pairing_score, "Pairing Score and GC-content", "GC-content % per 4kb window", "Pairing Score per 4kb window",0)

    # pairing score and te reads
    sp(data$avg_rna_te_coverage, data$pairing_score, "Pairing Score and TE Coverage", "TE-specific RNA reads per 4kb window", "Pairing Score per 4kb window",0)
    sp(data$avg_dna_te_coverage, data$pairing_score, "Pairing Score and TE Coverage", "TE-specific DNA reads per 4kb window", "Pairing Score per 4kb window",0)
    sp(data$avg_mat_te_coverage, data$pairing_score, "Pairing Score and TE Coverage", "TE-specific maternal DNA reads per 4kb window", "Pairing Score per 4kb window",0)
    sp(data$avg_pat_te_coverage, data$pairing_score, "Pairing Score and TE Coverage", "TE-specific paternal DNA reads per 4kb window", "Pairing Score per 4kb window",0)

    # ChIP plotting
    # scatter plot for degree of protein binding and pairing score for each protien
   sapply(
            names(data[grepl("SRR", names(data))]),
            function(name){
                sp(data[[name]], 
                    data$pairing_score, 
                    paste("Pairing Score and", gsub("_SRR.*", "", name), "ChIP-seq"), 
                    paste("Number of", gsub("_SRR.*", "", name), "ChIP-seq reads per 4kb window"), 
                    "Pairing Score per 4kb window", 5)
            }
        )
    # scatter plot of chip binding vs not for all proteins on chip
    sp(data$avg_chip_coverage, data$pairing_score, "Pairing Score and ChIP-seq (Multiple Proteins)", "Number of ChIP-seq reads per 4kb window", "Pairing Score per 4kb window", 5)



}

# new plots below here


# boxplot of tes stuff below here
tes<-unsplit_te(data$pairing_score, data$dna_te_family)
names(tes)<-c("pairing_score", "dna_te_family")

testing<-tetester(tes, "pairing_score", "dna_te_family", 0.05, "t", "bonferroni")
# subset based on being in this p value cutoff of 0.05 for plotting
plote<-as.data.frame(tes[tes$dna_te_family %in% c(names(testing), "No TE"),])
# section to count groups
levs<-levels(factor(plote$dna_te_family))
# first before messing with names i will get my p values using names
pvals<-as.vector(signif(testing[order(match(names(testing),levs))],3))
# then replace levs
levs2<-paste0(levs, " (n=", as.vector(table(plote$dna_te_family)), ")")
plote$dna_te_family <- factor(plote$dna_te_family, labels = levs2)
#bymedian <- with(plote, reorder("dna_te_family", "pairing_score", median, na.rm=TRUE))
bymedian <- with(plote, reorder(dna_te_family, pairing_score, median, na.rm=T))
# bigger plot for horizontal figs
pdf(file = paste0(figpath, "/", figname, "horizontal.", genome, ".pdf"),
    width = 11, height = 8) 

par(las=2) # perpendicular axis labels with las=2
par(mar=c(10,8,4,1))


# plot of het or homozygous tes


sp1<- sp(as.numeric(data[grepl(fstring, data[["dna_te_quality"]]),][["dna_te_bp_overlap"]]),
        as.numeric(data[grepl(fstring, data[["dna_te_quality"]]),][["pairing_score"]]), # this will be pairing score
        main="Pairing Score and number of Base Pairs\noverlapping with TE",
        xlab="Base Pairs of overlap with TE sequence",
        ylab="Pairing Score",
        outliers=0)



# transcript id is v13, family id is v15, class id is v19

data<-get_het_hom_te(data, "dna_te_family", "mat_te_family", "pat_te_family", "het_te_family", "het_mat_te_family", "het_pat_te_family", "hom_te_family", "undetected_te_family")# something to iterate over

# gather tes
te_list<-list(
        as.numeric(data[!grepl("\\.", data[["het_te_family"]]) & grepl(fstring, data[["dna_te_quality"]]),][["pairing_score"]]),
        as.numeric(data[!grepl("\\.", data[["het_mat_te_family"]]) & grepl(fstring, data[["dna_te_quality"]]),][["pairing_score"]]),
        as.numeric(data[!grepl("\\.", data[["het_pat_te_family"]]) & grepl(fstring, data[["dna_te_quality"]]),][["pairing_score"]]),
        as.numeric(data[!grepl("\\.", data[["hom_te_family"]]) & grepl(fstring, data[["dna_te_quality"]]),][["pairing_score"]]),
        as.numeric(data[grepl("\\.", data[["dna_te_family"]]),][["pairing_score"]])
        )

hetnames <- c("Heterozygous TEs", "Maternal Het. TEs", "Paternal Het. TEs", "Homozygous TEs", "No TEs")

    a<-boxplot(te_list,
        main = "Heterozygous TEs",
        ylim=c(-5,2),
        names = hetnames,
        ylab = "pairing score per 4kb window")
        axis(side = 1, at = seq_along(a$names), labels = a$names, tick = FALSE) # sets axis to be horizontal

make_het_te_summary <- function(data,names) {
    tmp <-data
    tmp[["het_te_family"]]<-replace(tmp[["het_te_family"]], which(tmp[["het_te_family"]] != "."), names[1])
    tmp[["het_mat_te_family"]]<-replace(tmp[["het_mat_te_family"]], which(tmp[["het_mat_te_family"]] != "."), names[2])
    tmp[["het_pat_te_family"]]<-replace(tmp[["het_pat_te_family"]], which(tmp[["het_pat_te_family"]] != "."), names[3])
    tmp[["hom_te_family"]]<-replace(tmp[["hom_te_family"]], which(tmp[["hom_te_family"]] != "."), names[4])
    tmp$te_summary<-"No TE"
    tmp$te_summary<-replace(tmp$te_summary, which(tmp[["het_te_family"]] != "."), tmp[["het_te_family"]]) 
    tmp$te_summary<-replace(tmp$te_summary, which(tmp[["het_mat_te_family"]] != "."), tmp[["het_mat_te_family"]]) 
    tmp$te_summary<-replace(tmp$te_summary, which(tmp[["het_pat_te_family"]] != "."), tmp[["het_pat_te_family"]]) 
    tmp$te_summary<-replace(tmp$te_summary, which(tmp[["hom_te_family"]] != "."), tmp[["hom_te_family"]]) 
    data$te_summary<-tmp$te_summary
    
    return(data)
}

  tes<-unsplit_te(data$pairing_score, data$dna_te_family)
    names(tes)<-c("pairing_score", "dna_te_family")
testing<-tetester(tes, "pairing_score", "dna_te_family", 0.05, "t", "bonferroni")

tegtf<-read.table("resources/genomes/dm3_rmsk_TE.gtf")

teclasses<-sapply(names(testing), function(name) { unique(tegtf[grepl(name, tegtf[[13]]),][c(16,19)])})
teclasses_withps<-rbind(testing, as.data.frame(teclasses))
teclasses_withps[["I_DM"]][2] <-"I" # these need to be fixed by hand
teclasses_withps[["I_DM"]][3] <-"LINE"
teclasses_withps <- apply(teclasses_withps,2,as.character)
teclasses_withps<-t(teclasses_withps)
colnames(teclasses_withps) <-c("P. Value", "Family", "Class")
write.csv(as.data.frame(teclasses_withps), file=paste0(figpath, "/enriched_te_table.", genome, ".csv"))


# plotting tes
b<-boxplot(pairing_score ~ bymedian, plote, main = "Pairing Score per Transposable Element (TE) Family in PnM DNA", ylab = "Pairing Score per 4kb window", xlab = "")

text(1:length(b$n), b$stats[5,]+1, paste("p=", c(pvals,"")), srt=90, cex=0.8)


fit <- lm(as.data.frame(na.omit(sapply(data[,c(97:136,138:142)], as.numeric))))
summary(fit)
plot(fit)
# save lm results to file
j <- 0
while (file.exists(paste0(figpath, "/", lmname, sprintf("%03d", j), ".txt"))) { j <- j+1 }
sink(paste0(figpath, "/", lmname, sprintf("%03d", j), ".txt"))
sink(paste0(figpath, "/", lmname, ".", genome, ".txt"))
print(summary(fit))
sink()

dev.off()

