#quality scores from Illumina pairs

base_scores <- read.table("..//Dunk_Illumina_mate_pairs_updated_Illumina_qualities.qual_per_base", sep="\t", fill=TRUE)
names(base_scores) <- c(1:76)
pdf("Illumina_dunk_100000_subset_quality_by_base.pdf")

summary_all <- c()
	
for(base_pos in 1:dim(base_scores)[2])
{
	summary_all <- cbind(summary_all, summary(base_scores[,base_pos]))
}

cex_sizes <- rep(1, 6) # c(2, 1.7, 1.4, 1.2,1)
shapes <- rep(16,6)
	
plot(1:dim(base_scores)[2], summary_all[1,], xlab="Base position", ylab="Quality Score (Sanger)", main="Illumina quality scores (by base position) (Dunk)", xlim=c(1,101), ylim=c(0,40),pch=shapes[1], cex=cex_sizes[1])
	
for(i in 2:6)
{
	if(i != 4)
	{
		points(1:dim(base_scores)[2], summary_all[i,], pch=shapes[i], col=i, cex=cex_sizes[i])
	}
}
	
legend("topright", legend=c("Minimum", "Q2", "Median","Q3", "Maximum"), col=c(1,2,3,5,6), bty="n", horiz=TRUE, pch=16)
dev.off()

#quality scores from Illumina pairs (Cleaned)
base_scores <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_mate_pairs.qual_per_base", sep="\t", fill=TRUE)
names(base_scores) <- c(1:76)
pdf("Illumina_dunk_100000_trimmed_mates_subset_quality_by_base.pdf")

summary_all <- c()

for(base_pos in 1:dim(base_scores)[2])
{
        summary_all <- cbind(summary_all, summary(base_scores[,base_pos]))
}

cex_sizes <- rep(1, 6) # c(2, 1.7, 1.4, 1.2,1)
shapes <- rep(16,6)

plot(1:dim(base_scores)[2], summary_all[1,], xlab="Base position", ylab="Quality Score (Phred)", main="Illumina quality scores (by base position) (Dunk - trimmed mates)", xlim=c(1,101), ylim=c(0,40),pch=shapes[1], cex=cex_sizes[1])

for(i in 2:6)
{
        if(i != 4)
        {
                points(1:dim(base_scores)[2], summary_all[i,], pch=shapes[i], col=i, cex=cex_sizes[i])
        }
}

legend("topright", legend=c("Minimum", "Q2", "Median","Q3", "Maximum"), col=c(1,2,3,5,6), bty="n", horiz=TRUE, pch=16)
dev.off()



#quality_scores from Illumina singletons (Cleaned)

base_scores <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_solitary_reads.qual_per_base", sep="\t", fill=TRUE)
names(base_scores) <- c(1:76)
pdf("Illumina_dunk_100000_singleton_trimmed_subset_quality_by_base.pdf")
summary_all <- c()

for(base_pos in 1:dim(base_scores)[2])
{
        summary_all <- cbind(summary_all, summary(base_scores[,base_pos]))
}

cex_sizes <- rep(1, 6) # c(2, 1.7, 1.4, 1.2,1)
shapes <- rep(16,6)

plot(1:dim(base_scores)[2], summary_all[1,], xlab="Base position", ylab="Quality Score (Phred)", main="Illumina quality scores (by base position) (Dunk - trimmed singletons)", xlim=c(1,101), ylim=c(0,40),pch=shapes[1], cex=cex_sizes[1])

for(i in 2:6)
{
        if(i != 4)
        {
                points(1:dim(base_scores)[2], summary_all[i,], pch=shapes[i], col=i, cex=cex_sizes[i])
        }
}

legend("topright", legend=c("Minimum", "Q2", "Median","Q3", "Maximum"), col=c(1,2,3,5,6), bty="n", horiz=TRUE, pch=16)
dev.off()

