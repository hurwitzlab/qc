#analysis of preliminary GBR data

#subset is 500,000 randomly sampled lines.

#something odd here...

newScheme <- heat.colors(17)


require(MASS)
#Illumina - paired uncleaned

#try to maintain the same bandwidth
#library(k2e2d)
	dunkIlluminaPaired <- read.table("../Dunk_Illumina_mate_pairs_updated_Illumina_qualities.size_and_gc.subset", header=FALSE, sep="\t")
	dim(dunkIlluminaPaired)
	dunkIlluminaPaired[1,]
	
	lengths <- dunkIlluminaPaired[,2]
	gc <- dunkIlluminaPaired[,3]

	overall_bandwidth <- c(0.1, 2)
	
	#took the h bandwidth out
	f1 <- kde2d(gc, lengths, h=overall_bandwidth,  n = 50, lims = c(0, 1, 0, 120))

	pdf("joint_distribution_gc_and_read_length_illumina_dunk.pdf", width=10, height=10)
	filled.contour(f1, col = rev(newScheme), nlevels = 11, xlab="G+C content", ylab="Read length (bp)", main="Joint distribution of G+C and read length (Dunk Illumina)")
	dev.off()

	pdf("hist_of_read_lengths_mate_pairs_raw.pdf")
	hist(lengths, main="Histogram of read-lengths (Dunk Illumina - raw)", breaks=101)
	dev.off()
		
	pdf("kernel_density_of_GC_random_100K_subset_dunk_Illumina.pdf")
	subset <- sample(gc, 100000, replace = FALSE, prob = NULL)
	mymean <- mean(subset)
	hist(subset, breaks=200)
	d <- density(subset)
	plot(d, main="Density plot of G+C content (Dunk Illumina)")
	polygon(d, col="red", border="black") 
	abline(v=mymean, col="blue", lty=2, lwd=2)
	text(x=mymean + .15, y=3, labels=c(paste("Mean", round(mymean,3))))
	dev.off()

	pdf("kernel_density_of_read_lengths_random_100K_subset_dunkIllumina.pdf")
	subset_l <- sample(lengths, 100000, replace = FALSE, prob = NULL)
	d2 <- density(subset_l,bw=0.5)
	mymean <- mean(subset_l)
	boxplot(subset_l,breaks=100)
	plot(d2, main="Density,  plot of read lengths (Dunk Illumina)")
	polygon(d2, col="red", border="black") 
	abline(v=mymean, col="blue", lty=2, lwd=2)
	text(x=mymean + 15, y=0.03, labels=c(paste("Mean", round(mymean,3))))
	dev.off()
	
	#length and GC distributions are done for both datasets now.
	#clear the decks for quality scores

	rm(dunkIlluminaPaired)
	rm(subset)
	rm(subset_l)
	rm(subset_ln)


#cleaned pairs
	dunkIlluminaPairedCleaned <-read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_mate_pairs.size_and_gc.subset",sep="\t", header=FALSE);
        dim(dunkIlluminaPairedCleaned)
        dunkIlluminaPairedCleaned[1,]

        lengths <- dunkIlluminaPairedCleaned[,2]
        gc <- dunkIlluminaPairedCleaned[,3]


        #took the h bandwidth out
        f1 <- kde2d(gc, lengths, h=overall_bandwidth, n = 50, lims = c(0, 1, 0, 120))

        pdf("joint_distribution_gc_and_read_length_illumina_cleaned_pairs_no_length_cutoff_dunk.pdf", width=10, height=10)
        filled.contour(f1, col = rev(newScheme), nlevels = 11, xlab="G+C content", ylab="Read length (bp)", main="Joint distribution of G+C and read length (Dunk Illumina - cleaned)")
        dev.off()
        
	pdf("hist_of_read_lengths_mate_pairs_cleaned.pdf")
        hist(lengths, main="Histogram of read-lengths (Dunk Illumina - cleaned MP)", breaks=101)
        dev.off()

	pdf("kernel_density_of_GC_random_100K_subset_cleaned_pairs_dunk_Illumina_no_length_cutoff.pdf")
	subset <- sample(gc, 100000, replace = FALSE, prob = NULL)
	mymean <- mean(subset)
	hist(subset, breaks=200)
	d <- density(subset)
	plot(d, main="Density plot of G+C content (Dunk Illumina)")
	polygon(d, col="red", border="black")
	abline(v=mymean, col="blue", lty=2, lwd=2)
	text(x=mymean + .15, y=3, labels=c(paste("Mean", round(mymean,3))))
	dev.off()

        pdf("kernel_density_of_read_lengths_random_100K_subset_cleaned_pairs_dunk_Illumina_no_length_cutoff.pdf")
        subset_l <- sample(lengths, 100000, replace = FALSE, prob = NULL)
        d2 <- density(subset_l,bw=0.5)
       	mymean <- mean(subset_l)
       	boxplot(subset_l,breaks=100, main="Histogram of read lengths (Dunk Illumina - cleaned)")
       	plot(d2, main="Density plot of read lengths (Dunk Illumina - cleaned)")
       	polygon(d2, col="red", border="black")
       	abline(v=mymean, col="blue", lty=2, lwd=2)
       	text(x=mymean + 15, y=0.03, labels=c(paste("Mean", round(mymean,3))))
       	dev.off()

        #length and GC distributions are done for both datasets now.
        #clear the decks for quality scores

        rm(dunkIllumina)
        rm(subset)
        rm(subset_l)
        rm(subset_ln)



	#singletons

       	dunkIlluminaSingle <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_solitary_reads.size_and_gc.subset", header=FALSE, sep="\t");
        dim(dunkIlluminaSingle)
        dunkIlluminaSingle[1,]

        lengths <- dunkIlluminaSingle[,2]
        gc <- dunkIlluminaSingle[,3]

       	 #took the h bandwidth out
       	 f1 <- kde2d(gc, lengths, n = 50, h=overall_bandwidth, lims = c(0, 1, 0, 120))
        
	pdf("joint_distribution_gc_and_read_length_illumina_cleaned_singletons_dunk_no_length_cutoff.pdf", width=10, height=10)
       	filled.contour(f1, col = rev(newScheme), nlevels = 11, xlab="G+C content", ylab="Read length (bp)", main="Joint distribution of G+C and read length (Dunk Illumina - cleaned singletons)")
       	dev.off()

       	pdf("kernel_density_of_GC_random_100K_subset_cleaned_singletons_dunk_Illumina_no_length_cutoff.pdf")
       	subset <- sample(gc, 100000, replace = FALSE, prob = NULL)
       	mymean <- mean(subset)
       	hist(subset, breaks=200, main="Histogram of singletons (Dunk Illumina - cleaned)")
       	d <- density(subset)
       	plot(d, main="Density plot of G+C content (Dunk Illumina) - cleaned singletons")
       	polygon(d, col="red", border="black")
       	abline(v=mymean, col="blue", lty=2, lwd=2)
       	text(x=mymean + .15, y=3, labels=c(paste("Mean", round(mymean,3))))
       	dev.off()
                
	pdf("hist_of_read_lengths_mate_pairs_singletons_no_length_cutoff.pdf")
        hist(lengths, main="Histogram of read-lengths (Dunk Illumina - trimmed singletons)", breaks=101)
        dev.off()
	        
		pdf("kernel_density_of_read_lengths_random_100K_subset_cleaned_singletons_dunk_Illumina_no_length_cutoff.pdf")
	        subset_l <- sample(lengths, 100000, replace = FALSE, prob = NULL)
	        d2 <- density(subset_l,bw=0.5)
	        mymean <- mean(subset_l)
	        boxplot(subset_l,breaks=100)
	        plot(d2, main="Density,  plot of read lengths (Dunk Illumina) - cleaned singletons")
	        polygon(d2, col="red", border="black")
	        abline(v=mymean, col="blue", lty=2, lwd=2)
	        text(x=mymean + 15, y=0.03, labels=c(paste("Mean", round(mymean,3))))
	        dev.off()

	        #length and GC distributions are done for both datasets now.
	        #clear the decks for quality scores

if(0)
{
### Quality score data

#start with average quality score (going to a different script for this, need torun on server

avg_qual <- read.table("Dunk_Illumina_QAQC.avg_qual", sep="\t")

#density plot 2D?

brewer.pal(11, "")

newScheme <- (heat.colors(14))[c(1:12, 14)]


f1 <- kde2d(dunk454[,3], dunk454[,4], h=c(bandwidth.nrd(dunk454[,3]), 2), n = 40, lims = c(0.1, 0.75, 0, 750))
pdf("joint_dist_gc_length_454_dunk.pdf", width=10, height=10)


filled.contour(f1, col = rev(newScheme), nlevels = 11, xlab="G+C content", ylab="Read length (bp)", main="Joint distribution of G+C and read length (Dunk 454)")
dev.off()

}

