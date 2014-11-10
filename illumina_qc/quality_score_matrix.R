

if(1)
{
	quality_scores_raw <- read.table("../Dunk_Illumina_mate_pairs_updated_Illumina_qualities.avg_qual.subset", sep="\t")

	tile <- 1
	x <- 2
	y <- 3
	readInPair <- 4
	quality_score_average <- 5

	pdf("hist_Illumina_dunk_raw_mean_quality_across_reads.pdf", width=10, height=10)

	hist(quality_scores_raw[,quality_score_average], breaks=250, xlab="Mean base quality across read", ylab="Frequency", main="Histogram of mean base-call quality across read (Illumina - Dunk)")

	dev.off()

	quality_scores_trimmed_pairs <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_mate_pairs.avg_qual.subset", sep="\t")
	pdf("hist_Illumina_dunk_trimmed_pairs_mean_quality_across_reads.pdf", width=10, height=10)
	hist(quality_scores_trimmed_pairs[,quality_score_average], breaks=250, xlab="Mean base quality across read", ylab="Frequency", main="Histogram of mean base-call quality across trimmed read (Illumina - Dunk)")
        dev.off()

	quality_scores_singletons <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_solitary_reads.avg_qual.subset", sep="\t")
	pdf("hist_Illumina_dunk_trimmed_singletons_mean_quality_across_reads.pdf", width=10, height=10)
	hist(quality_scores_singletons[,quality_score_average], breaks=250, xlab="Mean base quality across read", ylab="Frequency", main="Histogram of mean base-call quality across trimmed singleton read (Illumina - Dunk)")
        dev.off()
}



myColours <- heat.colors(11)
getColor <- function(x, myColours, breaks)
{
	for(i in 1: length(breaks))
	{
		if(breaks[i] >= x)
		{
			return(myColours[i])
		}
	}
}

#using correct qualities
quality_scores <- read.table("../Dunk_Illumina_mate_pairs_updated_Illumina_qualities.avg_qual.subset", sep="\t")
tile <- 1
x <- 2
y <- 3
readInPair <- 4
quality_score_average <- 5

pdf("Dunk_illumina_raw_avg_quality_score_distribution.pdf")
hist(quality_scores[,quality_score_average], breaks=250, main="Avg. quality score distribution (Dunk Illumina - trimmed mate-pairs)", xlab="Quality score")
dev.off()

minX <- min(quality_scores[,x]) -1
minY <- min(quality_scores[,y]) -1 
maxX <- max(quality_scores[,x]) + 2500
maxY <- max(quality_scores[,y]) + 1

#easiest way is to plot four points
pdf("spatial_distribution_qual_scores_across_flow_dunk.pdf", width=10, height=10)
mypointsX <- c(minX, maxX, minX, maxX)
myPointsY <- c(minY, minY, maxY, maxY)
plot(x=mypointsX, y=myPointsY, col="white", xlab="X-position", ylab="Y-position",
main="Spatial distribution of avg. quality score across the flow cell (Dunk)")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
"black") 

breaks <- c(5,10,15,20,21,23,25,27,30,35,40)

start <- 1
numRows <- dim(quality_scores)[1]

#dim(quality_scores)[1]

for(i in start:numRows)
{
	myColor <- getColor(quality_scores[i,quality_score_average], myColours, breaks);
	print(myColor)
	print(quality_scores[i,quality_score_average])
	points(x=quality_scores[i,x], y=quality_scores[i,y], col=myColor, pch=16, cex=0.5);
}

legend("bottomright", legend=rev(c("5", "10","15", "20","21","23","25","27","30","35","40")), col=rev(myColours), pch=16, cex=1, bty="n", text.col="white")
dev.off()

#plot trimmed reads, avg quality
quality_scores <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_mate_pairs.avg_qual.subset", sep="\t")

tile <- 1
x <- 2
y <- 3
readInPair <- 4
quality_score_average <- 5

pdf("Dunk_illumina_avg_quality_score_trimmed_distribution.pdf")
hist(quality_scores[,quality_score_average], breaks=250, main="Avg. quality score distribution (Dunk Illumina - quality trimmed)", xlab="Quality score")
dev.off()


minX <- min(quality_scores[,x]) -1
minY <- min(quality_scores[,y]) -1
maxX <- max(quality_scores[,x]) + 2500
maxY <- max(quality_scores[,y]) + 1


#easiest way is to plot four points
pdf("spatial_distribution_qual_scores_across_flow_trimmed_dunk.pdf", width=10, height=10)

mypointsX <- c(minX, maxX, minX, maxX)
myPointsY <- c(minY, minY, maxY, maxY)
plot(x=mypointsX, y=myPointsY, col="white", xlab="X-position", ylab="Y-position",
main="Spatial distribution of avg. quality score across the flow cell (Dunk - trimmed)")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="black")

#how many colours do I need?

breaks <- c(5,10,15,20,21,23,25,27,30,35,40)

start <- 1
numRows <- dim(quality_scores)[1]
#dim(quality_scores)[1]

for(i in start:numRows)
{
        myColor <- getColor(quality_scores[i,quality_score_average], myColours, breaks);
        print(myColor)
        print(quality_scores[i,quality_score_average])
        points(x=quality_scores[i,x], y=quality_scores[i,y], col=myColor, pch=16, cex=0.5);
}

legend("bottomright", legend=rev(c("5", "10","15", "20","21","23","25","27","30","35","40")), col=rev(myColours), pch=16, cex=1, bty="n", text.col="white")
dev.off()


#Plot singletons reads
quality_scores <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_solitary_reads.avg_qual.subset", sep="\t")
tile <- 1
x <- 2
y <- 3
readInPair <- 4
quality_score_average <- 5

pdf("Dunk_illumina_avg_quality_score_solitary_trimmed_distribution.pdf")
hist(quality_scores[,quality_score_average], breaks=250, main="Avg. quality score distribution (Dunk Illumina - quality trimmed singletons)", xlab="Quality score")
dev.off()

#Illumina_mate_pairs_updated_Illumina_qualities_trimmed_mate_pairs.avg_qual.subset
minX <- min(quality_scores[,x]) -1
minY <- min(quality_scores[,y]) -1
maxX <- max(quality_scores[,x]) + 2500
maxY <- max(quality_scores[,y]) + 1


#easiest way is to plot four points
pdf("spatial_distribution_singleton_qual_scores_across_flow_trimmed_dunk.pdf", width=10, height=10)
intsX <- c(minX, maxX, minX, maxX)
myPointsY <- c(minY, minY, maxY, maxY)
plot(x=mypointsX, y=myPointsY, col="white", xlab="X-position", ylab="Y-position",
main="Spatial distribution of avg. quality score across the flow cell (Dunk)")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
"black")


#how many colours do I need?

breaks <- c(5,10,15,20,21,23,25,27,30,35,40)

#whats this subsettin about??
start <- 1
numRows <- dim(quality_scores)[1]
#dim(quality_scores)[1]

for(i in start:numRows)
{
        myColor <- getColor(quality_scores[i,quality_score_average], myColours, breaks);
        print(myColor)
        print(quality_scores[i,quality_score_average])
        points(x=quality_scores[i,x], y=quality_scores[i,y], col=myColor, pch=16, cex=0.5);
}

legend("bottomright", legend=rev(c("5", "10","15", "20","21","23","25","27","30","35","40")), col=rev(myColours), pch=16, cex=1, bty="n", text.col="white")
dev.off()

#Pairs(forward - reverse) plots

#untrimmed
qpairs <- read.table("..//Dunk_Illumina_mate_pairs_updated_Illumina_qualities.paired_qualities.subset", sep="\t");
#wilcox.exact(qpairs[,5], qpairs[,6], paired=TRUE)

#says the means are different

#how can we show how the means are different on a plot??
tile <- 1
x <- 2
y <- 3
readInPair <- 4
quality_score_first <- 5
quality_score_second <- 6

minX <- min(qpairs[,x]) -1
minY <- min(qpairs[,y]) -1 
maxX <- max(qpairs[,x]) + 2500
maxY <- max(qpairs[,y]) + 1

mypointsX <- c(minX, maxX, minX, maxX)
myPointsY <- c(minY, minY, maxY, maxY)

pdf("spatial_distribution_of_paired_read_quality_differences.pdf", width=10, height=10)
plot(x=mypointsX, y=myPointsY, col="white", xlab="X-position", ylab="Y-position",
main="Spatial distribution: differences in paired-read quality scores (Dunk)")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
"#666666") 

getColor <- function(x, myColours, breaks)
{
	for(i in 1: length(breaks))
	{
		if(breaks[i] >= x)
		{
			return(myColours[i])
		}
	}
}



breaks <- c(-40, -30, -25, -20, -15, -10, -5, 0,5, 10, 15, 20, 25, 30, 40)

#colours not enough
colour_vec <- c("#FF0000", "#DB0000", "#B60000", "#920000", "#6D0000", "#490000", "#240000","#000000","#002400","#004900","#006D00","#009200","#00B600",
"#00DB00","#00FF00")

#colour_vec <- redgreen(13)

#colour_vec <- brewer.pal(11, "RdYlGn")
#colour_vec <- c(colour_vec[1:6], "#FFFFFF", colour_vec[7:11])

start <- 1
numRows <- dim(qpairs)[1]

for(i in start:numRows)
{
	myColor <- getColor( (qpairs[i,quality_score_first] - qpairs[i,quality_score_second]), colour_vec, breaks);
	#print(myColor)
	#print(qpairs[i,quality_score_average])
	points(x=qpairs[i,x], y=qpairs[i,y], col=myColor, pch=16, cex=0.5);
}

legend("bottomright", legend=rev(c("-40","-30","-25","-20","-15","-10", "-5", "0","5","10","15","20","25","30","40")), col=rev(colour_vec), pch=16, cex=1, bty="n", text.col="white")
dev.off()



# mates trimmed

qpairs <- read.table("Dunk_Illumina_updated_qualities_no_length_cutoff_trimmed_mate_pairs.avg_qual.paired_qualities.subset", sep="\t");
#wilcox.exact(qpairs[,5], qpairs[,6], paired=TRUE)

#says the means are different

#how can we show how the means are different on a plot??
tile <- 1
x <- 2
y <- 3
readInPair <- 4
quality_score_first <- 5
quality_score_second <- 6

minX <- min(qpairs[,x]) -1
minY <- min(qpairs[,y]) -1
maxX <- max(qpairs[,x]) + 2500
maxY <- max(qpairs[,y]) + 1

mypointsX <- c(minX, maxX, minX, maxX)
myPointsY <- c(minY, minY, maxY, maxY)

pdf("spatial_distribution_of_trimmed_paired_read_quality_differences.pdf", width=10, height=10)
plot(x=mypointsX, y=myPointsY, col="white", xlab="X-position", ylab="Y-position",
main="Spatial distribution: differences in trimmed paired-read quality scores (Dunk)")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
"#666666")

getColor <- function(x, myColours, breaks)
{
        for(i in 1: length(breaks))
        {
                if(breaks[i] >= x)
                {
                        return(myColours[i])
                }
        }
}
start <- 1
numRows <- dim(qpairs)[1]

for(i in start:numRows)
{
        myColor <- getColor( (qpairs[i,quality_score_first] - qpairs[i,quality_score_second]), colour_vec, breaks);
        #print(myColor)
        #print(qpairs[i,quality_score_average])
        points(x=qpairs[i,x], y=qpairs[i,y], col=myColor, pch=16, cex=0.5);
}

legend("bottomright", legend=rev(c("-40", "-30","-25","-20","-15","-10", "-5", "0","5","10","15","20","25","30","40")), col=rev(colour_vec), pch=16, cex=1, bty="n", text.col="white")
dev.off()

