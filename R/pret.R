#' PRedicting EThnicity and RELatedness (pretrel) from genotyping data.
#'
#' pret is used for predicting the ethnicity of the samples of interest. It uses Principal Component Analysis (PCA). GENESIS package is used as the base for appliying pca to a vcf file. Once the projections are obtained by pc-air function from GENESIS package, samples are assigned to one of the predefined clusters, from the reference set. 
#' @param inputVCF Merged vcf file that has samples from both reference and the project set with the selected set of SNPs. This file can be obtained by running the inputPrepare() function of the pretrel package. 
#' @return This script will output PopAssigned.txt, for each sample showing the predicted ethnicity, pca plot in pdf format, eigenvalues and eigenvectors as a plain text format.
#' @importFrom GENESIS pcair
#' @importFrom GWASTools GdsGenotypeReader GenotypeData getScanID
#' @importFrom SNPRelate snpgdsVCF2GDS snpgdsClose
#' @importFrom ggplot2 ggplot theme_bw theme geom_point aes labs guides guide_legend scale_color_manual ggtitle element_text ggsave
#' @importFrom utils write.table read.table data
#' @importFrom stats median mahalanobis cov
#' @author Salih Tuna, \email{st620@@cam.ac.uk} 
#' @references \url{https://bioconductor.org/packages/release/bioc/html/GENESIS.html}
#' @seealso selectSNPs(), inputPrepare(), prel(), pretrel()
#' @export

pret <- function(inputVCF){
	PC1 <- PC2 <- megapopMerged <- NULL
	if (!dir.exists("ethnicity")){
		dir.create("ethnicity")
	}	
	setwd("ethnicity")	

	snpgdsVCF2GDS(vcf.fn = inputVCF, out.fn = "genotype.gds")
	
	# Create a genotypeObject
	geno <- GdsGenotypeReader(filename = "genotype.gds")
	genoData <- GenotypeData(geno)

	# read individual IDs from GenotypeData object
	iids <- getScanID(genoData)

	# run PC-AiR
	ref_IDs <- iids[grep("^NA|HG", iids)]
	wgs_IDs <- iids[-grep("^NA|HG", iids)]


	#getGenotype(genoData)[1:5,1:5]
	mypcair <- pcair(genoData = genoData, v = 200, unrel.set = ref_IDs)

	write.table(mypcair$values, file = "EigenValues_Ethnicity.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names="EigenValues")

	pca.df <- data.frame(mypcair$vectors)
	colnames(pca.df) <- c(paste("PC",1:200,sep=""))
	write.table(pca.df, file = "EigenVectors_Ethnicity.txt",sep="\t",quote=FALSE)
	# hap.pops=read.table("integrated_call_samples_v3.20130502.ALL.panel",header=1)
	hap.pops=read.table("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",header=1)

	megapop=c()
	megapop[hap.pops$super_pop == c("EUR")] = "European"
	megapop[hap.pops$pop == c("FIN")] = "Finnish-European"
	megapop[hap.pops$super_pop == c("EAS")] = "East-Asian"
	megapop[hap.pops$super_pop == c("AFR")] = "African"
	megapop[hap.pops$super_pop == c("SAS")] = "South-Asian"

	pca.df$pop=hap.pops$pop[match(rownames(pca.df),hap.pops$sample)]
	pca.df$megapop=as.factor(megapop[match(rownames(pca.df),hap.pops$sample)])

	# 1000G
	pca.df.ref <- pca.df[grep("^NA|HG",rownames(pca.df)),]
	# test
	pca.df.wgs <- pca.df[-grep("^NA|HG",rownames(pca.df)),]

	npcaEtn <- 5

	my.pops.main <- c("European","Finnish-European","East-Asian","African","South-Asian")
	centers.main <- mat.or.vec(length(my.pops.main),npcaEtn)
	for (i in 1 : length(my.pops.main)){
       	    for (j in 1 : npcaEtn){
		centers.main[i,j] <- median(subset(pca.df.ref[,j],pca.df.ref$megapop==my.pops.main[i]))
	    }
	}

	colnames(centers.main) <- colnames(centers.main, do.NULL = FALSE, prefix = "PC")
	rownames(centers.main) <- my.pops.main


	my.pops.main <- c("European","Finnish-European","East-Asian","African","South-Asian")
	cov.list <- list()
	for(i in 1 : length(my.pops.main)){
	      cov.list[[my.pops.main[i]]] <- cov(subset(pca.df.ref[,1:5], pca.df.ref$megapop==my.pops.main[i]))
	}


	likelihoods <- matrix(NA,dim(pca.df.wgs)[1],5)
	pop <- c("European","Finnish-European","East-Asian","African","South-Asian")
	for (i in 1:5){
	    mah.dis<- mahalanobis(pca.df.wgs[,1:5], centers.main[i,],cov.list[[i]] )
	    likelihoods[,i] <- exp(-0.5 * mah.dis)/sqrt(det(cov.list[[i]]))
	}

	pop.assigned <- matrix()

	for (s in 1:dim(pca.df.wgs)[1]){
	    if(sum(as.numeric(likelihoods[s,] > 10000000)) >= 1){
	    	pop.assigned[s] <- which.max(likelihoods[s,] )
    	    }
	}

	pop.assigned[is.na(pop.assigned)] <- 0

	tmp <- cbind(pca.df.wgs,pop.assigned)
	if(any(unique(tmp$pop.assigned)==1)){
		tmp[which(tmp$pop.assigned == 1),]$megapop <- "European"
	}
	if(any(unique(tmp$pop.assigned)==2)){
		levels(tmp$megapop) <- c(levels(tmp$megapop), "Finnish-European")
		tmp[which(tmp$pop.assigned == 2),]$megapop <- "Finnish-European"
	}
	if(any(unique(tmp$pop.assigned)==3)){
		tmp[which(tmp$pop.assigned == 3),]$megapop <- "East-Asian"
	}
	if(any(unique(tmp$pop.assigned)==4)){
		tmp[which(tmp$pop.assigned == 4),]$megapop <- "African"
	}
	if(any(unique(tmp$pop.assigned)==5)){
		tmp[which(tmp$pop.assigned == 5),]$megapop <- "South-Asian"
	}
	if(any(unique(tmp$pop.assigned)==0)){
		levels(tmp$megapop) <- c(levels(tmp$megapop), "Other")
		tmp[which(tmp$pop.assigned == 0),]$megapop <- "Other"
	}

	write.table(likelihoods,"Pop_likelihood_score.txt",row.names=rownames(tmp),col.names= c("European","Finnish-European","East-Asian","African","South-Asian"),quote=F,sep="\t")

	pca.df.ref$megapopMerged <- pca.df.ref$megapop
	pca.df.ref$megapopMerged[which(pca.df.ref$megapopMerged == "Finnish-European")] <- "European"

	plot.pca <- ggplot() + ggplot2::theme_bw() + ggplot2::theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5)) + ggplot2::geom_point(data=pca.df.ref, aes(x=PC1, y=PC2,colour=megapopMerged),size=8, alpha=0.05) + ggplot2::labs(x="Principal Component 1", y="Principal Component 2") + ggplot2::guides(colour = guide_legend(title="Population",override.aes = list(size=3, alpha = 1))) + ggplot2::scale_color_manual(values=c("African"="#F8766D","East-Asian"="#7CAE00", "European"="#00BFC4", "South-Asian"="#C77CFF"),labels=c("African","East-Asian","European","South-Asian")) + ggplot2::geom_point(data=pca.df.wgs, aes(x=PC1, y=PC2),shape=21) + ggplot2::ggtitle("All samples with 1000G")
	
	ggsave("pca_ethnicity.pdf", plot = plot.pca, device = "pdf", path = "./", dpi=300)

	sample.labels <- tmp[,"megapop",drop=FALSE]
	sample.labels <- cbind(Sample = rownames(sample.labels), sample.labels)
	rownames(sample.labels) <- NULL
	write.table(sample.labels,file="PopAssigned.txt",sep="\t",quote=F,row.names=FALSE)

	save.image("ethnicity.RData")

	#snpgdsClose(genofile)
	#snpgdsClose()
	
	setwd("../")
}
