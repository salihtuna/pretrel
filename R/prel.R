#' prel is used for predicting the pairwise relationships.
#'
#' prel is a function for predicting pairwise relatedness using the GENESIS package. Furthermore, this function outputs a plain file with the type of relationship for each pair that is above the selected relatedness threshold (e.g. Sample1 sample2 Parent/Offspring). Using PRIMUS (see references) software, it outputs the the list of maximum unrelated set of samples. This function can be used on its own or as part of the pretrel package.
#'
#' @param inputVCF A vcf file with the selected set of SNPs and which includes samples of interest together with the reference set of samples 
#' @param maxUnrelSet Logical, A list of maximum unrelated set will be calculated using PRIMUS. This can be turned off by setting it to FALSE. This will call the Primus(https://primus.gs.washington.edu/primusweb/) software. Make sure that Primus is installed and it is in your path
#' @param t pi_hat threshold for pairwise relationships  
#' @return This script will output the predicted ethnicities based on reference set used (PopAssigned.txt), pca plot in pdf format, eigenvalues and eigenvectors as a plain text format, and the pairwise relatedness.
#' @importFrom GENESIS pcair pcrelate pcrelateMakeGRM pcrelateReadKinship pcrelateReadInbreed
#' @importFrom GWASTools GdsGenotypeReader GenotypeData getScanID
#' @importFrom SNPRelate snpgdsVCF2GDS snpgdsOpen snpgdsClose snpgdsIBDKING
#' @importFrom gdsfmt openfn.gds closefn.gds
#' @importFrom ggplot2 ggplot theme_bw theme geom_point aes labs guides guide_legend scale_color_manual ggtitle element_text
#' @importFrom methods as
#' @importFrom utils write.table read.table
#' @author Salih Tuna, \email{st620@@cam.ac.uk} 
#' @references \url{https://bioconductor.org/packages/release/bioc/html/GENESIS.html}
#' @references \url{https://primus.gs.washington.edu/primusweb/}
#' @seealso selectSNPs(), inputPrepare(), pret(), pretrel(), assignRel()
#' @export

prel <- function(inputVCF="../mergedInput.vcf.gz",maxUnrelSet=TRUE,t=0.0125){
     
     dir.create("relatedness")
     setwd("relatedness")
     
     snpgdsVCF2GDS(vcf.fn = inputVCF, out.fn = "genotype.gds")

     genofile <- snpgdsOpen("genotype.gds")
     kinMat <- snpgdsIBDKING(genofile, type="KING-robust")
     rownames(kinMat$kinship) <- kinMat$sample.id
     colnames(kinMat$kinship) <- kinMat$sample.id
     snpgdsClose(genofile)

     #Create a genotypeObject
     geno <- GdsGenotypeReader(filename = "genotype.gds")
     genoData <- GenotypeData(geno)

     #read individual IDs from GenotypeData object
     iids <- getScanID(genoData)
     ref_IDs <- iids[grep("^NA|HG", iids)]
     wgs_IDs <- iids[-grep("^NA|HG", iids)]


     #run PC-AiR
     mypcair <- pcair(genoData = genoData,  v = 200, kinMat = kinMat$kinship, divMat = kinMat$kinship, unrel.set = ref_IDs)

     pca.df <- data.frame(mypcair$vectors)
     colnames(pca.df) <- c(paste("PC",1:200,sep=""))
     write.table(pca.df, file = "EigenVectors.txt",sep="\t",quote=FALSE)

     # run PC-Relate
     npca <- 20
     if(!file.exists("tmp_pcrelate.gds")){
	mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:npca], training.set = mypcair$unrels,write.to.gds=TRUE)
     }
     mypcrelate <- openfn.gds("tmp_pcrelate.gds")

     kinship <- pcrelateReadKinship(pcrelObj = mypcrelate, scan.include = wgs_IDs)
     kin.file <- paste("kinship_analysis_results_all_",npca,"PC.txt",sep="")
     write.table(kinship,file=kin.file,row.names=F,col.names=T,sep="\t",quote=F)

     #Multiply the results with the scaleKin
     kinship_scaled <- pcrelateMakeGRM(pcrelObj = mypcrelate, scan.include = wgs_IDs, scaleKin = 2)
     kinship.scaled.file <- paste("kinshipScaled_analysis_results_all_",npca,"PC.txt",sep="")
     write.table(kinship_scaled,file=kinship.scaled.file,row.names=F,col.names=T,sep="\t",quote=F)

     kinship_plink <- data.frame(matrix(nrow=nrow(kinship),ncol=10))
     names(kinship_plink) <- c("FID1","IID1","FID2","IID2","NA","UN","Z0","Z1","Z2","PI_HAT")
     kinship_plink$FID1 <- kinship$ID1
     kinship_plink$IID1 <- kinship$ID1
     kinship_plink$FID2 <- kinship$ID2
     kinship_plink$IID2 <- kinship$ID2
     kinship_plink$Z0 <- kinship$k0
     kinship_plink$Z1 <- kinship$k1
     kinship_plink$Z2 <- kinship$k2
     kinship_plink$PI_HAT <- kinship$kin*2
     write.table(kinship_plink,file="kinship_analysis_results_all_20PC_X2.genome",row.names=F,col.names=T,sep="\t",quote=F)       
     

     #Assign relationship
     assignRel(kinship_plink)

     inbreed <- pcrelateReadInbreed(pcrelObj = mypcrelate, scan.include = wgs_IDs)
     inbreed.file <- paste("inbreed_analysis_results_all_",npca,"PC.txt",sep="")
     write.table(inbreed,file=inbreed.file,row.names=F,col.names=T,sep="\t",quote=F)
       
     closefn.gds(mypcrelate)
     save.image("relatedness.RData")
     

if(maxUnrelSet=='TRUE'){
     system("run_PRIMUS.pl --no_PR -p kinship_analysis_results_all_20PC_X2.genome -o primus_dir_X2_Genome")
}
     setwd("../")
}
