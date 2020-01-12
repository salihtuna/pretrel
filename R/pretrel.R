#' PRedicting EThnicity and RELatedness (pretrel) from genotyping data.
#'
#' pretrel function is a wrapper for four functions that starts with selectSNP() function for selecting SNPs that will be used throughout the analysis. If the user would like to user their own list of SNPs, then this selectSNPs() function can be ignored. Next function within pretrel is to run inputPrepare() which will prepare the input file using the reference set of samples, snp list and the raw vcf file. Next step is to predict the ethnicities(pret()) and predict the pairwise relatedness and also selects the set of maximum unrelated samples (prel()). All these functions can be used seperately.
#'
#' @param snpList List of SNPs that will be used to filter the inputRAW and reference set
#' @param inputRAW_VCF unfiltered vcf file(s) to be used for the analysis
#' @param RefDataDir Path to the 1000G reference set (ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz)
#' @param inputVCF A vcf file with the selected set of SNPs and which includes samples of interest together with the reference set of samples 
#' @return This script will output the predicted ethnicities based on reference set used (PopAssigned.txt), pca plot in pdf format, eigenvalues and eigenvectors as a plain text format, and the pairwise relatedness.
#' @importFrom GENESIS pcair
#' @importFrom GWASTools GdsGenotypeReader GenotypeData getScanID
#' @importFrom SNPRelate snpgdsVCF2GDS
#' @importFrom ggplot2 ggplot theme_bw theme geom_point aes labs guides guide_legend scale_color_manual ggtitle element_text
#' @importFrom utils globalVariables
#' @references \url{https://bioconductor.org/packages/release/bioc/html/GENESIS.html}
#' @references \url{https://primus.gs.washington.edu/primusweb/}
#' @seealso selectSNPs(), inputPrepare(), pret(), prel()
#' @export

pretrel <- function(snpList, inputRAW_VCF, RefDataDir, inputVCF){
	# SelectSNPs
	selectSNPs(inputRAW_VCF, RefDataDir)
	
	# Prepera the Input File
	inputPrepare(snpList, inputRAW_VCF, RefDataDir)

	# For ethnicity
	pret(inputVCF)

	# For relatedness
	prel(inputVCF)
}
