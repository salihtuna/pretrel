#' Prepares the input VCF file for ethnicity and relatedness analysis by using selected set of SNPs.
#'
#' InputPrepare function prepares the input VCF file for ethnicity and relatedness analysis by extracting the set of SNPs from the reference set and the set of samples and merge the two together. It is a series of bcftools commands and some of which run in parallel. So it relies on bcftools and GNU Parallel softwares. These two softwares need to be installed and be in your path. Please see 'References' section below where to get these softwares. The function makes sure that exact same of SNPs exist in both the test samples and the reference set. Currently this runs only on 1000Genome reference set(see references for details).
#'
#' @param snpList List of SNPs that will be used to filter the inputRAW and reference set. By default this is the output file from selectSNPs() function. You can enter a list of snps (ex 1	12345)
#' @param inputRAW_VCF unfiltered vcf file(s) to be used for the analysis
#' @param RefDataDir Path for the reference set
#' @param nnode number of nodes to be used for the parallel job execution. Default is two nodes.
#' 
#' @return This script will output a merged VCF file that contains both project and reference samples. The merged VCF file consists of only the SNPs selected either by the selectSNPs function or the user input.
#' 
#' @author Salih Tuna, \email{st620@@cam.ac.uk} 
#' @references \url{https://samtools.github.io/bcftools/bcftools.html}
#' @references \url{https://www.gnu.org/software/parallel/}
#' @references \url{https://www.internationalgenome.org/data/}
#'
#' @seealso selectSNPs(), pret(), prel(), pretrel()
#' @export


inputPrepare <- function(snpList, inputRAW_VCF, RefDataDir, nnode=2){
	excludeSamples_1000G <- NULL
	dir.create("inputPrepare")
	setwd("inputPrepare")
	data(excludeSamples_1000G,envir = environment())
	write.table(excludeSamples_1000G,"excludeSamples_1000G.txt", row.names=FALSE,col.names=FALSE,quote=F)    
	
	# Check bcftools can be run
	if(system("bcftools -v")!=0) {
    		stop("bcftools not in your PATH")
	}

		
	#Filter the samples of interest
	system(paste("bcftools view -R", snpList, "-T", snpList, inputRAW_VCF, "-O z -o snpList_genoTMP.vcf.gz", sep=" "))
	system("bcftools index snpList_genoTMP.vcf.gz")

	#Filter reference set
	refFiles.list <-  Sys.glob(file.path(RefDataDir, "*.vcf.gz"))

	if (length(refFiles.list) >= 22 || length(refFiles.list <= 25)){
	   refFile <- gsub("chr\\d+", "chr{}", refFiles.list)[1]
   	   system(paste("for c in {1..22}; do echo $c; done | parallel --verbose -j",nnode, " bcftools view -R", snpList, "-T", snpList, refFile, "-O z -o chr{}_refSamples_snplist_TMP.vcf.gz",sep=" "))
   	   refFiles.filtered <- Sys.glob(file.path("./", "*_refSamples_snplist_TMP.vcf.gz"))
   	   system(paste("bcftools concat", refFiles.filtered[1], refFiles.filtered[2], refFiles.filtered[3], refFiles.filtered[4], refFiles.filtered[5], refFiles.filtered[6], refFiles.filtered[7], refFiles.filtered[8], refFiles.filtered[9], refFiles.filtered[10], refFiles.filtered[11], refFiles.filtered[12], refFiles.filtered[13], refFiles.filtered[14], refFiles.filtered[15], refFiles.filtered[16], refFiles.filtered[17], refFiles.filtered[18], refFiles.filtered[19], refFiles.filtered[20], refFiles.filtered[21], refFiles.filtered[22], "-O z -o RefSamplesMerged_snplist_TMP.vcf.gz",sep=" "))
	}else if (length(refFiles.list) == 1){
   	   system(paste("bcftools view -R", snpList, "-T", snpList, refFiles.list, "-O z -o RefSamplesMerged_snplist_TMP.vcf.gz", sep=" "))
	}
	system("bcftools index RefSamplesMerged_snplist_TMP.vcf.gz")

	system("bcftools isec -p dir -n=2 -w1 snpList_genoTMP.vcf.gz RefSamplesMerged_snplist_TMP.vcf.gz")
	system("bcftools view dir/0000.vcf -O z -o dir/wgs10k_filtered.vcf.gz")
	system("bcftools index dir/wgs10k_filtered.vcf.gz")
	if (grepl('1000G', refFile)){
	   system("bcftools merge -m none RefSamplesMerged_snplist_TMP.vcf.gz dir/wgs10k_filtered.vcf.gz -O z -o Ref1000G-all_WGS_TMP.vcf.gz")
	   system("bcftools view -S ^excludeSamples_1000G.txt Ref1000G-all_WGS_TMP.vcf.gz -O z -o mergedInput.vcf.gz") 
	}else{
	   system("bcftools merge -m none RefSamplesMerged_snplist_TMP.vcf.gz dir/wgs10k_filtered.vcf.gz -O z -o mergedInput.vcf.gz")
	}
	
       	file.copy('mergedInput.vcf.gz', '../')
	setwd("../")
	unlink('inputPrepare', recursive=TRUE)
}
