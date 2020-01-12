#' Selects SNPs that will be used for the ethnicity and relatedness analysis
#'
#' selectSNPs is a series of bcftools commands. So it relies on bcftools and GNU Parallel. These two softwares need to be installed and be in your path. Please see below where to get these softwares. SelectSNPs function can be used on its own or as part of the pretrel() function. This approach requires that both the reference set of samples and project set of samples are mapped to the same reference version (i.e. either both is v37 or both is v38). 
#'
#' @param inputRAW_VCF unfiltered vcf file(s) to be used for the analysis.
#' @param RefDataDir Path for the reference set.
#' @param tmpDir Directory to write the temporary files.
#' @param nnode number of nodes to run the jobs in parallel.
#' 
#' @return This script will output a list of SNPs in plain text format. This file will be used as input for the inputPrepare() function.
#'
#' @author Salih Tuna, \email{st620@@cam.ac.uk} 
#' @references \url{https://samtools.github.io/bcftools/bcftools.html}
#' @references \url{https://www.gnu.org/software/parallel/}
#' @export

selectSNPs <- function(inputRAW_VCF, RefDataDir, tmpDir, nnode=2){
	   SNPsFiltered <- NULL
	   # Check bcftools can be run
           if(system("bcftools -v")!=0) {
                stop("bcftools not in your PATH")
           }
	   
	   dir.create("selectSNPs")
	   setwd("selectSNPs")
	   
	   system(paste("bcftools view -GH ",inputRAW_VCF, " | cut -f1-2 | sort -Vk1 | uniq -u > uniqPos",sep=""))
	   system(paste("bcftools view -R uniqPos -T uniqPos ", inputRAW_VCF, " -f .,PASS -O z -o uniq_PASS.vcf.gz",sep=""))
	   system(paste("bcftools sort -T ", tmpDir, " uniq_PASS.vcf.gz -O z -o uniq_PASS_sorted.vcf.gz",sep=""))
	   system("bcftools index uniq_PASS_sorted.vcf.gz")
	   system("bcftools view -v snps uniq_PASS_sorted.vcf.gz -O z -o uniq_PASS_sorted_snps.vcf.gz")
	   system("bcftools index uniq_PASS_sorted_snps.vcf.gz")
	   system("bcftools +missing2ref uniq_PASS_sorted_snps.vcf.gz -O z -o uniq_PASS_sorted_snps_noMiss.vcf.gz")
	   system("bcftools +fill-AN-AC uniq_PASS_sorted_snps_noMiss.vcf.gz -O z -o uniq_PASS_sorted_snps_noMiss_updAN.vcf.gz")

	   #Remove SNPs that have any missing genotypes and Remove SNPs whose MAF is <0.30 in full NIHR BR-RD dataset (minus QC excluded samples)
	   sampleNo <- as.numeric(system("bcftools query -l uniq_PASS_sorted_snps_noMiss.vcf.gz | wc -l", intern=TRUE))
	   system(paste("bcftools filter uniq_PASS_sorted_snps_noMiss_updAN.vcf.gz -i 'AN=",sampleNo*2," && MAF >= 0.30' -O z -o uniq_PASS_sorted_snps_noMiss_updAN_maf.vcf.gz",sep=""))

	   # "Indexing ..."
	   system("bcftools index uniq_PASS_sorted_snps_noMiss_updAN_maf.vcf.gz")
	   system("bcftools view -GH uniq_PASS_sorted_snps_noMiss_updAN_maf.vcf.gz | cut -f1-2 > maf.txt")

	   #Remove SNPs that are in multiallelic sites in 1000G set

	   refFiles.list <-  Sys.glob(file.path(RefDataDir, "*.vcf.gz"))
	   refFile <- gsub("chr\\d+", "chr{}", refFiles.list)[1]
	   system(paste("for c in {1..22}; do echo $c; done | parallel --verbose -j",nnode," bcftools view -R maf.txt -T maf.txt ", refFile, " -O z -o chr{}_refSamples_snplist_TMP.vcf.gz",sep=""))
    	   system(paste("for c in {1..22}; do echo $c; done | parallel --verbose -j",nnode," bcftools norm -m-both chr{}_refSamples_snplist_TMP.vcf.gz -O z -o All.chr{}_refSamples_snplist_norm_TMP.vcf.gz",sep=""))
	   system(paste("for c in {1..22}; do echo $c; done | parallel --verbose -j",nnode," bcftools view -GH All.chr{}_refSamples_snplist_norm_TMP.vcf.gz | cut -f1-2 | sort | uniq -u >> SNPsFiltered_tmp",sep=""))
	   

	   system("sort -Vk1 SNPsFiltered_tmp > SNPsFiltered")
	   system("bcftools view -R SNPsFiltered -T SNPsFiltered uniq_PASS_sorted_snps_noMiss_updAN_maf.vcf.gz -O z -o uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered.vcf.gz")
	   system("bcftools index uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered.vcf.gz")

	   #Check whether ref and alt are mathching in the ref and test set
	   
	   #Filter reference set
	   refFiles.list <-  Sys.glob(file.path(RefDataDir, "*.vcf.gz"))

	   if (length(refFiles.list) >= 22 || length(refFiles.list <= 25)){
	      refFile <- gsub("chr\\d+", "chr{}", refFiles.list)[1]
   	      system(paste("for c in {1..22}; do echo $c; done | parallel --verbose -j",nnode, " bcftools view -R SNPsFiltered -T SNPsFiltered", refFile, "-O z -o chr{}_refSamples_snplist_TMP.vcf.gz",sep=" "))
	      refFiles.filtered <- Sys.glob(file.path("./", "*_refSamples_snplist_TMP.vcf.gz"))
   	      system(paste("bcftools concat", refFiles.filtered[1], refFiles.filtered[2], refFiles.filtered[3], refFiles.filtered[4], refFiles.filtered[5], refFiles.filtered[6], refFiles.filtered[7], refFiles.filtered[8], refFiles.filtered[9], refFiles.filtered[10], refFiles.filtered[11], refFiles.filtered[12], refFiles.filtered[13], refFiles.filtered[14], refFiles.filtered[15], refFiles.filtered[16], refFiles.filtered[17], refFiles.filtered[18], refFiles.filtered[19], refFiles.filtered[20], refFiles.filtered[21], refFiles.filtered[22], "-O z -o RefSamplesMerged_snplist_TMP.vcf.gz",sep=" "))
	   }else if (length(refFiles.list) == 1){
   	      system(paste("bcftools view -R", SNPsFiltered, "-T", SNPsFiltered, refFiles.list, "-O z -o RefSamplesMerged_snplist_TMP.vcf.gz", sep=" "))
	   }
	   system("bcftools index RefSamplesMerged_snplist_TMP.vcf.gz")
			 
	   system("bcftools isec -p dir -n=2 -w1 uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered.vcf.gz RefSamplesMerged_snplist_TMP.vcf.gz")
	   system("bcftools view dir/0000.vcf -O z -o dir/wgs10k_filtered.vcf.gz")
	   system("bcftools index dir/wgs10k_filtered.vcf.gz")
	   system("bcftools view -GH dir/wgs10k_filtered.vcf.gz | cut -f1-2,4-5 > snpListTEMP.txt")


	   #Perform LD pruning using the PLINK command
	   system("bcftools view uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered.vcf.gz -Ou | bcftools annotate --set-id '%CHROM:%POS' -O z -o uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered_annot.vcf.gz")

	   system("plink --vcf uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered_annot.vcf.gz --double-id --chr 1-22 --indep-pairwise 50 5 0.2 --out uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered")

	   system("sed 's/:/\t/g' uniq_PASS_sorted_snps_noMiss_updAN_maf_filtered.prune.in > snpList.txt")

	   
	   file.copy('snpList.txt', '../')
	   setwd("../")
	   unlink('selectSNPs',recursive=TRUE)
}
