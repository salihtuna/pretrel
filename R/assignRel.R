#' Assigns the type of relationships for each pair that is above the selected relatedness threshold.
#'
#' Part of the pretrel package. This functions asssigns the type of relationship for each pair above the selected relatedness threshold. The input file needs to be in the plink pairwise relatedness format.
#'
#' @param kinshipFile plink format pairwise relatedness table. This can be obtrained from the prel function. 
#' @param t pi_hat threshold for pairwise relationships 
#' @return This script will output a plain txt file for the related pairs and the type of relationships the pair has.
#' @importFrom utils write.table globalVariables
#' @author Salih Tuna, \email{st620@@cam.ac.uk} 
#' @seealso selectSNPs(), inputPrepare(), pret(), prel(), pretrel()
#' @export


assignRel <- function(kinshipFile, t=0.0125){
		kin.ratios <- c(0,0.25,0.50,0.75,1)
		tmp <- kinshipFile[kinshipFile['PI_HAT'] >= t,]
		for (i in 1 : nrow(tmp)){
    
			z0 <- kin.ratios[which.min(abs(tmp[i,'Z0']-kin.ratios))]
			z1 <- kin.ratios[which.min(abs(tmp[i,'Z1']-kin.ratios))]
			z2 <- kin.ratios[which.min(abs(tmp[i,'Z2']-kin.ratios))]
			pihat <- kin.ratios[which.min(abs(tmp[i,'PI_HAT']-kin.ratios))]

    			# Assign relationship
    			# Idenical
    			if ((z0==0) && (z1==0) && (z2==1) && (pihat==1)){
       			  tmp[i,"relationship"] = "Duplicate or MZ twin"
    			# Parent/offspring
     			} else if ((z0==0) && (z1==1) && (z2==0) && (pihat==0.5)){
       			  tmp[i,"relationship"] = "Parent/Offspring" 
    			# Sibling
     			} else if ((z0==0.25) && (z1==0.5) && (z2==0.25) && (pihat==0.5)){
       			  tmp[i,"relationship"] = "Full Sibling" 
    			# Grandparent or half sibling or aunt or nephew
     			} else if ((z0==0.5) && (z1==0.5) && (z2==0) && (pihat==0.25)){
       			  tmp[i,"relationship"] = "Half sibling or grandparent or aunt or nephew" 
    			# Cousin
     			} else if ((z0==0.75) && (z1==0.25) && (z2==0)){
       			  tmp[i,"relationship"] = "First Cousin" 
    			# Unknown
      			} else{
      			  tmp[i,"relationship"] = "Unknown" 
      			}
		}
		write.table(tmp,file="relationship_assigned.txt",col.names=TRUE, row.names=FALSE, sep="\t",quote=F)
}
