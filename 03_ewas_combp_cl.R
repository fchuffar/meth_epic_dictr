####Exploitation donnees methylation CECILE 13/01/2022: extension auc CpGs autour de ceux identifies par Xu

#reprendre a partir d'ici pr les tests
library(MASS)
library(lmtest)
library(sandwich)
library(survey)
library(dplyr)

if(!exists("tabC_VF")){
  load(file="CL/tabC_VF.RData")
  names(tabC_VF)
  dim(tabC_VF) 
}


#si on veut faire sur les 381 CpG de Xu
#cpg_xu_cecile <- read.csv(file="//10.91.165.35/partage_commun$/USERS/CANCER_ENVIRONNEMENT/C_LEMARCHAND/Methylation/Version_finale/cpg_xu_cecile.csv")
#cpg_xu_cecile <- read.csv(file="//10.91.165.35/partage_commun$/USERS/CANCER_ENVIRONNEMENT/C_LEMARCHAND/Methylation/Version_finale/"nom_fichier".csv")
#head(cpg_xu_cecile)
#length(cpg_xu_cecile)
#dim(cpg_xu_cecile)
#class(cpg_xu_cecile)
#class(probe)


probes=readRDS("extended_probes_xu.rds")
length(probes)
probes=intersect(probes,colnames(tabC_VF))
length(probes)

#cpg_list_mojgan <- list()
#probes=cpg_xu_cecile$cpg
colnames(tabC_VF)[1:40]
sub_tabC_VF=tabC_VF[,c(colnames(tabC_VF)[1:29],probes)]
dim(sub_tabC_VF)
#for (i in 1:length(cpg_xu_cecile[[1]])) { #nb de lignes ds la colonne 1
#  cpg_list_mojgan[i] <- colnames(sub_tabC_VF[cpg_xu_cecile[[1]][i]]) #pr recuperer les noms des cpg
#  print(paste0("number of loop :", i))
#} 







#creation dun vecteur expo
expos=c(
  "no2_j0_6",
  "no2_j0_59",
  "conc_no2_GA",
  "moy_PM10_7j",
  "moy_PM10_60j",
  "moy_PM10_365j",
  "moy_PM25_7j",
  "moy_PM25_60j",
  "moy_PM25_365j"
)
expos



for (expo in expos) {
  print(paste0("************",expo))

  # ewas
  ewas = epimedtools::monitored_apply(t(t(probes)), 1, function(probe) { # rlm sur cpg selectionnes
    i=which (probe==probes)
    tmpm <- MASS::rlm(sub_tabC_VF[,probe] ~ sub_tabC_VF[,expo] + ageref + as.factor(dep) + as.factor(Fumeur) + as.factor(batch) + as.factor(row) + as.factor(col), data=sub_tabC_VF, maxit = 100)
    reg_beta <- coeftest(tmpm)[2, 1] 
    reg_p <- coeftest(tmpm)[2, 4] 
    return(c(beta=reg_beta, pval=reg_p))
  })
  ewas = t(ewas)
  rownames(ewas) = probes
  head(ewas)  
  dim(ewas)
  
  # export ewas results
  pval = ewas[,2]  
  bed = pf[names(pval),1:2] # Warning, here pf is global, it must be arg and indexed!!
  head(bed)
  bed$end = bed[,2] + 1
  bed$probes = names(pval)
  bed$pval = pval
  bed$strand = "+"
  colnames(bed) = c("chrom", "start", "end", "probe", "pval", "strand")
  head(bed)
  dim(bed)
  
  
  bed_ewas_filename = paste0("ewas4combp_", expo, ".bed")
  bed[,1] = as.character(bed[,1])
  bed = bed[order(bed[,1], bed[,2]),]
  write.table(bed,file=bed_ewas_filename , sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


for (expo in expos) {
  # Run comb-p, run!!
  pv_thresh = 0.01
  bed_ewas_filename = paste0("ewas4combp_", expo,".bed")
  cmd = "comb-p"
  arg = paste0("pipeline -c 5 --seed ", pv_thresh, " --dist 1000 -p dmr_", expo,"_", pv_thresh, " --region-filter-p 0.05 --region-filter-n 2 ", bed_ewas_filename)
  print(paste(cmd, arg))
  system2(cmd, arg)
}
