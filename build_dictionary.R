################################################################################
###    R script to preprocess Illumina Human Methylation EPIC annotations    ###
###                            Clovis Chabert                                ###
###                           Florent Chuffart                               ###
################################################################################

################################################################################
###                         Loading annotations                              ###
################################################################################

library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
pf = data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
pf$cg_id = rownames(pf)


################################################################################
###                Indexing probes using platform infos                      ###
################################################################################

if (!exists("gene_indexed_probes_df")) {
  gene_indexed_probes_df = apply(pf, 1, function(l) {
    # l = pf["cg05912063",]
    # l = pf["cg21301391",]
    # l = pf["cg01068023",]
    # l = pf["cg09383789",]
    # l = pf["cg15633699",]
    # l = pf["cg25294651", ]
    l
    # pf[1:100,c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
    if (is.na(l[["UCSC_RefGene_Name"]])) {
      return(NULL)
    } else if (l[["UCSC_RefGene_Name"]]=="") {
      return(NULL)      
    } else if (length(grep(";",l[["UCSC_RefGene_Name"]]))==0) {
      res = c(l[["UCSC_RefGene_Name"]], l[["UCSC_RefGene_Group"]], l[["cg_id"]], l[["chr"]])
    } else {
      res = cbind(
        strsplit(l[["UCSC_RefGene_Name"]], ";")[[1]],  
        strsplit(l[["UCSC_RefGene_Group"]], ";")[[1]]
      )
      res = cbind(res, l[["cg_id"]], l[["chr"]])
      res = res[!duplicated(paste0(res[,1], "_", res[,2])),]
    }
    
    res
  })
  
  length(gene_indexed_probes_df)
  gene_indexed_probes_df = gene_indexed_probes_df[!sapply(gene_indexed_probes_df, is.null)]
  length(gene_indexed_probes_df)
  
  gene_indexed_probes_df = do.call(rbind, gene_indexed_probes_df)
  rownames(gene_indexed_probes_df) = NULL
  
  dim(gene_indexed_probes_df)
  head(gene_indexed_probes_df)
  # gene_indexed_probes_df[,4] = paste0("chr", gene_indexed_probes_df[,4])
}

length(unique(gene_indexed_probes_df[,1]))
length(unique(gene_indexed_probes_df[,3]))
table(gene_indexed_probes_df[,4])
table(gene_indexed_probes_df[,2])
positions = unique(gene_indexed_probes_df[,2])
positions

# WARNING! check if there is not a gene on two chr

if (!exists("gene_indexed_probes_list")) {
  chrs = unique(gene_indexed_probes_df[,4])
  gene_indexed_probes_list = apply(t(t(chrs)), 1, function(chr) {
    #chr = "chr1"
    print(chr)
    sub_gene_indexed_probes_df = gene_indexed_probes_df[gene_indexed_probes_df[,4] %in% chr,]
    #print(dim(sub_gene_indexed_probes_df))
    genes = unique(sub_gene_indexed_probes_df[,1])
    # print(length(genes))
    foo = apply(t(t(genes)), 1, function(gene) {
      # gene = "TCEANC2"
      # gene = "LOC100133331"
      # print(gene)
      idx = sub_gene_indexed_probes_df[,1] %in% gene
      sub_sub = sub_gene_indexed_probes_df[idx,]
      if (sum(idx) == 1) {
        sub_sub = t(sub_sub)
      }
      bar = lapply(positions, function(pos) {
        sub_sub[sub_sub[,2] %in% pos,3]
      })
      names(bar) = positions
      bar$PROMOTER = sub_sub[sub_sub[,2] %in% c("TSS200", "TSS1500", "5'UTR"),3]
      length(bar)
      bar
    })
    #foo
    names(foo) = genes
    foo
  })
  
  names(gene_indexed_probes_list) = chrs
}


################################################################################
###                                 Results                                  ###
################################################################################


# Probes are indexed by chromosome names...
print(head(names(gene_indexed_probes_list)))
# then by gene names... 
print(head(names(gene_indexed_probes_list[[1]])))
# and finally by relative postion according to gene.
print(gene_indexed_probes_list[[1]][[1]])










