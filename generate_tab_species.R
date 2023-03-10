options( stringsAsFactors = F, scipen = 999 )

path = "/home/fbenitiere/data/"

metazoa_species_clade_lht = rbind( read.table(paste(path , "Projet-SplicedVariants/Fichiers-data/metazoa_species_clade_lht.tab",sep=""),header=T) ,
                                   read.table(paste(path , "Projet-SplicedVariants/Fichiers-data/metazoa_species2_clade_lht.tab",sep=""),header=T),
                                   read.table(paste(path , "Projet-SplicedVariants/Fichiers-data/embryophyta_species_clade_lht.tab",sep=""),header=T) )
rownames(metazoa_species_clade_lht) = metazoa_species_clade_lht$species

table(metazoa_species_clade_lht$clade)
metazoa_species_clade_lht[,"clade_group"] = "Other Invertebrates"
metazoa_species_clade_lht[metazoa_species_clade_lht$clade %in% c("Amphibia","Mammalia","Aves","Lepidosauria","Testudines","Crocodylia","Teleostei","Anura","Chondrichthyes"),"clade_group"] = "Other Vertebrates"
metazoa_species_clade_lht[metazoa_species_clade_lht$clade %in% c("Blattodea","Coleoptera","Ephemeroptera","Hemiptera","Thysanoptera"),"clade_group"] = "Other Insecta"
metazoa_species_clade_lht[metazoa_species_clade_lht$clade %in% c("Diptera","Lepidoptera"),"clade_group"] = "Lepido Diptera"
metazoa_species_clade_lht[metazoa_species_clade_lht$clade %in% c("Nematoda","Hymenoptera","Mammalia","Aves","Teleostei","Embryophyta"),"clade_group"] = metazoa_species_clade_lht[metazoa_species_clade_lht$clade %in% c("Hymenoptera","Nematoda","Mammalia","Aves","Teleostei","Embryophyta"),"clade"]
metazoa_species_clade_lht$clade_group = factor(metazoa_species_clade_lht$clade_group, levels = c("Lepido Diptera","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Mammalia","Aves","Teleostei","Other Vertebrates","Embryophyta"))
table(metazoa_species_clade_lht$clade_group)

# restrain = metazoa_species_clade_lht[!sapply(metazoa_species_clade_lht$species,function(x) file.exists(paste(path, "Projet-SplicedVariants/Analyses/",x,"/summary_characteristics.tab",sep="")) ),]
# 
# restrain$got_overlap = sapply(restrain$species,function(x) file.exists(paste(path, "Projet-SplicedVariants/Analyses/",x,"/by_intron_major_overlap.tab",sep="")) )
# restrain$got_gene = sapply(restrain$species,function(x) file.exists(paste(path, "Projet-SplicedVariants/Analyses/",x,"/by_gene_analysis.tab",sep="")) )
# restrain$got_gene = sapply(restrain$species,function(x) file.exists(paste(path, "Projet-SplicedVariants/Analyses/",x,"/by_gene_analysis.tab",sep="")) )
# restrain = restrain[restrain$got_gene,]
# restrain$got_gc = sapply(restrain$species,function(x) file.exists(paste(path, "Projet-SplicedVariants/Annotations/",x,"/GC_content.tab",sep="")) )
# restrain$got_intron = sapply(restrain$species,function(x) file.exists(paste(path, "Projet-SplicedVariants/Analyses/",x,"/by_intron_cds.tab",sep="")) )
# 
# paste(restrain[restrain$got_overlap,]$species,collapse="','")


get_CM_dNdS_true <- function(D) {
  cum_KS=sum(D$num_dS)
  cum_KN=sum(D$num_dN)
  cum_OS=sum(D$den_dS/D$branch_length)
  cum_ON=sum(D$den_dN/D$branch_length)
  
  cum_dS=cum_KS/cum_OS
  cum_dN=cum_KN/cum_ON
  cum_dNdS=cum_dN/cum_dS
  return(c(cum_dN,cum_dS,cum_dNdS))
}


# species = 'Cervus_elaphus'
# RECUPERE COLONNE
all_dt = c()
for (species in metazoa_species_clade_lht$species){
  if (file.exists(paste(path, "Projet-SplicedVariants/Analyses/",species,"/summary_characteristics.tab",sep=""))){
    dt = read.delim(paste(path, "Projet-SplicedVariants/Analyses/",species,"/summary_characteristics.tab",sep=""))
    if( "longevity;quant" %in% dt$label){print(species)
      # file.remove(paste(path, "Projet-SplicedVariants/Analyses/",species,"/summary_characteristics.tab",sep=""))
      }
    all_dt = c(dt$label,all_dt)
  }
}


columns = unique(all_dt)



all_dt = data.frame()
for (species in metazoa_species_clade_lht$species){print(species)
  if (file.exists(paste(path, "Projet-SplicedVariants/Analyses/",species,"/summary_characteristics.tab",sep=""))){
    dt = read.delim(paste(path, "Projet-SplicedVariants/Analyses/",species,"/summary_characteristics.tab",sep=""))
    rownames(dt) = dt$label
    if (metazoa_species_clade_lht[species,]$clade == "Embryophyta"){
      dt = dt[!grepl("_metazoa",(rownames(dt))),]
    } else {
      dt = dt[!grepl("_embryophyt",(rownames(dt))),]
    }
    
    all_dt = rbind(all_dt,dt[columns,]$value)
  }
}
colnames(all_dt) = columns


method_list = c("Embryophyta_v1"="Embryophyta_v1/subset_200_ksites_GC3/data_calculation.tab",
                "Eukaryota_v6"="Eukaryota_v6/subset_200_ksites_GC3/data_calculation.tab",
                "Eukaryota_v7"="Eukaryota_v7/subset_200_ksites_GC3/data_calculation.tab",
                "Metazoa_v11"="Metazoa_v11/subset_200_ksites_GC3/data_calculation.tab",
                "Metazoa_v9"="Metazoa_v9/subset_200_ksites_GC3/data_calculation.tab",
                "Metazoa_v6"="Metazoa_v6/subset_200_ksites_GC3/data_calculation.tab",
                "Metazoa_v5"="Metazoa_v5/subset_200_ksites_GC3/data_calculation.tab",
                "Metazoa_v4"="Metazoa_v4/subset_200_ksites/data_calculation.tab",
                "Metazoa_species_filtered_v53.2"="Metazoa_species_filtered_v53.2/subset_200_ksites_GC3/data_calculation.tab")

for (method in names(method_list)){ print(method)
  data_dNdS = read.delim(paste(path,"Projet-SplicedVariants/DnDs/",method_list[method],sep=""))
  data_dNdS = data_dNdS[!is.na(data_dNdS$branch_length),]
  
  for (species in unique(all_dt$species)){
    data_sequence = data_dNdS[ data_dNdS$species == species,]
    dnds_values = get_CM_dNdS_true(data_sequence)
    all_dt[all_dt$species == species,paste(method,"_dn;quant",sep="")] = dnds_values[1]
    all_dt[all_dt$species == species,paste(method,"_ds;quant",sep="")] = dnds_values[2]
    all_dt[all_dt$species == species,paste(method,"_dnds;quant",sep="")] = dnds_values[3]
  }
}

all_dt$`longevity;quant` = metazoa_species_clade_lht[all_dt$species,]$lifespan
all_dt$`length;quant` = metazoa_species_clade_lht[all_dt$species,]$length
all_dt$`weight;quant` = metazoa_species_clade_lht[all_dt$species,]$weight


write.table(all_dt,paste("/home/fbenitiere/Desktop/ICASDB_app/www/species_informations_tables/data_by_species.tab",sep=""),sep="\t",quote=F,row.names = F, col.names = T)

