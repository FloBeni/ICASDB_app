######################################## R generate www/name_species.tab ##########################################################################################################################

pathData <- "~/data/Projet-SplicedVariants/"
options(stringsAsFactors = F, scipen = 999) 
library(stringr)
library(readxl)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # fonction de lecture de fichier excel
  sheets <- readxl::excel_sheets(filename)
  x <-
    lapply(sheets, function(X)
      readxl::read_excel(filename, sheet = X))
  if (!tibble)
    x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
setwd(pathData)
dataSpecies = read.table("Fichiers-data/dNdS_investigation", sep="\t",header=T)

setwd('~/LBBE-Projects/Projet SplicedVariants/analyses/Interface Graphique') #pour tourner en local

mysheets <- read_excel_allsheets(paste("/home/fbenitiere/data/Projet-SplicedVariants/Fichiers-data/metazoa_species.xls",sep=""))
speciesStudied=c()
for (Species in names(mysheets)) {# recupere les especes a analyser
    speciesStudied = append(speciesStudied, Species)
}

speciesStudied=sort(speciesStudied) #ordonne les espèces par ordre alphabétique

listNomSpecies = str_replace_all(speciesStudied,"_"," ")
names(speciesStudied)=listNomSpecies

dt_species = data.frame(
  full_names=speciesStudied,
  species = names(speciesStudied),
  clades=sapply(speciesStudied,function(x) mysheets[[x]]$Clade[1]),
  nb_rnaseq = sapply(speciesStudied,function(x) length(mysheets[[x]]$`Assetion RNA-Seq`)))

mysheets <- read_excel_allsheets(paste("/home/fbenitiere/data/Projet-SplicedVariants/Fichiers-data/embryophyta_species.xls",sep=""))
speciesStudied=c()
for (Species in names(mysheets)) {# recupere les especes a analyser
    speciesStudied = append(speciesStudied, Species)
  # 
}
speciesStudied=sort(speciesStudied) #ordonne les espèces par ordre alphabétique
listNomSpecies = str_replace_all(speciesStudied,"_"," ")
names(speciesStudied)=listNomSpecies

dt_species = rbind(dt_species,data.frame(
  full_names=speciesStudied,
  species = names(speciesStudied),
  clades=sapply(speciesStudied,function(x) mysheets[[x]]$Clade[1]),
  nb_rnaseq = sapply(speciesStudied,function(x) length(mysheets[[x]]$`Assetion RNA-Seq`)))
)

table(dt_species$clades)
vector = rep("Other Invertebrates",nrow(dt_species))
vector[dt_species$clade %in% c("Amphibia","Mammalia","Aves","Lepidosauria","Testudines","Crocodylia","Teleostei","Anura","Chondrichthyes")] = "Other Vertebrates"
vector[dt_species$clade %in% c("Blattodea","Coleoptera","Ephemeroptera","Hemiptera","Thysanoptera")] = "Other Insecta"
vector[dt_species$clade %in% c("Diptera","Lepidoptera")] = "Lepido Diptera"
vector[dt_species$clade %in% c("Nematoda","Hymenoptera","Mammalia","Aves","Teleostei", "Embryophyta")] = dt_species[dt_species$clade %in% c( "Embryophyta","Hymenoptera","Nematoda","Mammalia","Aves","Teleostei"),"clades"]
dt_species$clades= vector
dt_species$clades = factor(dt_species$clades, levels = c("Embryophyta","Lepido Diptera","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Mammalia","Aves","Teleostei","Other Vertebrates"))
table(dt_species$clades)

write.table(dt_species,"name_species", sep="\t", quote=F)


######################################## Upload Image species ##########################################################################################################################


list=$(awk '{print $2}' www/name_species.tab)
for query in ${list}; do
echo "${query}"

count=1
find=true

while ${find} = true:
do
find=false
imagelink=$(wget --user-agent 'Mozilla/5.0' -qO - "www.google.be/search?q=${query}\&tbm=isch&tbs=il:cl" | sed 's/</\n</g' | grep '<img' | head -n"$count" | tail -n1 | sed 's/.*src="\([^"]*\)".*/\1/')

if wget $imagelink -O www/species_images/${query}.png ; then
    echo "Command succeeded"
    wget $imagelink -O www/species_images/${query}.png
else
    echo "Command failed"
    count=2
    find=true
fi

done
done
######################################## Upload data per species ##########################################################################################################################


list=$(awk '{print $2}' www/name_species.tab)
for query in ${list}; do
echo "${query}"
mkdir -p www/per_species_data/${query}
cp /home/fbenitiere/data/Projet-SplicedVariants/Annotations/${query}/busco_analysis/busco_to_gene_id_* www/per_species_data/${query}/
done


cp /home/fbenitiere/data/Projet-SplicedVariants/Analyses/${query}/by_intron_cds.tab www/per_species_data/${query}/by_intron_cds.tab
cp /home/fbenitiere/data/Projet-SplicedVariants/Annotations/${query}/busco_analysis/busco_to_gene_id_* www/per_species_data/${query}/
cp /home/fbenitiere/data/Projet-SplicedVariants/Analyses/${query}/by_intron_major_overlap.tab www/per_species_data/${query}/by_intron_major_overlap.tab

cp /home/fbenitiere/data/Projet-SplicedVariants/Analyses/${query}/by_gene_analysis.tab www/per_species_data/${query}/by_gene_analysis.tab



######################################## R upload phylogenetic trees ##########################################################################################################################


pathData <- "/home/fbenitiere/data/Projet-SplicedVariants/"
phylogenetic_trees = c(
  "257 species on 731 genes" = paste(pathData,"/DnDs/Metazoa_v11/RAxML/concatenatAAS_cons150.aln.raxml",sep=""),
  "227 species on 870 genes" = paste(pathData,"/DnDs/Metazoa_v9/RAxML/concatenatAAS_cons150.aln.raxml",sep=""),
  "190 species on 866 genes" = paste(pathData,"/DnDs/Metazoa_v6/RAxML/concatenatAAS_cons100.aln.raxml",sep=""),
  "175 species on 194 genes" = paste(pathData,"/DnDs/Eukaryota_v7/RAxML/concatenatAAS.aln.raxml.rooted",sep=""),
  "140 species on 236 genes" = paste(pathData,"/DnDs/Eukaryota_v6/RAxML/concatenatAAS.aln.raxml.rooted",sep=""),
  "128 species on 862 genes" = paste(pathData,"/DnDs/Metazoa_v5/RAxML/concatenatAAS.aln.raxml.rooted",sep=""),
  "118 species on 854 genes" = paste(pathData,"/DnDs/Metazoa_v4/RAxML/concatenatAAS.aln.raxml.rooted",sep=""),
  "67 species on 896 genes" = paste(pathData,"/DnDs/Metazoa_v2/RAxML/concatenatAAS.aln.raxml.rooted",sep=""),
  "61 species on 909 genes" = paste(pathData,"/DnDs/Metazoa/RAxML/concatenatAAS.aln.raxml.rooted",sep=""),
  "53 species on 917 genes" = paste(pathData,"/DnDs/Metazoa_53sp_Busco_alignment_53sp_prank/DnDs_replicate_aas_tree_917genes/concatenatAAS.aln.raxml.rooted",sep=""),
  "53 species on 16 genes" = paste(pathData,"/DnDs/Metazoa_53sp_genes_in_every_species/prank_concatenat/concatenatAAS.phy_phyml_tree_rooted.txt.rooted",sep=""),
  "69 species on 6 genes" = paste(pathData,"/DnDs/Metazoa_69sp/pgls_tree.rooted",sep=""))

dt = data.frame()
for (name in names(phylogenetic_trees)){
  i = phylogenetic_trees[name]
  split = str_split(i,"/")[[1]]
  file.copy(i, paste("/home/fbenitiere/ICASDB_app/UI_v2/www/phylogenetic_tree/",split[length(split)-2],"_", split[length(split)-1],"_",split[length(split)],sep=""), overwrite = TRUE)

  dt = rbind(dt,data.frame(
    name = paste(split[length(split)-2],"_", split[length(split)-1],"_",split[length(split)],sep=""),
    description = name
  ))
}
write.table(dt,"/home/fbenitiere/ICASDB_app/UI_v2/www/phylogenetic_tree/tree_description.tab",col.names = T,row.names = F,quote=F,sep="\t")

