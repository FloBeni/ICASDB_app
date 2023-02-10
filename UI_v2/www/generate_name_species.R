


options(stringsAsFactors = F, scipen = 999) 
library(readxl)
library(stringr)


read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <-
    lapply(sheets, function(X)
      readxl::read_excel(filename, sheet = X))
  if (!tibble)
    x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
list_excel_documents = c("/home/fbenitiere/data/Projet-SplicedVariants/Fichiers-data/metazoa_species.xls",
                         "/home/fbenitiere/data/Projet-SplicedVariants/Fichiers-data/embryophyta_species.xls")

dt_species = data.frame()
for (excel_path in list_excel_documents){
  mysheets <- read_excel_allsheets(excel_path)
  for (species in names(mysheets)){
    
    dt_species = rbind(dt_species,data.frame(
      full_names = species,
      species = str_replace_all(species,"_",""),
      clades = mysheets[[species]]$Clade[1],
      nb_rnaseq = length(mysheets[[species]]$`SRA Accession ID`)
    ))
    
  }
}

table(dt_species$clades)
vector = rep("Other Invertebrates",nrow(dt_species))
vector[dt_species$clade %in% c("Amphibia","Mammalia","Aves","Lepidosauria","Testudines","Crocodylia","Teleostei","Anura","Chondrichthyes")] = "Other Vertebrates"
vector[dt_species$clade %in% c("Blattodea","Coleoptera","Ephemeroptera","Hemiptera","Thysanoptera")] = "Other Insecta"
vector[dt_species$clade %in% c("Diptera","Lepidoptera")] = "Lepido Diptera"
vector[dt_species$clade %in% c("Nematoda","Hymenoptera","Mammalia","Aves","Teleostei", "Embryophyta")] = dt_species[dt_species$clade %in% c( "Embryophyta","Hymenoptera","Nematoda","Mammalia","Aves","Teleostei"),"clades"]
dt_species$clades= vector
dt_species$clades = factor(dt_species$clades, levels = c("Embryophyta","Lepido Diptera","Hymenoptera","Other Insecta","Nematoda","Other Invertebrates","Mammalia","Aves","Teleostei","Other Vertebrates"))
table(dt_species$clades)

rownames(dt_species)=dt_species$full_names
dt_species = dt_species[sort(dt_species$full_names),]

write.table(dt_species,"name_species.tab", sep="\t", quote=F)


