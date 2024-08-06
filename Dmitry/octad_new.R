#BiocManager::install("octad")

library('octad')

phenoDF=get_ExperimentHub_data('EH7274')

breast_ph=subset(phenoDF, biopsy.site=='BREAST')
breast_primary = breast_ph[grepl('BRCA', breast_ph$loss_list),]
case_id=breast_primary$sample.id

breast_adjacent=subset(phenoDF,cancer=='breast invasive carcinoma'&sample.type == 'adjacent')
control_id=breast_adjacent$sample.id

res=diffExp(case_id,control_id) 

sRGES = runsRGES(res,max_gene_size=100,permutations=1000)

sRGESf = sRGES[sRGES$sRGES < -0.2,]
sRGESf = sRGESf[sRGESf$n > 3,]
sRGESf = sRGESf[!grepl("BRD-", sRGESf$pert_iname),]

current_directory <- getwd() 

folderEnrich = 'octad_demo'
dir.create(folderEnrich)
octadDrugEnrichment(sRGES=sRGES, target = c('chembl_targets','mesh','ChemCluster'), outputFolder = folderEnrich)
list.files(paste0(current_directory, "/octad_demo/enrichFolder")

