library(topGO)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"

source(paste0(project_dir, "scripts/general/R/GO_mappings_to_list.R"))

get_gsea_results <- function(stats, gene2GO, mode = "ks", alpha = 0.05, 
                             nodeSize = 20, ontology = "BP"){
  
  if(mode == "fisher"){
    topGOobject <- new("topGOdata", ontology = ontology,
                       allGenes = stats,
                       nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
    
    resultsFisher <- runTest(topGOobject, algorithm = "classic", statistic = "fisher")
    #resultsFisher.weight01 <- runTest(topGOobject, algorithm = "weight01", statistic = "fisher")
    
    topNodes = length(score(resultsFisher))
    
    allRes <- GenTable(topGOobject, 
                       classicFisher = resultsFisher, orderBy = "classicFisher", 
                       ranksOf = "classicFisher", topNodes = topNodes, numChar = 500)
    
    allRes$classicFisher <- as.numeric(str_remove(allRes$classicFisher, "< "))
    #allRes$weight01FisherPadj = p.adjust(allRes$weigh01Fisher, method = "BH")
    allRes$classicFisherPadj <- p.adjust(allRes$classicFisher, method = "BH")
    
    
  } else if(mode == "ks"){
    topGOobject <- new("topGOdata", ontology = ontology,
                       allGenes = stats, geneSel = function(x) x < alpha,
                       nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
    
    resultsKS <- runTest(topGOobject, algorithm = "classic", statistic = "ks")
    resultsKS.weight01 <- runTest(topGOobject, algorithm = "weight01", statistic = "ks")
    
    topNodes = length(score(resultsKS.weight01))
    
    allRes <- GenTable(topGOobject, weight01KS = resultsKS.weight01, 
                       classicKS = resultsKS, orderBy = "weight01KS", 
                       ranksOf = "weight01KS", topNodes = topNodes, numChar = 500)
    
    allRes$weight01KSpadj = p.adjust(allRes$weight01KS, method = "BH")
    allRes$classicKSpadj = p.adjust(allRes$classicKS, method = "BH")
  }
  
  return(allRes)
}
