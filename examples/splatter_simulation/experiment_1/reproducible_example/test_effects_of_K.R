#### Function for testing the effects of k.param and dimension parameter on the KNN-based similarity matrix ####
test_effects_of_K <- function(k.param, dimensions,seuratobject) {
  processedseurat <- processSeuratObject(seuratobject = seuratobject,dimensions = dimensions,resolution=0.1,K=k.param)
  KNNgraph <- processedseurat[[2]]
  celltypelabels <- Idents(processedseurat[[3]])
  KNN_based_graph <- KNN_transition(KNNgraph,celltypelabels)
  return(KNN_based_graph)
}