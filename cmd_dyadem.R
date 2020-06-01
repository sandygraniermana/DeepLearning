library(flowCore)


data = read.FCS("/home/workspace/data_annote/IMPC1_V5_A2_S SSC_001_celltype.fcs", emptyValue = FALSE) #FCS 1 pour training

edata = exprs(logiclTransformCIPHE(compensate(data,data@description[["SPILL"]]))) #compensation + transformation des data
edata = edata[,-which(is.na(data@parameters@data$desc))] #donnees sans les NA


reticulate::source_python("/home/workspace/SG/neuralnet.py")

train = training(annot[[1]], edata, use_ungated = TRUE)  #1er fichier FCS nouvellement annote #fct training
#train = training(annot[[1]], edata, use_ungated = FALSE)

pred = prediction(train, annot)  #prediction des annotations

