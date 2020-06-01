#perceptron



training <- function(annot, edata, use_ungated = FALSE){
  
  
  if (!use_ungated) {   #on garde que les gated
    ss = which(annot != "ungated")   #position des booleens avec ungated = FALSE
    print(length(ss))
    print(dim(ss))
    
    edata = edata[ss,]   #data des gated
    print(dim(edata))
    
    Pop = annot[ss]   #pop des gated
    
  } else {  #avec les ungated
    
    Pop = annot   #pop gated + pop ungated
    
  }
  
  
  tt = taille_pop1(edata, as.factor(as.character(Pop))) #pop avec chacune 20000 cellules
  
  mm = model.matrix(~as.factor(tt$annot)-1) #dummy variable 
  
  reticulate::source_python("/home/workspace/SG/neuralnet.py") #lien vers script python 
  
  n = neuralnet(tt$data, mm) #training
  
  
  return(list(n))
  
}




prediction <- function(train, annot, dossier = "/home/workspace/data_annote/"){
  
  print(str(train))
  
  
  ann = percept(dossier = dossier, n = train[[1]], seuil = 0.5)   #prediction annotations
  
  visualisation = vis(ann, annot)   #visualisation des predictions
  
  
  return(list(ann, visualisation))
  
}

