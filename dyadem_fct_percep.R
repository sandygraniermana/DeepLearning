
###################################################### PIPELINE POUR DATA CYTOMETRIE EN FLUX #####################################################



################# FONCTIONS ###################



#NA

na_data <- function(flow.frame){
  na = flow.frame[,-which(is.na(data@parameters@data$desc))]
  return(na)
}



#matrix

matrix_data <- function(flow.frame){
  mm = model.matrix(~as.factor(data@exprs[, dim(flow.frame)[2]]-1))
  return(mm)
}



#compensation des data

compensate_data <- function(flow.frame){
  comp = compensate(flow.frame, data@description[["SPILL"]])
  return(comp)
}



#transformation des data

logiclTransformCIPHE <- function(flow.frame, value = NULL, marker = NULL){
  
  if(is.null(marker)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- marker
  }
  
  list.index <- names(unlist(lapply(markers.transform, function(x)
    return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(is.null(value) || is.na(value)){
    if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])){
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
      )
    } else
      if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
        r.values <- unlist(lapply(list.index, function(x)
          as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
        )
      } else {
        r.values <- rep(90, length(list.index))
      }
  }
  else {
    r.values <- rep(value, length(list.index))
  }
  
  w.values <- (4.5-log10(262143/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5
  
  for(t in 1:length(markers.transform)){
    lgcl <- logicleTransform(w=w.values[t])
    flow.frame <- transform(flow.frame,
                            transformList(markers.transform[t],lgcl))
  }
  
  return(flow.frame)
}




#traitement des data compensation + transformation

preprocess_data = function(dossier){
  fcss = dir(path = dossier, pattern = 'fcs$', full.names = TRUE)
  efcs = NULL
  ann = lefcs = list()[1:length(fcss)]
  #ann1=lefcs=list()[1:length(fcss)]
  j = 0
  for (fcs_name in fcss){
    j = j+1
    print(fcs_name)
    fcs = flowCore::read.FCS(fcs_name, emptyValue = FALSE)
    ann[[j]]=exprs(fcs)[,dim(fcs)[2]] #derniere colonne des fichiers (annotations)
    #ann1[[j]]=exprs(fcs)[,dim(fcs)[2]] #derniere colonne des fichiers (annotations)
    fcs = compensate_data(fcs)
    fcs = logiclTransformCIPHE(fcs)
    fcs = na_data(fcs)
    #fcs=matrix_data(fcs)
    
    
    efcs = rbind(efcs,exprs(fcs)) 
    lefcs[[j]] = exprs(fcs)
  }
  return(list(efcs = efcs, lefcs = lefcs, ann = ann))
  #return(list(efcs=efcs,lefcs=lefcs,ann1=ann1)) #lefcs = data pour les 6 fichiers dont chacun 12 colonnes et ann1 = pop des 6 fichiers et efcs = data mais 3 listes(?)
}




#UMAP

Umap = function(flow.frame){
  u = umap::umap(flow.frame, verbose = TRUE)
  return(u)
}





#pop 20000 cellules pour la nouvelle annotation

taille_pop1 <- function(edata, Pop){ #car connait edata
  
  res = NULL
  for (i in sort(unique(Pop))){ #ds chaque pop
    inds = sample(x = 1:length(which(Pop == i)), size = 20000, replace = TRUE) #x = longueur de pop.   replace = T -> des valeurs peuvent etre reprises
    Inds = which(Pop == i) #position
    taille = Inds[inds]
    res = c(res, taille) #vecteur qui recupere resultat et taille à chq tour de boucle
    
  }
  res <- sample(x = res, size = length(res), replace = FALSE) #pour que la derniere population ne soit pas enlevee avec la validation_split=0.1, qui prend 10% de la derniere population (colonne)
  
  return(list(data = edata[res,], annot = Pop[res])) #fichier avec res + pop avec res
  
}



#classification avec les autres fichiers

percept = function(dossier, n, seuil = 0.0){  
  
  edata = preprocess_data(dossier = dossier) #fct de traitement des data
  
  um = annot = list()[1:length(edata$lefcs)]
  
  for (i in 1:length(edata$lefcs)){ #prediction : avec autres fichiers
    
    print(dim(edata$lefcs[[i]]))
    # u=Umap(edata$lefcs[[i]])
    # um[[i]]=u$layout
    
    print(i) # 6 fichiers
    
    # t=taille_pop(edata$lefcs[[i]])
    
    #res_loc=n[[2]](k$centers)  -> kmeans
    #res_loc=n[[2]](f$codes)   #-> flowsom
    #res_loc=n[[2]](edata$lefcs[[i]]) # -> sans clusters
    res_loc = n[[2]](edata$lefcs[[i]]) # -> sans clusters + taille
    
    center_annot = apply(res_loc$numpy(), MARGIN = 1, FUN = function(x, seuil){
      ii = which.max(x)
      ifelse(x[ii] > seuil, ii, length(x)+1)
    }, seuil = seuil) #voir arrangement des annotations : le max
    

    #cell_annot=center_annot # -> sans clusters
    annot[[i]] = center_annot  
    
  }
  
  return(list(annot, edata, um, res_loc$numpy()))
}




#visualisation des annotations

vis <- function(ann, annot){ #annot = nouvelle annotation
  #population des fichiers (du 1er au 6e)
  res = list()[1:length(annot)]
  for (i in 1:length(annot)){ #les 6 fichiers
    col1 = diag(table(factor(ann[[1]][[i]], levels = 1:35), annot[[i]]))/colSums(table(factor(ann[[1]][[i]], levels = 1:35), annot[[i]]))
    col2 = colSums(table(factor(ann[[1]][[i]], levels = 1:35), annot[[i]]))
    col3 = rowSums(table(factor(ann[[1]][[i]], levels = 1:35), annot[[i]]))
    col4 = table(factor(ann[[1]][[i]],levels = 1:35), annot[[i]])[,35]
    res[[i]] = data.frame(precision = col1, taille = col2, taille_annot = col3, ungated = col4, prop_ungated = col4/col3)
  }
  
  return(res)
}








# #comparaison vraie pop et faux ungated
# 
# comp <- function(ann, annot){ #FCS 4
#   
#   data3 = read.FCS("/home/workspace/data_annote/IMPC1_V5_A5_S SSC_004_celltype.fcs",emptyValue = FALSE)
#   edata3=exprs(logiclTransformCIPHE(compensate(data3,data3@description[["SPILL"]]))) #compensation + transformation des donnees
#   edata3=edata3[,-which(is.na(data@parameters@data$desc))] #donnees sans les NA
#   
#   #umap_edata3=umap::umap(edata3,verbose=TRUE) #umap
#   res_umap_popV= res_umap_ungatedF = list()[4] #resultats fcs 4
#   p = plot.default(umap_edata3[[1]], pch='.', col=annot[[4]]) #umap des nouvelles annotations du fcs 4
#   
#   for (i in 1:length(unique(annot[[4]]))){
#     print(i)
#     png(paste("pop_", i,".png", sep = ""), width = 4*480, height = 4*480)
#     par(mfrow = c(1, 2))
#     p1 = plot(umap_edata3[[1]],col = rep("grey",nrow(umap_edata3[[1]])), pch = ".")
#     P1 = points(umap_edata3[[1]][annot[[4]] == levels(annot[[4]])[i],],col = "gold",pch = 19,cex = 0.5)
#     p2 = plot(umap_edata3[[1]],col = rep("grey",nrow(umap_edata3[[1]])), pch=".")
#     P2 = points(umap_edata3[[1]][annot[[4]] == "ungated" & ann[[1]][[4]] == i,], col = "blue", pch = 19,cex = 0.5)
#     dev.off()
#     
#   }
#   
#   return(0)
# }
# 
# 
# c = comp(ann, annot)
# 
# 
# # 
# # #umap_edata3=umap::umap(edata3,verbose=TRUE) #faire umap du fcs 4
# # 
# # ann = percept(dossier = "/home/workspace/data_annote/", n = train[[1]], seuil = 0.5)   #return(list(annot, edata, um, res_loc$numpy()))
# # 
# # plot.default(umap_edata3[[1]], pch='.', col=annot[[4]]) #umap des nouvelles annotations du fcs 4
# # 
# # plot(umap_edata3[[1]],col=rep("grey",nrow(umap_edata3[[1]])), pch=".") #umap des datas du fcs 4 (gris)
# # ann[[1]][[4]]==35 #annotation predite des ungated du fcs 4
# # ind=which(ann[[1]][[4]]==35)    
# # points(umap_edata3[[1]][ind,], col=annot[[4]],pch=19,cex=0.5) #umap des mauvais ungated predits avec couleurs annotations correctes
# # 
# # 
# # 
# # 
# # points(umap_edata3[[1]][annot[[4]]=="ungated",],col="gold",pch=".")
# # 
# # 
# # par(mfrow=c(1,2))
# # plot.default(umap_edata3[[1]], pch='.', col=annot[[4]]) #umap des nouvelles annotations du fcs 4
# # plot(umap_edata3[[1]],col=rep("grey",nrow(umap_edata3[[1]])), pch=".") #"verification : tout est gris"
# # points(umap_edata3[[1]][annot[[4]]=="ungated",],col=ann[[1]][[4]][annot[[4]]=="ungated"],pch=".") #umap des ungated predits annotes
# # 
# # # #a faire pour les 35 pop
# # # col=rep("grey",nrow(umap_edata3[[1]]))
# # # plot(umap_edata3[[1]],col=col, pch=".") #umap verification : tout est gris
# # # col[annot[[4]]==levels(annot[[4]])[2]]="gold" #pop coloree
# # 
# # #vraie pop
# # plot(umap_edata3[[1]],col=rep("grey",nrow(umap_edata3[[1]])), pch=".")
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[1],],col="gold",pch=19,cex=0.5)
# # #ungated annotés comme pop (ungated faux)
# # plot(umap_edata3[[1]],col=rep("grey",nrow(umap_edata3[[1]])), pch=".")
# # points(umap_edata3[[1]][annot[[4]]=="ungated" & ann[[1]][[4]]==1,], col="blue",pch=19,cex=0.5)
# # 
# # 
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[2],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[3],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[4],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[10],],col="gold",pch=19)
# # 
# # 
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[10],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[10],] ,col=ann[[1]][[4]][annot[[4]]=="ungated"],pch=".")
# # 
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[13],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[18],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[22],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[26],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[29],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[30],],col="gold",pch=19)
# # points(umap_edata3[[1]][annot[[4]]==levels(annot[[4]])[35],],col="gold",pch=19)
# # 
# 
# 
# 
# comp <- function(ann, annot){ #FCS 1
#   
#   data = read.FCS("/home/workspace/data_annote/IMPC1_V5_A2_S SSC_001_celltype.fcs",emptyValue = FALSE)
#   edata=exprs(logiclTransformCIPHE(compensate(data,data@description[["SPILL"]]))) #compensation + transformation des donnees
#   edata=edata[,-which(is.na(data@parameters@data$desc))] #donnees sans les NA
#   
#   umap_edata=umap::umap(edata,verbose=TRUE) #umap
#   # res_umap_popV= res_umap_ungatedF = list()[4] #resultats fcs 4
#   p = plot.default(umap_edata[[1]], pch='.', col=annot[[1]]) #umap des nouvelles annotations du fcs 4
#   
#   for (i in 1:length(unique(annot[[1]]))){
#     print(i)
#     png(paste("pop_",i,".png",sep=""),width=4*480,height=4*480) #resultat en png
#     par(mfrow=c(1,2)) #2 images sur 1 ligne et 2 colonnes
#     p1 = plot(umap_edata[[1]],col=rep("grey",nrow(umap_edata[[1]])), pch=".") #verification, tout est gris
#     P1 = points(umap_edata[[1]][annot[[1]]==levels(annot[[1]])[i],],col="gold",pch=19,cex=0.5) #annotations vraies
#     p2 = plot(umap_edata[[1]],col=rep("grey",nrow(umap_edata[[1]])), pch=".") #verification, tout est gris
#     P2 = points(umap_edata[[1]][annot[[1]]=="ungated" & ann[[1]][[1]]==i,], col="blue",pch=19,cex=0.5) #faux ungated : ungated annotes comme une pop
#     dev.off()
#     
#   }
#   
#   return(0)
# }
# 
# 
# c = comp(ann, annot)




