
######################################### fonctions pour annotation perceptron ###################################


#les commandes sont dans le fichier : "/home/data/SG/cmd_annot_percep.R"


library("flowCore")
library("FlowCIPHE")


### traitement des data ###


#NA

na_data <- function(flow.frame){
  na = flow.frame[,-which(is.na(flow.frame@parameters@data$desc))]
  return(na)
}



#compensation des data

compensate_data <- function(flow.frame){
  comp = compensate(flow.frame,flow.frame@description[["SPILL"]])
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



#bonnes annotations

fct<-function(dir_fcs){
  
  flow.frames <- list.files(dir_fcs, pattern =".fcs", full.names = TRUE, recursive = FALSE)#liste des fichiers
  
  pop<-as.character(gsub(".fcs$","",gsub("^.*//","",flow.frames)))
  
  flow.frames <- lapply(flow.frames, function(i) {read.FCS(i)}) #infos des fichiers
  
  #names(flow.framesB) <- basename(list.files("/home/data/SG/Concat_gated_clean/Camille/", pattern = ".fcs", recursive = FALSE)) #noms des fichiers dans le chemin
  fcs <- concatenate.FCS.CIPHE(flow.frames, params = "Annotation") #infos + concatenation des annotations (annotations regroupées dans col 33)
  
  pop_GOOD<-read.table("/home/data/SG/Good NODE.csv",header=TRUE,sep=";",stringsAsFactors = FALSE)
  pop<-data.frame(pop=pop,file_annot=1:nrow(pop_GOOD))
  pop_GOOD<-data.frame(pop_GOOD)
  corresp<-merge(pop,pop_GOOD,by.y="pop.Names",by.x="pop")
  corresp<-corresp[order(corresp$file_annot),]
  
  fcs = compensate_data(fcs)
  fcs = logiclTransformCIPHE(fcs)
  #fcs = na_data(fcs)
  
  efcs<-exprs(fcs)
  efcs[,"Annotation"]<-corresp$pop.ID[efcs[,"Annotation"]]
  efcs
}



#augmentation des pop

taille_pop_dyadem <- function(efcsALL){

  res = NULL
  Pop=efcsALL[,"Annotation"]

  K<-ncol(efcsALL)
  efcsALL<-efcsALL[,-K]
  for (i in sort(unique(Pop))){ #ds chaque pop
    inds = sample(x = 1:length(which(Pop == i)), size = 10000, replace = TRUE) #x = longueur de pop.   replace = T -> des valeurs peuvent etre reprises
    Inds = which(Pop == i) #position
    taille = Inds[inds]
    res = c(res, taille) #vecteur qui recupere resultat et taille à chq tour de boucle
  }
  res <- sample(x = res, size = length(res), replace = FALSE) #pour que la derniere population ne soit pas enlevee avec la validation_split=0.1, qui prend 10% de la derniere population (colonne)

  return(list(data = efcsALL[res,], annot = Pop[res])) #fichier avec res + pop avec res

}




#Kmeans

k_means = function(flow.frame, K = 5){
  k = kmeans(flow.frame, centers = K)
  return(k)
}


# training avec modification taille pop

train_taille_modif <- function(K=5){

  efcsB<-fct("/home/data/SG/Manual_Gated/Basal/")
  efcsC<-fct("/home/data/SG/Manual_Gated/CD3e/")
  efcsP<-fct("/home/data/SG/Manual_Gated/polyIC/")

  efcsALL<-rbind(efcsB,efcsC,efcsP)

  dataset<-taille_pop_dyadem(efcsALL = efcsALL)

  reticulate::source_python('/home/workspace/toys/dyadem_percep_JS.py')

  scalefactor<-rep(NA,6) #pour les col 1 à 6
  for (i in 1:6){
    scalefactor[i]<-quantile(dataset$data[,i],probs = 0.9)
    dataset$data[,i]<-dataset$data[,i]/scalefactor[i]
  }
  k <- k_means(dataset$data[which(dataset$annot==11),1:24],K)$cluster
  k[k!=1]=k[k!=1]+max(dataset$annot)-1
  k[k==1]=11
  dataset$annot[dataset$annot==11]=k

  mm = model.matrix(~as.factor(dataset$annot)-1)

  #n = neuralnet(dataset$data[,7:24], mm) #training : colonnes de 7 à 24 car de 1 à 6 : scatters à échelles différentes et col 25 - 26 : inutiles (flag et time)

  n = neuralnet(dataset$data[,1:24], mm) #training

  return( list(model=n,scalefactor=scalefactor, nonId = sort(unique(k))))

}



compensate_data1 <- function(flow.frame){
  comp = compensate(flow.frame,flow.frame@description[["$SPILLOVER"]])
  return(comp)
}


logiclTransformCIPHE1 <- function(flow.frame, value = NULL, marker = NULL){
  
  if(is.null(marker)){
    if(is.null(flow.frame@description[["$SPILLOVER"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["$SPILLOVER"]])
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



#concatenation des FMO

fct_FMO<-function(dir_fcs){
  
  flow.frames <- list.files(dir_fcs, pattern =".fcs", full.names = TRUE, recursive = FALSE)#liste des fichiers
  
  pop<-as.character(gsub(".fcs$","",gsub("^.*//","",flow.frames)))
  
  flow.frames <- lapply(flow.frames, function(i) {read.FCS(i)}) #infos des fichiers
  
  fcs <- concatenate.FCS.CIPHE(flow.frames, params = "Annotation") #infos + concatenation des annotations (annotations regroupées dans col 33)
  
  pop_GOOD<-read.table("/home/data/SG/Good NODE2.csv",header=TRUE,sep=";",stringsAsFactors = FALSE)
  pop<-data.frame(pop=pop,file_annot=1:nrow(pop_GOOD))
  pop_GOOD<-data.frame(pop_GOOD)
  corresp<-merge(pop,pop_GOOD,by.y="pop.Names",by.x="pop")
  corresp<-corresp[order(corresp$file_annot),]
  
  fcs = compensate_data1(fcs)
  fcs = logiclTransformCIPHE1(fcs)
  
  efcs<-exprs(fcs)
  efcs[,"Annotation"]<-corresp$pop.ID[efcs[,"Annotation"]]
  efcs
}




#classification


percept = function(efcsALL_FMO, n, seuil = 0.0){  
  
  #  scalefactor<-rep(NA,6) #pour les col 1 à 6
  for (i in 1:6){
    #   scalefactor[i]<-quantile(efcsALL_FMO[,i],probs = 0.9)
    efcsALL_FMO[,i]<-efcsALL_FMO[,i]/n[[2]][i]
  }
  
  # print(system.time(res_loc <- n[[1]][[2]](efcsALL_FMO[,1:24])))
  #   
  # center_annot = apply(res_loc$numpy(), MARGIN = 1, FUN = function(x, seuil){
  #   ii = which.max(x)
  #   ifelse(x[ii] > seuil, ii, length(x)+1)
  #  }, seuil = seuil)  #voir arrangement des annotations : le max
  #   
  
  center_annot<-n[[1]][[3]](efcsALL_FMO[,1:24])+1 # python is 0-base
  
  new_pop = n[[3]]  #31 29 30 28 37
  #cell_pop_11 = dataset$data[which(dataset$annot==11),1:24]
  center_annot[center_annot %in% new_pop]<-new_pop[1]
  
  
  return(center_annot)
}




#annotations fichiers fcs avec perceptron

annot_perceptron <- function(dir_fcs){
  
  fcss <- list.files(dir_fcs, pattern ="fcs$", full.names = TRUE, recursive = FALSE)  #chemin des fichiers à annoter
 # fcss<-dir(path = "/home/data/SG/All_FMO_Annoted_detrans_decomp/",pattern = "fcs$",full.names = TRUE)
  TT<-length(fcss)                                #les fichiers fcs
  tab<-matrix(0,27,TT)                            #27 pop
  colnames(tab)<-gsub("^.*/","",fcss[1:TT])
  rownames(tab)<-pop_GOOD$pop.Names
  
  system.time({
    for (i in 1:length(fcss[1:TT])){
      print(fcss[i])
      fcso<-fcs<-read.FCS(fcss[i])
      fcs = compensate_data(fcs)                  #compensation des data
      fcs = logiclTransformCIPHE(fcs)             #transformation des data
      ann<-percept(exprs(fcs),training)           #annotation perceptron
      tab[, i] <- table(factor(ann,levels=1:27))/sum(table(factor(ann,levels=1:27)))
      
      print(tab[,i])
      fcso<-enrich.FCS.CIPHE(original=fcso,new.column = ann,nw.names = "Pop_ann_percept")  #nouvelle colonne d'annotations
      write.FCS(fcso,filename = gsub(".fcs$","_enrich.fcs",fcss[i]))

    }
   }
   )
  #scaff = r1@exprs[,"popIDscaffold"]
 # p = table(scaff)/sum(table(scaff))
  #t_p = cbind(tab,p)
  write.table(tab,file="/home/data/SG/pop_proportions_FMO_Camille_Landmarks_Detrans_ann_percep_logTrans500.csv",sep=",")  #csv avec proportions
  # ref=table(annot)/sum(table(annot))
  # cbind(tab,ref=ref)
}

