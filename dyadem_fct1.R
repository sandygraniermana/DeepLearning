

library("flowCore")
library("FlowCIPHE")

########################################### traitement des data #############################################################

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




#entrainement

train <- function(efcs){
  
  efcsB<-fct("/home/data/SG/Manual_Gated/Basal/")
  efcsC<-fct("/home/data/SG/Manual_Gated/CD3e/")
  efcsP<-fct("/home/data/SG/Manual_Gated/polyIC/")
  
  efcsALL<-rbind(efcsB,efcsC,efcsP)
  
 # dataset<-taille_pop_dyadem(efcsALL = efcsALL)

  
  #mm = model.matrix(~as.factor(dataset$annot)-1)  
  MM = model.matrix(~as.factor(efcsALL[,"Annotation"])-1)
  
  reticulate::source_python('/home/data/SG/dyadem_percep.py')
  
  #n = neuralnet(dataset$data, mm) #training
  n = neuralnet(efcsALL, MM)
  return(n)
  
}


#Kmeans

k_means = function(flow.frame, K = 5){
  k = kmeans(flow.frame, centers = K)
  return(k)
}


# avec modification taille pop

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


# efcs_FMO_B<-fct_FMO("/home/workspace/SG/Concat_gated_clean/FMO_2019_populations cells/Basal/")
# efcs_FMO_C<-fct_FMO("/home/workspace/SG/Concat_gated_clean/FMO_2019_populations cells/CD3e/")
# efcs_FMO_P<-fct_FMO("/home/workspace/SG/Concat_gated_clean/FMO_2019_populations cells/PolyIC/")
# 
# efcsALL_FMO = rbind(efcs_FMO_B, efcs_FMO_C, efcs_FMO_P)


#classification

# percept = function(efcsALL_FMO, n, seuil = 0.0){  
#   
# #  scalefactor<-rep(NA,6) #pour les col 1 à 6
#   for (i in 1:6){
#  #   scalefactor[i]<-quantile(efcsALL_FMO[,i],probs = 0.9)
#     efcsALL_FMO[,i]<-efcsALL_FMO[,i]/n[[2]][i]
#   }
#   
#   res_loc = n[[1]][[2]](efcsALL_FMO[,1:24])
#     
#   center_annot = apply(res_loc$numpy(), MARGIN = 1, FUN = function(x, seuil){
#     ii = which.max(x)
#     ifelse(x[ii] > seuil, ii, length(x)+1)
#    }, seuil = seuil)  #voir arrangement des annotations : le max
#     
#   
#   new_pop = n[[3]]  #31 29 30 28 37
#   #cell_pop_11 = dataset$data[which(dataset$annot==11),1:24]
#   center_annot[center_annot %in% new_pop]<-new_pop[1]
#   
#   
#   return(center_annot)
# }


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





#visualisation des annotations :

#table(efcsALL_FMO[,"Annotation"]) = 27 mais il manque pop 19, 21, 22, 23
#table(efcsALL_FMO[,"popIDscaffold"]) = 27

vis_ann_annotation <- function(ann){ 
  
  annot = efcsALL_FMO[,"Annotation"]
  #annot = efcsALL_FMO[,"popIDscaffold"]
  res = list()
  col1 = diag(table(factor(ann, levels = 1:27), annot))/colSums(table(factor(ann, levels = 1:27), annot))
  col2 = colSums(table(factor(ann, levels = 1:27), annot))
  col3 = rowSums(table(factor(ann, levels = 1:27), annot))
  col4 = table(factor(ann,levels = 1:27), annot)[,11]
  res = data.frame(precision = col1, taille = col2, taille_annot = col3, nonIdentifiees = col4, prop_nonIdentifiees = col4/col3)
  
  return(res)
}


# vis_annotation_scaffold <- function(){ 
#   
#   annotation = efcsALL_FMO[,"Annotation"]
#   scaffold = efcsALL_FMO[,"popIDscaffold"]
#   res = list()
#   col1 = diag(table(factor(annotation, levels = 1:27), scaffold))/colSums(table(factor(annotation, levels = 1:27), scaffold))
#   col2 = colSums(table(factor(annotation, levels = 1:27), scaffold))
#   col3 = rowSums(table(factor(annotation, levels = 1:27), scaffold))
#   col4 = table(factor(annotation,levels = 1:27), scaffold)[,11]
#   res = data.frame(precision = col1, taille = col2, taille_annot = col3, nonIdentifiees = col4, prop_nonIdentifiees = col4/col3)
#   
#   return(res)
# }



vis <- function(annot,a){ 
  
  #annot = efcsALL_FMO[,"Annotation"]
  #a = efcsALL_FMO[,"popIDscaffold"]
  res = list()
  col1 = diag(table(factor(a, levels = 1:27), annot))/colSums(table(factor(a, levels = 1:27), annot))
  col2 = colSums(table(factor(a, levels = 1:27), annot))
  col3 = rowSums(table(factor(a, levels = 1:27), annot))
  col4 = table(factor(a,levels = 1:27), annot)[,11]
  res = data.frame(precision = col1, taille = col2, taille_annot = col3, nonIdentifiees = col4, prop_nonIdentifiees = col4/col3)
  
  return(res)
}



precision_recall <- function(ann, a){
  res = list()
  col_names = pop_GOOD$pop.Names
  col_precision = diag(table(factor(ann, levels = 1:27), a))/rowSums(table(factor(ann, levels = 1:27), a))
  col_recall = diag(table(factor(ann, levels = 1:27), a))/colSums(table(factor(ann, levels = 1:27), a))
  res = data.frame(pop_names = col_names, precision = col_precision, recall = col_recall)

  return(res)
}



precision_recall_ann_annot <- function(ann, annot){
  res = list()
  col_names = pop_GOOD$pop.Names
  col_precision = diag(table(factor(ann, levels = 1:27), annot))/rowSums(table(factor(ann, levels = 1:27), annot))
  col_recall = diag(table(factor(ann, levels = 1:27), annot))/colSums(table(factor(ann, levels = 1:27), annot))
  res = data.frame(pop_names = col_names, precision = col_precision, recall = col_recall)
  
  return(res)
}


# precision_recall_annot_a <- function(annot, a){
#   res = list()
#   col_names = pop_GOOD$pop.Names
#   col_precision = diag(table(factor(annot, levels = 1:27), a))/colSums(table(factor(annot, levels = 1:27), a))
#   col_recall = diag(table(factor(annot, levels = 1:27), a))/rowSums(table(factor(annot, levels = 1:27), a))
#   res = data.frame(pop_names = col_names, precision = col_precision, recall = col_recall)
#   
#   return(res)
# }


# precision_recall_ann_a <- function(ann, a){
#   col_names = pop_GOOD$pop.Names
#   tab = table(ann, a)
#   for (colonne in tab){
#     Tab = tab[,colonne]
#     pop_present = row.names(Tab)
#     pop_present = as.numeric(pop_present)
#     pop_n = Tab
#     i = which(pop_present == colonne)
#     col_recall = pop_n[i]/sum(pop_present[,colonne])
#     col_precision = tab[i,colonne] / sum(tab[i,])
#   }
#   
#   return(list(pop_names = col_names, recall = col_recall, precision = col_precision))
# }


precision_recall_ann_a_E <- function(ann, a){
  v_names = pop_GOOD$pop.Names
  tab = table(ann, a)
  v_precision=v_recall=rep(NA,ncol(tab))
  res = list()
  for (j in 1:ncol(tab)){
    pop_present = row.names(tab)
    pop_present = as.numeric(pop_present)
    pop_n = tab[,j] #colonne
    i = which(pop_present == j) #nom ligne = nom col
    if (length(pop_n[i]/sum(tab[, j]))>0) v_recall[j] <- pop_n[i]/sum(tab[, j])
 #   v_recall[j] = pop_n[i]/sum(tab[,j]) #pop_n[i] = diag = num col correspondant au num ligne. sum(pop_present[,j]) = sum de la col j
   # v_precision[j] = tab[i,j] / sum(tab[i,])
    if (length(tab[i,j] / sum(tab[i,]))>0) v_precision[j] <-tab[i,j] / sum(tab[i,])
   #v_names[j] = pop_GOOD$pop.Names
    res = data.frame(recall = v_recall, precision = v_precision)
  }
  return(res)
}




Compute.F1 <- function(ref.vec, run.vec)
{
  ref.IDS <- as.integer(names(table(ref.vec))[rev(order(table(ref.vec)))])
  run.IDS <- as.integer(names(table(run.vec))[rev(order(table(run.vec)))])
  
  F1.score.mat <- function(ref.vec, run.vec)
  {
    ref.IDS <- as.integer(names(table(ref.vec))[rev(order(table(ref.vec)))])
    run.IDS <- as.integer(names(table(run.vec))[rev(order(table(run.vec)))])
    f1.mat <- matrix(nrow=length(run.IDS),ncol=length(ref.IDS))
    
    for(i in 1:length(run.IDS))
    {
      for(j in 1:length(ref.IDS))
      {
        current.cluster.ids <- unlist(unlist(which(run.vec==run.IDS[i])))
        current.ref.ids <- unlist(unlist(which(ref.vec==ref.IDS[j])))
        TP <- length(intersect(current.cluster.ids, current.ref.ids))
        FP <- length(current.cluster.ids) - length(intersect(current.cluster.ids, current.ref.ids))
        FN <- length(current.ref.ids) - length(intersect(current.cluster.ids, current.ref.ids))
        
        f1.mat[i,j] <- 2*TP / (2*TP + FN + FP)
      }
    }
    return(f1.mat)
  }
  
  f1.mat <- F1.score.mat(ref.vec, run.vec)
  
  
  f1.score <- matrix(0, nrow=ncol(f1.mat), ncol=3)
  colnames(f1.score) <- c("Real Annotation", "Tested Annotation", "Score")
  f1.score[,1] <- 1:ncol(f1.mat)
  ref.list <- 1:ncol(f1.mat)
  run.list <- 1:nrow(f1.mat)
  for(i in 1:nrow(f1.score))
  {
    if(!is.matrix(f1.mat))
    {
      f1.mat <- matrix(f1.mat, ncol=length(ref.list))
    }
    
    score.id <- which.max(f1.mat)
    run.id <- score.id%%nrow(f1.mat)
    if(run.id==0)
    {
      run.id <- nrow(f1.mat)
    }
    run.id <- as.integer(run.IDS[ run.list[[run.id]] ])
    ref.id <- as.integer(ref.IDS[ ref.list[[ceiling(score.id / nrow(f1.mat))]] ])
    tmp.score <- f1.mat[score.id]
    
    f1.score[ref.id,] <- c(ref.id, run.id, tmp.score)
    
    if(nrow(f1.mat)>1 && ncol(f1.mat)>1)
    {
      run.id <- score.id%%nrow(f1.mat)
      if(run.id==0)
      {
        run.id <- nrow(f1.mat)
      }
      
      ref.list <- ref.list[-ceiling(score.id / nrow(f1.mat))]
      
      run.list <- run.list[-run.id]
      
      f1.mat <- f1.mat[,-ceiling(score.id / nrow(f1.mat))]
      if(is.matrix(f1.mat))
      {
        f1.mat <- f1.mat[-run.id,]
      }
    }
    else
    {
      i <- nrow(f1.score)
    }
  }
  f1.score[duplicated(f1.score[,2]),c(2,3)] <- c(0,0) 
  f1.score <- data.frame(f1.score,
                         stringsAsFactors = F)
  
  return(f1.score)
}


pre_rec_cellAnnot_F1_ann_annot <- function(ann, annot){
  res = list()
  col_names = pop_GOOD$pop.Names
  col_precision = diag(table(factor(ann, levels = 1:27), annot))/rowSums(table(factor(ann, levels = 1:27), annot))
  col_recall = diag(table(factor(ann, levels = 1:27), annot))/colSums(table(factor(ann, levels = 1:27), annot))
  col2 = colSums(table(factor(ann, levels = 1:27), annot))                                                         #nombre de cellules annotées manuellement
  col3 = rowSums(table(factor(ann, levels = 1:27), annot))                                                         #nombre de cellules annotées prédites
  col_F1 = Compute.F1(ref.vec = annot, run.vec = ann)
  moy_harmo = 2*((col_precision*col_recall)/(col_precision+col_recall))
  mean_r = mean(col_recall)
  mean_p = mean(col_precision)
  res = data.frame(pop_names = col_names, cell_annot_manu = col2, cell_annot_pred = col3, precision = col_precision, recall = col_recall, score_F1 = col_F1, moyenne_harmonique = moy_harmo, mean_precision_sans_nonId = mean_p, mean_recall_sans_nonId = mean_r)
  
  return(res)
}



pre_rec_cellAnnot_ann_a <- function(ann, a){
  v_names = pop_GOOD$pop.Names
  tab = table(ann, a)
  v_precision=v_recall=cell_pred=cell_scaffold=moy_harmo=rep(NA,ncol(tab))
  res = list()
  for (j in 1:ncol(tab)){
    pop_present = row.names(tab)
    pop_present = as.numeric(pop_present)
    pop_n = tab[,j] #colonne
    i = which(pop_present == j) #nom ligne = nom col
    cell_pred[j] <- sum(tab[i,])
    cell_scaffold[j] <- sum(tab[,j])
    if (length(pop_n[i]/sum(tab[,j]))>0) v_recall[j] <- pop_n[i]/sum(tab[, j])
    # v_recall[j] = pop_n[i]/sum(tab[,j]) #pop_n[i] = diag = num col correspondant au num ligne. sum(pop_present[,j]) = sum de la col j
    # v_precisifon[j] = tab[i,j] / sum(tab[i,])
    if (length(tab[i,j] / sum(tab[i,]))>0) v_precision[j] <-tab[i,j] / sum(tab[i,])
    v_names[j] = pop_GOOD$pop.Names[j]
    #F1[j] = Compute.F1(ref.vec = tab[,j], run.vec = tab[i,])
    moy_harmo[j] = 2*((v_precision[j]*v_recall[j])/(v_precision[j]+v_recall[j]))
    mean_r = mean(v_recall[j!=11])
    mean_p = mean(v_precision[j!=11])
    res = data.frame(pop_names = v_names,cell_annot_scaffold = cell_scaffold, 
                     cell_annot_pred = cell_pred, recall = v_recall,
                     precision = v_precision, moyenne_harmonique = moy_harmo, 
                     mean_precision_ = mean_p, mean_recall = mean_r)
  }
  return(res)
}



pre_rec_cellAnnot_ann_annot <- function(ann, annot){
  v_names = pop_GOOD$pop.Names
  tab = table(ann, annot)
  v_precision=v_recall=cell_pred=cell_scaffold=moy_harmo=rep(NA,ncol(tab))
  res = list()
  for (j in 1:ncol(tab)){
    pop_present = row.names(tab)
    pop_present = as.numeric(pop_present)
    pop_n = tab[,j] #colonne
    i = which(pop_present == j) #nom ligne = nom col
    cell_pred[j] <- sum(tab[i,])
    cell_scaffold[j] <- sum(tab[,j])
    if (length(pop_n[i]/sum(tab[, j]))>0) v_recall[j] <- pop_n[i]/sum(tab[, j])
    # v_recall[j] = pop_n[i]/sum(tab[,j]) #pop_n[i] = diag = num col correspondant au num ligne. sum(pop_present[,j]) = sum de la col j
    # v_precision[j] = tab[i,j] / sum(tab[i,])
    if (length(tab[i,j] / sum(tab[i,]))>0) v_precision[j] <-tab[i,j] / sum(tab[i,])
    v_names[j] = pop_GOOD$pop.Names[j]
    #F1[j] = Compute.F1(ref.vec = tab[,j], run.vec = tab[i,])
    moy_harmo[j] = 2*((v_precision[j]*v_recall[j])/(v_precision[j]+v_recall[j]))
    # mean_r = mean(v_recall[j!=11])
    # mean_p = mean(v_precision[j!=11])
    res = data.frame(pop_names = v_names,cell_annot_scaffold = cell_scaffold, cell_annot_pred = cell_pred, recall = v_recall, precision = v_precision, moyenne_harmonique = moy_harmo)
  }
  return(res)
}


pre_rec_cellAnnot_a_annot <- function(a, annot){
  v_names = pop_GOOD$pop.Names
  tab = table(a, annot)
  v_precision=v_recall=cell_pred=cell_scaffold=moy_harmo=rep(NA,ncol(tab))
  res = list()
  for (j in 1:ncol(tab)){
    pop_present = row.names(tab)
    pop_present = as.numeric(pop_present)
    pop_n = tab[,j] #colonne
    i = which(pop_present == j) #nom ligne = nom col
    cell_pred[j] <- sum(tab[i,])
    cell_scaffold[j] <- sum(tab[,j])
    if (length(pop_n[i]/sum(tab[, j]))>0) v_recall[j] <- pop_n[i]/sum(tab[, j])
    # v_recall[j] = pop_n[i]/sum(tab[,j]) #pop_n[i] = diag = num col correspondant au num ligne. sum(pop_present[,j]) = sum de la col j
    # v_precision[j] = tab[i,j] / sum(tab[i,])
    if (length(tab[i,j] / sum(tab[i,]))>0) v_precision[j] <-tab[i,j] / sum(tab[i,])
    v_names[j] = pop_GOOD$pop.Names[j]
    #F1[j] = Compute.F1(ref.vec = tab[,j], run.vec = tab[i,])
    moy_harmo[j] = 2*((v_precision[j]*v_recall[j])/(v_precision[j]+v_recall[j]))
    # mean_r = mean(v_recall[j!=11])
    # mean_p = mean(v_precision[j!=11])
    res = data.frame(pop_names = v_names,cell_annot_scaffold = cell_scaffold, cell_annot_pred = cell_pred, recall = v_recall, precision = v_precision, moyenne_harmonique = moy_harmo)
  }
  return(res)
}

