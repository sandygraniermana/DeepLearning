library("flowCore")
library("FlowCIPHE")


###################################################################################################################
                                                      
                                                      # Concated_gated_clean :

###################################################################################################################


flow.frames <- list.files("/home/workspace/SG/Concat_gated_clean/Camille/", pattern =".fcs", full.names = TRUE, recursive = FALSE)#liste des fichiers
flow.frames <- lapply(flow.frames, function(i) {read.FCS(i)}) #infos des fichiers

names(flow.frames) <- basename(list.files("/home/workspace/SG/Concat_gated_clean/Camille/", pattern = ".fcs", recursive = FALSE)) #noms des fichiers dans le chemin
fcs <- concatenate.FCS.CIPHE(flow.frames, params = "Annotation") #infos + concatenation des annotations (annotations regroupées dans col 33)

#annotations = col 33

table(fcs@exprs[,33]) #27 pop

rain <- rainbow(length(unique(fcs@exprs[,"Annotation"])))
color <- rain[fcs@exprs[,"Annotation"]]
plot(fcs@exprs[,c(31,32)],pch=".",cex=1.5, col=color) # 31 et 32 = tSNE1 et tSNE2 (red dim)

MM=model.matrix(~as.factor(fcs@exprs[,33])-1) #dummy variable pour les annotations (64e colonne)
reticulate::source_python('/home/workspace/SG/dyadem_percep.py') #lien avec script python red dim
#n = neuralnet(fcs@exprs[,c(1:24)], MM) 


#####################################################################################################################

                                                      # Manuel_Gated

#####################################################################################################################


flow.framesB <- list.files("/home/workspace/SG/Manual_Gated/Basal/", pattern =".fcs", full.names = TRUE, recursive = FALSE)#liste des fichiers
flow.framesC <- list.files("/home/workspace/SG/Manual_Gated/CD3e/", pattern =".fcs", full.names = TRUE, recursive = FALSE)#liste des fichiers
flow.framesP <- list.files("/home/workspace/SG/Manual_Gated/polyIC/", pattern =".fcs", full.names = TRUE, recursive = FALSE)#liste des fichiers

flow.framesB <- lapply(flow.framesB, function(i) {read.FCS(i)}) #infos des fichiers
flow.framesC <- lapply(flow.framesC, function(i) {read.FCS(i)}) #infos des fichiers
flow.framesP <- lapply(flow.framesP, function(i) {read.FCS(i)}) #infos des fichiers

names(flow.framesB) <- basename(list.files("/home/workspace/SG/Manual_Gated/Basal/", pattern = ".fcs", recursive = FALSE)) #noms des fichiers dans le chemin
names(flow.framesC) <- basename(list.files("/home/workspace/SG/Manual_Gated/CD3e/", pattern = ".fcs", recursive = FALSE))
names(flow.framesP) <- basename(list.files("/home/workspace/SG/Manual_Gated/polyIC/", pattern = ".fcs", recursive = FALSE))

fcsB <- concatenate.FCS.CIPHE(flow.framesB, params = "Annotation") #infos + concatenation des annotations (annotations regroupées dans col 33)
fcsC <- concatenate.FCS.CIPHE(flow.framesC, params = "Annotation")
fcsP <- concatenate.FCS.CIPHE(flow.framesP, params = "Annotation")

table(fcs@exprs[,27])


MM=model.matrix(~as.factor(fcs@exprs[,27])-1) #dummy variable pour les annotations (64e colonne)

reticulate::source_python('/home/workspace/SG/dyadem_percep.py') #lien avec script python red dim

#n = neuralnet(fcs.trans@exprs[,marker], MM)
n = neuralnet(fcs......, MM) 




# flow.framesB <- list.files("/home/workspace/SG/Manual_Gated/Basal/", pattern =".fcs", full.names = TRUE, recursive = FALSE)#liste des fichiers
# 
# pop<-gsub(".fcs$","",gsub("^.*//","",flow.framesB)) #expression reguliere : nom des pop
# 
# pop<-data.frame(pop=pop,file_annot=1:nrow(pop_GOOD)) #nom des pop  = pop et numero d'annotation (num de ligne) = file_annot
# pop_GOOD<-read.table("/home/workspace/SG/Good NODE.csv",header=TRUE,sep=";",stringsAsFactors = FALSE) #lecture du csv où y a les vrais num d'annotation
# pop_GOOD<-data.frame(pop_GOOD) #mettre sous forme d'un tableau dans R
# corresp<-merge(pop,pop_GOOD,by.y="pop.Names",by.x="pop") #pop + tab csv : merger pop Names du csv avec pop des pop(flow.framesB)
# corresp<-corresp[order(corresp$file_annot),] #dans le bon ordre : noms en fct des annotations(file_annot)
# 
# fcsB <- concatenate.FCS.CIPHE(flow.framesB, params = "Annotation") #infos + concatenation des annotations (annotations regroupées dans col 33)
# fcsC <- concatenate.FCS.CIPHE(flow.framesC, params = "Annotation")
# fcsP <- concatenate.FCS.CIPHE(flow.framesP, params = "Annotation")
# 
# # fcs1 = compensate_data(fcsB)
# # fcs1 = logiclTransformCIPHE(fcs1)
# 
# efcsB<-exprs(fcsB)
# efcsB[,"Annotation"]<-corresp$pop.ID[efcsB[,"Annotation"]]

 
