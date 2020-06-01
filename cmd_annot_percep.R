
######################################### commandes pour annotations perceptron #################################


# les fonctions sont dans le fichier : "/home/data/SG/annotation_perceptron.R"




training = train_taille_modif() #entraînement avec taille pop modifiées (10000 cellules) (~ 13 min et accuracy = 0.92)


pop_GOOD<-read.table("/home/data/SG/Good NODE2.csv",header=TRUE,sep=";",stringsAsFactors = FALSE) #noms des pop


annotation = annot_perceptron("/home/data/SG/FMO_Camille_Landmarks&Landmarks_Detrans_ann_percep_logTrans500") #annotation des fichiers à annoter avec perceptron

#l = read.FCS("/home/data/SG/DYADEM_Data_annoted_scaffold_ann_percep/enriched_NotCST_CD008b.2_53-5.8_IgG1k_Rat_basal_19-JUL-2018_1208_mice#41_filter-N.fcs_enrich.fcs")
# r = read.FCS("/home/data/SG/All_FMO_Annoted_detrans_decomp_ok/enriched_NotCST_FMO_FMO_FMO_FMO_basal_03-OCT-2018_4325_mice#99_filter-y.fcs_enrich.fcs")
# 
# R = exprs(r)
# 
# a = r@exprs[,28]
# scaff = r@exprs[,"popIDscaffold"]
# 
# 
# table(r@exprs[,28])  # Pop_ann_percept
# table(r@exprs[,27])  # popIDscaffold
# 

# r1 = read.FCS("/home/data/SG/All_FMO_Annoted_detrans_decomp_pres/enriched_NotCST_FMO_FMO_FMO_FMO_polyIC_29-JUN-2018_871_mice#29_filter-N.fcs.fcs")
# scaff = r1@exprs[,"popIDscaffold"]
# p = table(scaff)/sum(table(scaff))
# n = pop_GOOD$pop.Names
# c = cbind(n,p)
# write.table(c,file="/home/data/SG/pop_proportions_ess4.csv",sep=",")


#n[[2]]$save("n_10000cell.h5")
