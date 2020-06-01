
######################################### commandes pour annotations perceptron #################################


# les fonctions sont dans le fichier : "/home/data/SG/annotation_perceptron.R"




training = train_taille_modif() #entraînement avec taille pop modifiées (10000 cellules) (~ 13 min et accuracy = 0.92)


pop_GOOD<-read.table("/home/data/SG/Good NODE2.csv",header=TRUE,sep=";",stringsAsFactors = FALSE) #noms des pop


annotation = annot_perceptron("/home/data/SG/FMO_Camille_Landmarks&Landmarks_Detrans_ann_percep_logTrans500") #annotation des fichiers à annoter avec perceptron

