

#train(efcs) #taille pop non modifiée

training = train_taille_modif() #taille pop modifiée (10000 cellules)


efcs_FMO_B<-fct_FMO("/home/data/SG/Concat_gated_clean/FMO_2019_populations cells/Basal/")
efcs_FMO_C<-fct_FMO("/home/data/SG/Concat_gated_clean/FMO_2019_populations cells/CD3e/")
efcs_FMO_P<-fct_FMO("/home/data/SG/Concat_gated_clean/FMO_2019_populations cells/PolyIC/")
 
efcsALL_FMO = rbind(efcs_FMO_B, efcs_FMO_C, efcs_FMO_P)



ann = percept(efcsALL_FMO=efcsALL_FMO, n=training, seuil = 0.0)


visualisation = vis_ann_annotation(ann) # prediction vs manuelle


annot = efcsALL_FMO[,"Annotation"]
a = efcsALL_FMO[,"popIDscaffold"]


visua = vis(annot,a) # scaffold vs manuelle


#pre_rec = precision_recall(ann, a)

#pre_rec1 = precision_recall_ann_a(ann, a)

pre_recE = precision_recall_ann_a_E(ann, a)

pre_rec2 = precision_recall_ann_annot(ann, annot)

pre_rec3 = precision_recall_ann_annot(annot, a)


ref.vec = a
run.vec = ann

F1 = Compute.F1(ref.vec, run.vec)

prcf = pre_rec_cellAnnot_F1_ann_annot(ann, annot)

prcf2 = pre_rec_cellAnnot_ann_a(ann,a)


prcf3 = pre_rec_cellAnnot_ann_annot(ann,annot)

prcf4 = pre_rec_cellAnnot_a_annot(a,annot)

c = colMeans(prcf3[-11,4:5])
c1 = colMeans(prcf4[-11,4:5])
