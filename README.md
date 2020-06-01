# DeepLearning

Stage de Master 2 _ CIPHE : développement d'un perceptron de classification pour annoter des cellules issues de cytométrie en flux.

Voici quelques uns des scripts :

-- Pour des données IMPC :

- fct_perceptron.R : contient toutes les fonctions 
- train_prediction.R : contient toutes les commandes
- neuralnet.py : contient l'entraînement

-- Pour les données DYADEM :

 #prédiction + tests de performance :
 
  ###Entraînement avec un fichier déjà annoté avec une étape de preprocessing : augmentation de chaque population à 10000 cellules après concaténation des 3 conditions inflammatoires + clustering sur la population des cellules non identifiées. Puis prédiction des nnotations cellulaires des fichiers DYADEM 2018 et 2019 puis les comparer :
- dyadem_fct.R : contient toutes les fonctions
- dyadem_2018.R et dyadem_2019 : contient toutes les commandes
- dyadem_percep.py : contient l'entraînement 

  ###Entraînement avec un fichier déjà annoté avec une étape de préprocessing : augmentation de chaque population à 3000 cellules par condition inflammatoire + clustering sur la population des cellules non identifiées. Puis prédire les annotations cellulaires des fichiers DYADEM 2018 et 2019 puis les comparer :
- annotation_percept_3000cell.R : contient toutes les fonctions
- cmd_annot_percept_stats_3000cell.R : contient toutes les commandes
- dyadem_percep.py : contient l'entraînement 

 #prédiction + nouvelle colonne d'annotation dans des FCS :

- annotation_perceptron.R : contient toutes les fonctions
- cmd_annot_percept.R : contient toutes les commandes
- dyadem_percep.py : contient l'entraînement 

