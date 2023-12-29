from connect import *
import numpy as np
import datetime

#Script créé par Christine BOUTRY (Décembre 2023) pour calculer l'indice de complexité MCSV de chaque faisceau d'un plan de traitement sur l'Halcyon
#l'indice est calculé pour chaque MLC (Upper et Lower) et pour le faisceau complet suivant la méthode classique MCVS5
#Il faudra intégrer à l'avenir la méthode pondérée pour calculer MCSVw
#référence biblio: Tamura 2020, https://doi.org/10.1007/s13246-020-00891-2


#--------------------Fonctions-------------------#

# Récupération du nombre de lames ouvertes et dont l'écart entre les bancs droit et gauche est inférieur à 10mm

# Collimateur Upper 
def get_Na_u(cp):
    retval_u = 0.0
    for i in range(leafPairCount_u):
        aLeaf_u = cp.LeafPositions[0][i]   # position des lames du banc gauche du MLC Upper
        bLeaf_u = cp.LeafPositions[1][i]   # position des lames du banc droit du MLC Upper
        
        if abs(aLeaf_u - bLeaf_u) <= 0.07:
           continue
        
        if abs(aLeaf_u - bLeaf_u) > 1:
           continue
        
        retval_u += 1
        
    return retval_u

# Collimateur Lower 
def get_Na_l(cp):
    retval_l = 0.0
    for i in range(leafPairCount_l):
        aLeaf_l = cp.LeafPositions[2][i]   # position des lames du banc gauche du MLC Lower
        bLeaf_l = cp.LeafPositions[3][i]   # position des lames du banc droit du MLC Lower
        
        if abs(aLeaf_l - bLeaf_l) <= 0.07:
           continue
        
        if abs(aLeaf_l - bLeaf_l) > 1:
           continue
        
        retval_l += 1
        
    return retval_l





########## Calcul de l'indice SAS10 ##########


CurrBS = get_current('BeamSet')
CurrPat = get_current('Patient')

plan=get_current('Plan')


# initialisation des variables pour les collimateurs Upper (29 lames) et Lower (28 lames)
leafPairCount_u = 29
leafPairCount_l = 28

results_u = {}
results_l = {}

planMU = 0
SAS10_u = 0
SAS10_l = 0

list_SAS10_beam=[]

# Boucle pour chaque faisceau du beamSet
    
for beam in CurrBS.Beams:
    cp_module = beam.Segments

    SASj_u = 0.0
    results_u[beam.Name] = []
    prev_normalizedMUj_u = 0.0
    
    SASj_l = 0.0
    results_l[beam.Name] = []
    prev_normalizedMUj_l = 0.0

    # Boucle sur les CP de chaque faisceau du BeamSet
    for cp in cp_module:

        normalizedMUj_u = cp.RelativeWeight
        Naj_u = get_Na_u(cp)
        Naj_u = Naj_u / leafPairCount_u
        SASj_u += (normalizedMUj_u * Naj_u)
        
        normalizedMUj_l = cp.RelativeWeight
        Naj_l = get_Na_l(cp)
        Naj_l = Naj_l / leafPairCount_l
        SASj_l += (normalizedMUj_l * Naj_l)
        
    # Calcul et affichage détaillé du SAS10 pour chaque faisceau
    SAS10_beam = (SASj_u + SASj_l) / 2
    list_SAS10_beam.append(SAS10_beam)
    print (' Indices de Modulation faisceau: ',beam.Name)
    print (' SAS10_beam_u : \t\t{0:.4f}'.format(SASj_u))
    print (' SAS10_beam_l : \t\t{0:.4f}'.format(SASj_l))
    print (' SAS10_beam : \t\t{0:.4f}'.format(SAS10_beam))
    print (' -')


    planMU = planMU + beam.BeamMU

    
    SAS10_u += SASj_u * beam.BeamMU
    SAS10_l += SASj_l * beam.BeamMU



if planMU != 0 :
    SAS10_u = SAS10_u/planMU
    SAS10_l = SAS10_l/planMU
else :
    SAS10_u = 'NaN'
    SAS10_l = 'NaN'

# Affichage condensé du SAS10 pour chaque faisceau
print (' ')
print ('--------------------------------------------')
print ('ID faisceau , SAS10_beam')

i=0
for beam in CurrBS.Beams:
    SAS10_beam = list_SAS10_beam[i]
    print (beam.Name , ' ', round(SAS10_beam,4))
    i += 1    
print ('--------------------------------------------') 


# Calcul et affichage du SAS10 du BeamSet    
SAS10_plan = (SAS10_u + SAS10_l) / 2.0 

print (' ')  
print (' Indices de Modulation du beamset: ',CurrBS.DicomPlanLabel)
print (' ')
print (' SAS10_plan_u : \t\t{0:.2f}'.format(SAS10_u))
print (' SAS10_plan_l : \t\t{0:.2f}'.format(SAS10_l))
print (' SAS10_plan : \t\t{0:.2f}'.format(SAS10_plan))
print ('ID plan , SAS10_plan')
print (CurrBS.DicomPlanLabel , ' ', round(SAS10_plan,4))

