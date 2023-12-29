import os, sys, System, clr, random, math
from connect import *


#Script créé par Christine BOUTRY (Décembre 2023) pour calculer les indices de complexité PA, PI et PM de chaque faisceau d'un plan de traitement sur l'Halcyon
#les indices sont calculés pour chaque MLC (Upper et Lower) et pour le faisceau complet suivant la méthode classique PA5, PI5 et PM5
#Il faudra intégrer à l'avenir la méthode pondérée pour calculer MCSVw
#référence biblio: Tamura 2020, https://doi.org/10.1007/s13246-020-00891-2


            



#--------------------------- Functions ---------------------------------------#

# initialisation des variables pour les collimateurs Upper (29 lames) et Lower (28 lames)
leafPairCount_u = 29
#leafThicknessesInMM_u = []

leafPairCount_l = 28
#leafThicknessesInMM_l = []

max_positions_A_u = []
min_positions_B_u = []


# Calcul de l'union des aires définies par les lames de chaque collimateur (U(AAij))
# Dénominateur de l'indice BMi

def get_area_union(beam):
    area_u = 0.0

# Initialisation des positions max des lames ( ouverture max du champ=28*28cm²)
    max_positions_A_u = [-14.0] * len(leafThicknessesInMM_u)
    min_positions_B_u = [+14.0] * len(leafThicknessesInMM_u)
    area_l = 0.0
    max_positions_A_l = [-14.0] * len(leafThicknessesInMM_l)
    min_positions_B_l = [+14.0] * len(leafThicknessesInMM_l)


    cp_module = beam.Segments

# recherche de la position max de chaque lame sur l'ensemble des CP    
    for cp in cp_module:
        aLeaf_u = list(cp.LeafPositions[1])   # position des lames du banc droit du MLC Upper
        bLeaf_u = list(cp.LeafPositions[0])   # position des lames du banc gauche du MLC Upper
        for j in range(len(leafThicknessesInMM_u)):
            max_positions_A_u[j] = max(aLeaf_u[j], max_positions_A_u[j])
            min_positions_B_u[j] = min(bLeaf_u[j], min_positions_B_u[j])
        
        
        aLeaf_l = list(cp.LeafPositions[3])   # position des lames du banc droit du MLC Lower
        bLeaf_l = list(cp.LeafPositions[2])   # position des lames du banc gauche du MLC Lower
        for j in range(len(leafThicknessesInMM_l)):
            max_positions_A_l[j] = max(aLeaf_l[j], max_positions_A_l[j])
            min_positions_B_l[j] = min(bLeaf_l[j], min_positions_B_l[j])
                
        
# Calcul de l'union des aires pour le collimateur Upper en ne considérant que les lames ouvertes
    for i in range(len(max_positions_A_u)):
        aLeaf_u = max_positions_A_u[i]
        bLeaf_u = min_positions_B_u[i]
        
        if abs(aLeaf_u - bLeaf_u) <= 0.07:
           continue
           
        area_u += abs(aLeaf_u - bLeaf_u) * leafThicknessesInMM_u[i]


# Calcul de l'union des aires pour le collimateur Lower en ne considérant que les lames ouvertes
    for i in range(len(max_positions_A_l)):
        aLeaf_l = max_positions_A_l[i]
        bLeaf_l = min_positions_B_l[i]
        
        if abs(aLeaf_l - bLeaf_l) <= 0.07:
           continue
           
        area_l += abs(aLeaf_l - bLeaf_l) * leafThicknessesInMM_l[i]


    return area_u, area_l
    

    
# Calcul de l'aire définie par les lames pour chaque CP de chaque collimateur (AAij)
# Dénominateur de l'indice AIij, numérateur des indices BAi et BMi
def get_area(cp):

# Collimateur Upper
    retval_u = 0.0
    for i in range(leafPairCount_u):
        aLeaf_u = cp.LeafPositions[1][i]
        bLeaf_u = cp.LeafPositions[0][i]
        retval_u += leafThicknessesInMM_u[i] * abs(aLeaf_u - bLeaf_u)

# Collimateur Lower
    retval_l = 0.0
    for i in range(leafPairCount_l):
        aLeaf_l = cp.LeafPositions[3][i]
        bLeaf_l = cp.LeafPositions[2][i]
        retval_l += leafThicknessesInMM_l[i] * abs(aLeaf_l - bLeaf_l)


    return retval_u, retval_l


    
# Calcul du périmètre irradié défini par les lames pour chaque CP de chaque collimateur (APij)
# numérateur de l'indice AIij    
def get_perimeter(cp):

# Collimateur Upper - on ne garde que les lames ouvertes
    perimeter_u = 0.0
    for n in range(leafPairCount_u):
        an_u = cp.LeafPositions[1][n]   # position des lames du banc droit du MLC Upper
        bn_u = cp.LeafPositions[0][n]   # position des lames du banc gauche du MLC Upper

        if n == 0 and an_u - bn_u > 0.07:
            perimeter_u += abs(an_u - bn_u)

        if n == leafPairCount_u - 1:
            perimeter_u += abs(an_u - bn_u)
            continue

        ap_u = cp.LeafPositions[1][n + 1]
        bp_u = cp.LeafPositions[0][n + 1]
        
        if abs(an_u-bn_u) <= 0.07:
           continue

        adj_pos_diff_B_u = bn_u - bp_u
        adj_pos_diff_A_u = an_u - ap_u

        if adj_pos_diff_B_u < 0.0:
            perimeter_u += min(abs(adj_pos_diff_B_u), abs(bn_u - an_u))
        else:
            perimeter_u += min(abs(adj_pos_diff_B_u), abs(bp_u - ap_u))

        if adj_pos_diff_A_u < 0.0:
            perimeter_u += min(abs(adj_pos_diff_A_u), abs(an_u - bn_u))
        else:
            perimeter_u += min(abs(adj_pos_diff_A_u), abs(ap_u - bp_u))


        perimeter_u += 2.0 * leafThicknessesInMM_u[n]

# Collimateur Lower - on ne garde que les lames ouvertes
    perimeter_l = 0.0
    for n in range(leafPairCount_l):
        an_l = cp.LeafPositions[3][n]   # position des lames du banc droit du MLC Lower
        bn_l = cp.LeafPositions[2][n]   # position des lames du banc gauche du MLC Lower

        if n == 0 and an_l - bn_l > 0.07:
            perimeter_l += abs(an_l - bn_l)

        if n == leafPairCount_l - 1:
            perimeter_l += abs(an_l - bn_l)
            continue

        ap_l = cp.LeafPositions[3][n + 1]
        bp_l = cp.LeafPositions[2][n + 1]
        
        if abs(an_l-bn_l) <= 0.07:
           continue

        adj_pos_diff_B_l = bn_l - bp_l
        adj_pos_diff_A_l = an_l - ap_l

        if adj_pos_diff_B_l < 0.0:
            perimeter_l += min(abs(adj_pos_diff_B_l), abs(bn_l - an_l))
        else:
            perimeter_l += min(abs(adj_pos_diff_B_l), abs(bp_l - ap_l))

        if adj_pos_diff_A_l < 0.0:
            perimeter_l += min(abs(adj_pos_diff_A_l), abs(an_l - bn_l))
        else:
            perimeter_l += min(abs(adj_pos_diff_A_l), abs(ap_l - bp_l))


        perimeter_l += 2.0 * leafThicknessesInMM_l[n]


    return perimeter_u, perimeter_l




#--------------------- Call of funtions --------------------------------------#



#    plan = self.context.PlanSetup

CurrBS = get_current('BeamSet')
CurrPat = get_current('Patient')

plan=get_current('Plan')
results_u = {}
results_l = {}

planMU = 0
PA_u = 0
PA_l = 0
PI_u = 0
PI_l = 0
PM_u = 0
PM_l = 0

list_BA_beam = []
list_BI_beam = []
list_BM_beam = []
list_BA_u = []
list_BI_u = []
list_BM_u = []
list_BA_l = []
list_BI_l = []
list_BM_l = []
    
for beam in CurrBS.Beams:
    leafThicknessesInMM_u = [1] * 29  # 29 lames de 1cm d'épaisseur
    leafThicknessesInMM_l = [1] * 28  # 28 lames de 1cm d'épaisseur


    cp_module = beam.Segments
        

    BI_u = 0.0
    BA_u = 0.0
    BM_u = 0.0



    prev_normalizedMUj_u = 0.0
    total_area_u = 0.0
    total_square_area_u = 0.0
    
    BI_l = 0.0
    BA_l = 0.0
    BM_l = 0.0



    prev_normalizedMUj_l = 0.0
    total_area_l = 0.0
    total_square_area_l = 0.0
    



    for cp in cp_module:

# récupération du poids relatif de chaque CP (MUij/MUi)
        normalizedMUj = cp.RelativeWeight

# Calcul des aires irradiées pour chaque CP de chaque collimateur
        AAj_u,AAj_l = get_area(cp)

        total_area_u += AAj_u
        total_square_area_u += AAj_u ** 2
        BA_u += (normalizedMUj * AAj_u)

        total_area_l += AAj_l
        total_square_area_l += AAj_l ** 2
        BA_l += (normalizedMUj * AAj_l)

# Calcul des périmètres irradiés pour chaque CP de chaque collimateur
        perimeter_u, perimeter_l = get_perimeter(cp)
        BI_u += (normalizedMUj * (perimeter_u ** 2 / (4 * math.pi * AAj_u)))
        BI_l += (normalizedMUj * (perimeter_l ** 2 / (4 * math.pi * AAj_l)))


# Calcul de l'union des aires irradiées pour chaque faisceau de chaque collimateur
    areaunion_u, areaunion_l = get_area_union(beam)
    


    avgArea_u = total_area_u / len(cp_module)
    stdArea_u = math.sqrt(total_square_area_u / len(cp_module) - avgArea_u ** 2)
    BM_u = 1 - BA_u / areaunion_u
    
    avgArea_l = total_area_l / len(cp_module)
    stdArea_l = math.sqrt(total_square_area_l / len(cp_module) - avgArea_l ** 2)
    BM_l = 1 - BA_l / areaunion_l
    

# Calcul et affichage des indices BA, BI et BM pour chaque faisceau   
    BA_beam = (BA_u + BA_l) / 2.0
    BI_beam = (BI_u + BI_l) / 2.0
    BM_beam = (BM_u + BM_l) / 2.0
    list_BA_beam.append(BA_beam)
    list_BI_beam.append(BI_beam)
    list_BM_beam.append(BM_beam)
    list_BA_u.append(BA_u)
    list_BI_u.append(BI_u)
    list_BM_u.append(BM_u)
    list_BA_l.append(BA_l)
    list_BI_l.append(BI_l)
    list_BM_l.append(BM_l)
    print (' Indices de Modulation faisceau: ',beam.Name)
    print (' BA_u : \t\t{0:.4f}'.format(BA_u), "  ", ' BA_l : \t\t{0:.4f}'.format(BA_l))
    print (' BI_u : \t\t{0:.4f}'.format(BI_u), "  ", ' BI_l : \t\t{0:.4f}'.format(BI_l))
    print (' BM_u : \t\t{0:.4f}'.format(BM_u), "  ", ' BM_l : \t\t{0:.4f}'.format(BM_l))
    print (' ')
    print ('ID faisceau , BA_beam , BA_beam , BM_beam')
    print (beam.Name , ' ', round(BA_beam,4) , ' ', round(BI_beam,4) , ' ', round(BM_beam,4))
    print (' -')

# Affichage condensé des indices BA, BI et BM pour chque faisceau  
print ('-----------------------------------------------------------')
print ('ID faisceau , BA_beam , BI_beam , BM_beam')
i=0
for beam in CurrBS.Beams:
    BA_beam = list_BA_beam[i]
    BI_beam = list_BI_beam[i]
    BM_beam = list_BM_beam[i]

    print (beam.Name , ' ', round(BA_beam,4) , ' ', round(BI_beam,4) , ' ', round(BM_beam,4))
    i += 1
    
print ('-----------------------------------------------------------')


# Calcul et affichage des indices PA, PI et PM pour l'ensemble du BeamSet sélectionné  
i=0
for beam in CurrBS.Beams:
    BA_u = list_BA_u[i]
    BI_u = list_BI_u[i]
    BM_u = list_BM_u[i]
    BA_l = list_BA_l[i]
    BI_l = list_BI_l[i]
    BM_l = list_BM_l[i]
    

    planMU = planMU + beam.BeamMU
    PA_u += BA_u * beam.BeamMU
    PA_l += BA_l * beam.BeamMU
    
    PI_u += BI_u * beam.BeamMU
    PI_l += BI_l * beam.BeamMU
    
    PM_u += BM_u * beam.BeamMU
    PM_l += BM_l * beam.BeamMU

    i += 1

if planMU != 0 :
    PA_u = PA_u/planMU
    PI_u = PI_u/planMU
    PM_u = PM_u/planMU
    
    PA_l = PA_l/planMU
    PI_l = PI_l/planMU
    PM_l = PM_l/planMU
    
else :
    PA_u = 'NaN'
    PI_u = 'NaN'
    PM_u = 'NaN'
    
    PA_l = 'NaN'
    PI_l = 'NaN'
    PM_l = 'NaN'
    
PI_beamSet = (PI_u + PI_l) / 2.0
PA_beamSet = (PA_u + PA_l) / 2.0
PM_beamSet = (PM_u + PM_l) / 2.0

print (' ')
print (' Indices de Modulation du beamset: ',CurrBS.DicomPlanLabel)
print (' ')
print (' PA_u : \t\t{0:.2f}'.format(PA_u), "  ", ' PA_l : \t\t{0:.2f}'.format(PA_l))
print (' PI_u : \t\t{0:.2f}'.format(PI_u), "  ", ' PI_l : \t\t{0:.2f}'.format(PI_l))
print (' PM_u : \t\t{0:.2f}'.format(PM_u), "  ", ' PM_l : \t\t{0:.2f}'.format(PM_l))
print ('ID BeamSet , PA_beamset , PI_beamset , PM_beamset')
print (CurrBS.DicomPlanLabel , ' ', round(PA_beamSet,4) , ' ', round(PI_beamSet,4) , ' ', round(PM_beamSet,4))
