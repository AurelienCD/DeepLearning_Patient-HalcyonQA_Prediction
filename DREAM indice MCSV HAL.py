from connect import *
#from os.path import join 

#from System.Windows import Window
#clr.AddReference("System.Windows.Forms")
#from System.Windows.Forms import MessageBox

#Script créé par Christine BOUTRY (Décembre 2023) pour calculer l'indice de complexité MCSV de chaque faisceau d'un plan de traitement sur l'Halcyon
#l'indice est calculé pour chaque MLC (Upper et Lower) et pour le faisceau complet suivant la méthode classique MCVS5
#Il faudra intégrer à l'avenir la méthode pondérée pour calculer MCSVw
#référence biblio: Tamura 2020, https://doi.org/10.1007/s13246-020


            

  

#------- définition de la fonction de calcul du MCSV


def MCSv(iplan) : 
    
    MCSv_u_dict = {}
    MCSv_l_dict = {}
    

# boucle sur chaque faisceau du plan
    ib=1
    for k, b in zip(iplan['Beams'].keys(),iplan['Beams']) :

        #print k        

# inialisation des variables        
        max_leftleafs_u = []
        max_rightleafs_u = []
        max_aperture_u = []
        
        max_leftleafs_l = []
        max_rightleafs_l = []
        max_aperture_l = []
        
        MCSv_arc_u=0
        MCSv_arc_l=0

# Boucle sur chaque point de controle du faisceau pour récupérer les positions max de chacune des lames du MLC pour l'ensemble des CP.
# Données nécessaires pour calculer les termes AAV et LSV
               
        for cp_index, c in enumerate(iplan['Beams'][b]['ControlPoints']) :
            

            if cp_index==0: 
                max_leftleafs_u = list(c['pLeafPosLeft'])
                max_rightleafs_u = list(c['pLeafPosRight'])
                max_leftleafs_l = list(c['dLeafPosLeft'])
                max_rightleafs_l = list(c['dLeafPosRight'])
                continue

            # MLC Upper          
            for lx_index_u, lx_left_u in enumerate(c['pLeafPosLeft']):           
                if lx_left_u < max_leftleafs_u[lx_index_u]:  max_leftleafs_u[lx_index_u] = lx_left_u
                
            for lx_index_u, lx_right_u in enumerate(c['pLeafPosRight']):
                if lx_right_u > max_rightleafs_u[lx_index_u]:  max_rightleafs_u[lx_index_u] = lx_right_u

            # MLC Lower  
            for lx_index_l, lx_left_l in enumerate(c['dLeafPosLeft']):           
                if lx_left_l < max_leftleafs_l[lx_index_l]:  max_leftleafs_l[lx_index_l] = lx_left_l
                
            for lx_index_l, lx_right_l in enumerate(c['dLeafPosRight']):
                if lx_right_l > max_rightleafs_l[lx_index_l]:  max_rightleafs_l[lx_index_l] = lx_right_l    
                   
        # calcul du dénominateur de AAV - "max(pos)"
        max_aperture_u = [abs(x_u-y_u) for x_u,y_u in zip(max_rightleafs_u,max_leftleafs_u)]
        aav_denominator_u = sum(max_aperture_u)
        
        max_aperture_l = [abs(x_l-y_l) for x_l,y_l in zip(max_rightleafs_l,max_leftleafs_l)]
        aav_denominator_l = sum(max_aperture_l)
        

        # calcul du dénominateur de LSV - "posmax"
        pos_max_left_u = pos_max_right_u = abs(max(max_rightleafs_u)-min(max_leftleafs_u))
        pos_max_left_l = pos_max_right_l = abs(max(max_rightleafs_l)-min(max_leftleafs_l))

        

        
        lsv_cp_u = []
        aav_cp_u = []
        
        lsv_cp_l = []
        aav_cp_l = []
        
        for cp_index, c in enumerate(iplan['Beams'][b]['ControlPoints']) :  
                        
            leaf_pos_leftbank_u = [x_u for x_u,y_u in zip(c['pLeafPosLeft'], max_aperture_u) if y_u>0.]
            leaf_pos_rightbank_u = [x_u for x_u,y_u in zip(c['pLeafPosRight'], max_aperture_u) if y_u>0.]
            
            leaf_pos_leftbank_l = [x_l for x_l,y_l in zip(c['dLeafPosLeft'], max_aperture_l) if y_l>0.]
            leaf_pos_rightbank_l = [x_l for x_l,y_l in zip(c['dLeafPosRight'], max_aperture_l) if y_l>0.]
            
            
            # Calcul du AAV et LSV pour chaque CP
            
            lsv_left_u = 0
            lsv_right_u = 0
            aav_u = 0
            
            lsv_left_l = 0
            lsv_right_l = 0
            aav_l = 0
            
            il=0
            for lx_left_u, lx_right_u, lx_left_l, lx_right_l in zip(leaf_pos_leftbank_u, leaf_pos_rightbank_u,leaf_pos_leftbank_l, leaf_pos_rightbank_l):
                
                # Calcul du numérateur du AAV pour chaque CP

                # Collimateur Upper
                aav_u = aav_u + (lx_right_u-lx_left_u)
                # Collimateur Lower
                aav_l = aav_l + (lx_right_l-lx_left_l)

            
                # Calcul des numérateurs de LSV pour chaque CP et chaque banc droit et gauche

                if il==0: 
                    lx_left_n_u = lx_left_u
                    lx_right_n_u = lx_right_u
                    lx_left_n_l = lx_left_l
                    lx_right_n_l = lx_right_l
                    il=-1
                    continue
                
           
                
                lsv_left_u = lsv_left_u + (pos_max_left_u - abs(lx_left_n_u - lx_left_u))
                lsv_left_l = lsv_left_l + (pos_max_left_l - abs(lx_left_n_l - lx_left_l))

                lx_left_n_u = lx_left_u
                lx_left_n_l = lx_left_l
                
                lsv_right_u = lsv_right_u + (pos_max_right_u - abs(lx_right_n_u - lx_right_u))
                lx_right_n_u = lx_right_u
                lsv_right_l = lsv_right_l + (pos_max_right_l - abs(lx_right_n_l - lx_right_l))
                lx_right_n_l = lx_right_l
                

                                               
            # Calcul du LSV pour chaque CP (LSVij)
            lsv_left_u = lsv_left_u/((len(leaf_pos_leftbank_u)-1)*pos_max_left_u)
            lsv_right_u = lsv_right_u/((len(leaf_pos_rightbank_u)-1)*pos_max_right_u)
            lsv_left_l = lsv_left_l/((len(leaf_pos_leftbank_l)-1)*pos_max_left_l)
            lsv_right_l = lsv_right_l/((len(leaf_pos_rightbank_l)-1)*pos_max_right_l)


            # rangement des données dans un tableau pour chaque MLC            
            lsv_cp_u.append(lsv_left_u*lsv_right_u)
            lsv_cp_l.append(lsv_left_l*lsv_right_l)
                        

            
            # Calcul du AAV pour chaque CP (AAVij)
            aav_u = aav_u/aav_denominator_u
            aav_l = aav_l/aav_denominator_l
            
            # rangement des données dans un tableau pour chaque MLC  
            aav_cp_u.append(aav_u)
            aav_cp_l.append(aav_l)
            

# Boucle sur le CP pour calculer Le MCSV de chaque faisceau (MCSi)            
        for cp_index, c in enumerate(iplan['Beams'][b]['ControlPoints']) :
            
            if cp_index==0:
                aav_n_u = aav_cp_u[cp_index]
                lsv_n_u = lsv_cp_u[cp_index]
                aav_n_l = aav_cp_l[cp_index]
                lsv_n_l = lsv_cp_l[cp_index]
                w = c['Weight']

                continue
            


            MCSv_arc_u = MCSv_arc_u + ((aav_n_u+aav_cp_u[cp_index])/2*(lsv_n_u+lsv_cp_u[cp_index])/2*w)
            MCSv_arc_l = MCSv_arc_l + ((aav_n_l+aav_cp_l[cp_index])/2*(lsv_n_l+lsv_cp_l[cp_index])/2*w)
            w = c['Weight']

            
            aav_n_u = aav_cp_u[cp_index]
            lsv_n_u = lsv_cp_u[cp_index]
            aav_n_l = aav_cp_l[cp_index]
            lsv_n_l = lsv_cp_l[cp_index]

            
# Rangement des données MCS pour chaque faisceau dans un tableau (MCSi)
        MCSv_u_dict[k]=MCSv_arc_u
        MCSv_l_dict[k]=MCSv_arc_l
        ib=ib+1
        

    
    return MCSv_u_dict , MCSv_l_dict
 



#------------------- Plan Info -----------------------------------------------#

plan = {}
beams = {}
control_point = {}

CurrBS = get_current('BeamSet')
CurrPat = get_current('Patient')
CurrPlan = get_current('Plan')

#plan['PatientID'] = CurrPat.PatientID
plan['PlanName'] = CurrPlan.Name
plan['PrescriptionDose'] = CurrBS.Prescription.PrimaryPrescriptionDoseReference.DoseValue/100 #/100 to convert to Gy 
plan['Fractions'] = CurrBS.FractionationPattern.NumberOfFractions

if CurrBS.Beams.Count != 0 :
    for b in CurrBS.Beams:
        
        Name=b.Name
        beamMU = b.BeamMU
        beams[b.Name] = {'TotalMU':beamMU,'ControlPoints':[],'Nom_faisceau':Name}
             
        if  b.Segments.Count > 1  :

            for s in b.Segments :

                dr = s.DoseRate                     #MU/min
                rel_weight = s.RelativeWeight       # Poids relatif du CP (MUij/MUi)
                dose_segment = rel_weight*beamMU    #MU

                
                if dose_segment!=0 :
                    time_segment = dose_segment/dr*60   #seconds

                    gs = 2/time_segment                 #deg/sec  * control points are spaced 2 deg
                else :
                    gs = 0.
                    time_segment = 0.
                '''
                for mlcb in s.LeafPositions :

                    
                    for lpos in mlcb :
                        print lpos
                '''        
                leaf_pos_leftbank_u = list(s.LeafPositions[0])   # banc gauche du collimateur Upper
                leaf_pos_rightbank_u = list(s.LeafPositions[1])  # banc droit du collimateur Upper
                leaf_pos_leftbank_l = list(s.LeafPositions[2])   # banc gauche du collimateur Lower
                leaf_pos_rightbank_l = list(s.LeafPositions[3])  # banc droit du collimateur Lower
    

                                
                control_point={'Weight':rel_weight, 'DoseRate':dr, 'GantrySpeed':gs, 'Time': time_segment}     
                control_point['pLeafPosLeft'] = leaf_pos_leftbank_u
                control_point['pLeafPosRight'] = leaf_pos_rightbank_u
                control_point['dLeafPosLeft'] = leaf_pos_leftbank_l
                control_point['dLeafPosRight'] = leaf_pos_rightbank_l
                

# Rangement dans un tableau de l'ensemble des positions de lames et du poids relatif de chaque CP                                
                beams[b.Name]['ControlPoints'].append(dict(control_point))
                
            
    plan['Beams'] = beams



#--------------------- Call of funtions --------------------------------------#



mcsv_u,mcsv_l = MCSv(plan)   #récupération des valeurs de l'indice MCS pour chaque collimateur et chaque faisceau (MCSi)
planMU = 0
mcsv_mean_u = 0
mcsv_mean_l = 0

# Calcul du MCS pour l'ensemble du plan
for arc in plan['Beams'].keys():
    planMU = planMU + plan['Beams'][arc]['TotalMU']
    mcsv_mean_u = mcsv_mean_u + mcsv_u[arc]*plan['Beams'][arc]['TotalMU']
    mcsv_mean_l = mcsv_mean_l + mcsv_l[arc]*plan['Beams'][arc]['TotalMU']
if planMU != 0 : mcsv_mean_u = mcsv_mean_u/planMU
else : mcsv_mean_u = 'NaN'
if planMU != 0 : mcsv_mean_l = mcsv_mean_l/planMU
else : mcsv_mean_l = 'NaN'

list_MCSV_beam = []

# Affichage des informations détaillées pour chaque faisceau
for arc in plan['Beams'].keys():
    MCSV_beam = (mcsv_u[arc] + mcsv_l[arc]) / 2   # Calcul du MCS pour chaque faisceau suivant la méthode classique (MCS5)
    list_MCSV_beam.append(MCSV_beam)
    print (' ')
    print (' Indices de Modulation faisceau: ',plan['Beams'][arc]['Nom_faisceau'])
    print (' pMCSv : \t\t{0:.4f}'.format(mcsv_u[arc]))
    print (' dMCSv : \t\t{0:.4f}'.format(mcsv_l[arc]))
    print ('ID faisceau , MCSV_beam')
    print (plan['Beams'][arc]['Nom_faisceau'] , ' ', round(MCSV_beam,4))
    print (' -')


# Affichage des MCS5 pour chaque faisceau du plan
print (' ')
print ('-----------------------------------------------')
print ('ID faisceau , MCSV_beam')
i=0
for beam in CurrBS.Beams:
    MCSV_beam = list_MCSV_beam[i]
    print (beam.Name , ' ', round(MCSV_beam,4))
    i += 1
print ('-----------------------------------------------')
