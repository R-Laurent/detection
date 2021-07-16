import numpy as np 
import matplotlib.pyplot as plt
import molecule as mol



# 1/ J'ai refactorisé le programme en faisant un programme principal (ci-dessous) et en faisant les modifications suivantes
# 2/ Pour les variables, de mon poitn de vue, il vaut mieux garder uniquement les variables les plus importantes
# comme globales le reste, on les definit à la volée
# 3/ Je propose de créer une classe molecule (une sorte d'équivalent d'un objet) pour conserver tout ce qui concerne une
# molecule: lire un fichier xyz, écrire un fichier xyz, trouiver les voisins, etc.
# 4/ J'ai mis cette classe dans un fichier à part
def main():
    print("---------------------detection--------------------------")

    molecule = mol.Molecule("geoms/geom02.xyz")
    z = molecule.preBarycenterZig_Zag(molecule.zig_zag())
    k = molecule.preBarycenterZig_Zag(molecule.kRegion())
    print()
    print("ZIG-ZAG : ")
    print()
    for el in molecule.zig_zag():
        el[0].print()
    print()
    print("K REGIONS : ")
    for el in molecule.kRegion():
        el[0].print()

    print()
    molecule.preBarycenterK_Regions(molecule.kRegion())
    molecule.preBarycenterZig_Zag(molecule.zig_zag())
    print("BARYCENTERS : ")
    print()
    print("Barycenters Zig-Zag regions : ")
    print(molecule.getBarycenter(molecule.preBarycenterZig_Zag(molecule.zig_zag())))
    a = molecule.getBarycenter(molecule.preBarycenterZig_Zag(molecule.zig_zag()))
    print()
    print("nombre de zig-zag : ", len(a))
    print()
    print("Barycenters K Regions : ")
    print(molecule.getBarycenter(molecule.preBarycenterK_Regions(molecule.kRegion())))
    b = molecule.sort_barycenters(molecule.getBarycenter(molecule.preBarycenterK_Regions(molecule.kRegion())))
    print()
    print("nombre de K : ", len(b))
    print()
    print("premier élément du premier élément de la K : ", b[0][1])
    print()

    print("logueur liste bay_dist : ", len(molecule.bay_dist_Hydrogene()))
    c1 = molecule.getBarycenter(molecule.bay_dist_Hydrogene())
    print()
    print(c1, "la longueur de la liste est : ", len(c1))
    c = molecule.sort_barycenters(c1)
    print("Barycenters bay regions Regions : ")
    print(c, "nombre de barycentres : ", len(c))
    print()
    #print("le premier élément : ", c[0][0])

    Barycenters = []

    for i in range(len(b)):
        Barycenters.append(b[i])
    for i in range(len(a)):
        Barycenters.append(a[i])
    for i in range(len(c)):
        Barycenters.append(c[i])
    
    print()
    molecule.printVoisins_L(k)

    print(molecule.scatter1(Barycenters))
    
   
    print("---------------------fin detection--------------------------")

if __name__=="__main__":
    print("Appel en tant que programme principal")
    main()

#YC

# Je ne comprends pas ce que ça fait:
# Je comprends que les boucles sur i et j cherchent qui est voisin de qui mais les autres, je ne sais pas.
# À mon avis il faut d'abord faire la liste des voisins, puis la traiter au lieu de vouloir tout faire en même temps.
# Mon algo serait


# faire la liste des premiers voisins                           # RL : j'ai essayé de faire cela pour  l'instant 
# faire une boucle sur tous les atomes:
#   si un atome a 2 voisins qui portent chacun un hydrogène, alors c'est un site que je recherche
#   donc je calcule le barycentre

