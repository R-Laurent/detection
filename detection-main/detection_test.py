from os import remove
import numpy as np 
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import sort
import molecule as mol
import math as mt


def main():
    print("---------------------detection--------------------------")

    molecule = mol.Molecule("geoms/geom01.xyz")
    print("------------------detection_test-------------------------------")
    
    z = molecule.preBarycenterZig_Zag(molecule.zigzag)
    k = molecule.preBarycenterZig_Zag(molecule.kRegion)
    
    molecule.printVoisins_L(z)

    print()
    print(z)

    print()
    print("barycentres de zig zag :")
    print(molecule.getBarycenter(molecule.preBarycenterZig_Zag(molecule.zigzag)))


    print()
    h = molecule.hydrogene_list()
    print("voisins hydrogènes : ",h , "la longueur de la liste est : ", len(h))
    print()
    hd = molecule.dist_hydrogene()
    hd.sort()
    print(" la liste des distances entre hydrogènes : ", hd)

    """
    print()
    print("Bay REGION")
    a = molecule.bay_region
    print(a, "nombre de bay : ", len(a))
    print()
    print("premier élément de bay region")
    print(a[0])
    m = a[0]
    print()
    print("voici le premier élément plus en profondeur : ")
    print(a[0][0])


    
    a2 = molecule.bay_region_label()
    a2_1 = a2
    print("bay region label: ")
    print()
    
    for i in range(len(a2)):
        print(a2[i])


    print(len(a2))







    for i in range(3):
        del a[i]
    
    print("voici la liste bay après del ??")
    print()
    print(a, len(a))


    for i in range(len(a)): 
        for j in range(len(a[i])):
            plt.scatter(a[i][j][0],a[i][j][1], c = 'b' ) 
    plt.title("DETECTION_TEST")
    plt.show()
    


    for i in range(len(m)):
        plt.scatter(m[i][0], m[i][1], c = 'b')
    plt.title("DETECTION_TEST")
    plt.show()

    print()
    print()
    print()

    h = molecule.sort_bay(molecule.bay_has_same_direction())
    print("la bay region après la comparaison des vecteurs : ")
    print(h, "nombre d'éléments de la liste : ", len(h))

    for i in range(len(h)): 
        for j in range(len(h[i])):
            plt.scatter(h[i][j][0],h[i][j][1], c = 'b' ) 
    plt.title("DETECTION_TEST_geom01")
    plt.show()
    """
    print("------------------------ fin detection-test--------------------------------")


if __name__=="__main__":
    print("Appel en tant que programme principal")
    main()
