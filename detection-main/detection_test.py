from os import remove
import numpy as np 
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import sort
import molecule as mol
import math as mt


def main():
    print("---------------------detection--------------------------")

    molecule = mol.Molecule("geoms/geom02.xyz")
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
    print()
    hd = molecule.dist_hydrogene()
    hd.sort()
    print(" la liste des distances entre hydrogènes : ")
    for i in range(len(hd)):
        print(hd[i])
    


    a = molecule.bay_dist_Hydrogene()
    c = molecule.get_coords_bay()
    #c1 = molecule.sort_bay(molecule.get_coords_bay())
    
    print("nombre d'éléments de la liste a : ", len(a))
    print()
    print(c, "nombre d'éléments : ", len(c))
    """for i in range(len(c)):
        print(c[i])"""
    for i in range(len(c)):
        for j in range(len(c[i])):
            plt.scatter(c[i][j][0],c[i][j][1])
    plt.show()


    


    
    print("------------------------ fin detection-test--------------------------------")


if __name__=="__main__":
    print("Appel en tant que programme principal")
    main()
