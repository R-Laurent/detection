import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import atome 
import pandas as pd
import math as mt


def read_xyz(filename):
    # RL
    # à la base je voulais faire des sous-programmes qui s'appelaient entre eux mais je suis assez rouillé dans cette notion de python, du coup j'ai fait tout dans un
    # je vais me renseigner plus sur cela
    # pour l'instant j'ai fait la liste des carbones voisins, je préfère vous demander si c'est bon avant de passez aux hydrogènes
    # RL
    atomlist = []
    xyz = open(filename)
    n_atoms = int(xyz.readline())
    title = xyz.readline()
    for line in xyz:
        atom, x, y, z = line.split()
        atomlist.append(atome.Atome(atom, float(x), float(y), float(z)))
    xyz.close()

    if n_atoms != len(atomlist):
        raise ValueError("File says %d atoms but read %d points." %
                         (n_atoms, len(atomlist)))
    return atomlist


class Molecule:
    def __init__(self, filename):
        self.neighbour = []
        self.atomlist = read_xyz(filename)
        self.voisin()
        self.printVoisins()
        self.zigzag = self.zig_zag()
        #self.L_barycenter = self.preBarycenter()
        #self.barycenters = self.getBarycenter()
        self.kRegion = self.kRegion()
        self.bay_region = self.bay_region_coords()
        
        

    """
    Pour chaque atom de carbone, on fait une liste de voisins
    QUI DEBUTE PAR L'ATOME CONSIDERE
    Donc un atome a 3 voisins aura une liste de taille 4:
    lui-même puis les 3 voisins quel que soit leur type
    """

    def voisin(self):
        """
        Génère la liste des voisins pour chaque atome de carbone
        """
        for at1 in self.atomlist:
            if at1.getlabel() == 'C':
                pos1 = np.array(at1.getCoords())
                voisins_1 = []
                voisins_1.append(at1)
                for at2 in self.atomlist:
                    pos2 = np.array(at2.getCoords())
                    V12 = pos2 - pos1
                    dist = np.linalg.norm(V12)
                    if (dist > 1e-6) and (dist < 2):
                        voisins_1.append(at2)
                self.neighbour.append(voisins_1)

    def printVoisins(self):
        for el in self.neighbour:
            print("Voisins de {} :".format(el[0]))
            for at in el:
                at.print()
    
    def printVoisins_L(self, L):
        for el in self.neighbour:
            print("Voisins de {} :".format(el[0]))
            for at in el:
                at.print()



    """
    Renvoit la liste des carbones ayant 3 voisins carbone
    """

    def zig_zag(self):
        v = []
        neighbour2 = []
        neighbour2 = self.neighbour
        print()
        for el in self.neighbour:
            store = True
            for i in range(4):
                if el[i].getlabel() == 'H':                                                                
                    store = False
            if store == True: 
                v.append(el)

        return v
    
    def zig_zag_boolean(self, atom):
        z = True
        for el in self.neighbour:
            if el[0] == atom:
                for i in range(4):
                    if el[i].getLabel() == 'H':
                        z = False
        return z
    
    def kRegion(self):
        v = []
        neighbour2 = []
        neighbour2 = self.neighbour
        print()
        for el in self.neighbour:
            store = True
            for i in range(4):
                if el[i].getlabel() == 'H':                                                                
                    store = False
                    v.append(el)
                    
        return v
    
    
    def kRegion_boolean(self, atom):
        v = []
        neighbour2 = []
        neighbour2 = self.neighbour
        print()
        for el in self.neighbour:
            store = True
            if el[0] == atom:
                for i in range(4):
                    if el[i].getlabel() == 'H':                                                                
                        store = False
          
        return store
    
    def H_in_Kregion(self, a):
        L = self.kRegion
        b = False
        for el in L:
            if a in el:
                b = True
        return b 
                
    def bay_has_boolean_H_neighbor(self, a):
        L = self.preBarycenterZig_Zag(self.zigzag)
        for el in L:
            b = True
            for i in range(len(el)):
                if self.hasHneighbor(el[i]) == True and self.getHneighbor(el[i]) == a:
                    b = False
                    return b
        return b

    def bay_has_boolean_H_neighbor(self, a):
        M = self.preBarycenterZig_Zag(self.zigzag)
        b = False
        for i in range(len(M)):
            for el in range(len(M[i])):
                if M[i][el] == a:
                    b = True
                    return b
        return b


    def bay_region_no_coords(self):
        v = []
        zig_zag = self.zig_zag()
        for el in zig_zag:
            for el2 in zig_zag:
                if self.are_neighbour_bay(el[0],el2) and el2 != el:
                    b = []
                    print("ils sont voisins")
                    ok1, a1, h1 = self.bay_C_H(el[0])
                    ok2, a2, h2 = self.bay_C_H(el2[0])
                    if ok1 and ok2:
                        print(self.bay_has_boolean_H_neighbor(h1))
                        print(self.bay_has_boolean_H_neighbor(h2))
                        print("ils sont vrais")
                        b.append(el[0])
                        b.append(a1)
                        b.append(h1)
                        b.append(el2[0])
                        b.append(a2)
                        b.append(h2)
                        v.append(b)
        
        return v
    
    
    def bay_region_label(self):
        v = []
        zig_zag = self.zig_zag()
        for el in zig_zag:
            for el2 in zig_zag:
                if self.are_neighbour_bay(el[0],el2) and el2 != el:
                    b = []
                    print("ils sont voisins")
                    ok1, a1, h1 = self.bay_C_H(el[0])
                    ok2, a2, h2 = self.bay_C_H(el2[0])
                    if ok1 and ok2:
                        print(self.bay_has_boolean_H_neighbor(h1))
                        print(self.bay_has_boolean_H_neighbor(h2))
                        print("ils sont vrais")
                        b.append(h1.getlabel())
                        b.append(a1.getlabel())
                        b.append(el[0].getlabel())
                        b.append(el2[0].getlabel())
                        b.append(a2.getlabel())
                        b.append(h2.getlabel())
                        v.append(b)
        
        return v

    def bay_region_coords(self):
        v = []
        zig_zag = self.zig_zag()
        for el in zig_zag:
            for el2 in zig_zag:
                if self.are_neighbour_bay(el[0],el2) and el2 != el:
                    b = []
                    print("ils sont voisins")
                    ok1, a1, h1 = self.bay_C_H(el[0])
                    ok2, a2, h2 = self.bay_C_H(el2[0])
                    if ok1 and ok2:
                        print(self.bay_has_boolean_H_neighbor(h1))
                        print(self.bay_has_boolean_H_neighbor(h2))
                        print("ils sont vrais")
                        b.append(h1.getCoords())
                        b.append(a1.getCoords())
                        b.append(el[0].getCoords())
                        b.append(el2[0].getCoords())
                        b.append(a2.getCoords())
                        b.append(h2.getCoords())
                        v.append(b)
        
        return v
    
    def bay_has_same_direction(self):
        L = self.bay_region
        v = []
        b = True
        for i in range(len(L)):
            print("marqueur : ", i)
            pos1 = np.array(L[i][0]) 
            pos2 = np.array(L[i][1])
            vect1 = pos1 - pos2
            
            pos3 = np.array(L[i][1]) 
            pos4 = np.array(L[i][2])
            vect2 = pos3 - pos4
            n1 = np.cross(vect1, vect2)

            pos5 = np.array(L[i][3]) 
            pos6 = np.array(L[i][4])
            vect3 = pos5 - pos6
            
            pos7 = np.array(L[i][4]) 
            pos8 = np.array(L[i][5])
            vect4 = pos7 - pos8
            n2 = np.cross(vect3, vect4)

            z = np.vdot(n1,n2)
            if z > 0:
                b = True
                v.append(L[i])
        return v

    def sort_bay(self, L):
        for i in range(len(L)):
            for j in range(len(L)):
                if L[i] == L[j] and i != j:
                    del L[j]
        return L

    def are_neighbour_bay(self,a,b):
        k = False
        if a in b:
            k = True            # dit si deux atome sont voisins 
        return k
    
    def bay_C_H(self, a):
        Z = self.preBarycenterZig_Zag(self.zigzag)
        b = False
        for el in self.neighbour:
            if el[0] == a:
                for i in range(len(el)):
                    if self.hasHneighbor(el[i]) == True:    
                        b = True 
                        return b, el[i], self.getHneighbor(el[i])
        return b, None, None                                                                                                                                        # en gros le pb c'est que deux atomes de zig zag se compare, 
                                                                                                                                        # il y en a un qui est vraiment dans un zig zag mais l'autre non
                                                                                                                                        # peut être qu'il faut mettre dans une des conditions que 
                                                                             
                                                                                                                                                # l'atome d'hydrogène ne dois pas doit appartenir à une K region et a une zig zag

    

    def hasHneighbor(self, atom):
        H=False
        Hy = []
        for el in self.neighbour:
            if el[0] == atom: 
                for i in range(len(el)):
                    if el[i].getlabel() == 'H':
                        H = True
                        Hy.append(el[i])
        return H
    
    def getHneighbor(self, atom):
        H=False
        a = atom
        for el in self.neighbour:
            if el[0] == atom: 
                for i in range(len(el)):
                    if el[i].getlabel() == 'H':
                        H = True
                        a = el[i]
        return a

    def preBarycenterZig_Zag(self, M):
        Lbarycentre= []
        for el in M:
            L=[]
            L.append(el[0])
            k = 0
            for i in range(len(el)):
                if self.hasHneighbor(el[i]) == True:      # peut être si il est vrai deux fois 
                    k = k + 1                                      # alors on ajoute 
                    a = el[i]
                    b = self.getHneighbor(el[i])
                    L.append(el[i])
                    L.append(self.getHneighbor(el[i]))
            if k == 2:
                Lbarycentre.append(L)
        return Lbarycentre
    
    def preBarycenterK_Regions(self, M):
            Lbarycentre= []
            for el in M: 
                for i in range(1,3):
                    if self.hasHneighbor(el[i]):
                        K = []
                        K.append(el[0])
                        K.append(el[i])
                        K.append(self.getHneighbor(el[i]))
                        K.append(self.getHneighbor(el[0]))
                        Lbarycentre.append(K)
            
            return Lbarycentre
        
    
    
    def getBarycenter(self, L):
        barycenters = []
        for el in L: 
            nat = len(el)
            barycenter = np.array([0,0,0])
            for i in range(len(el)):
                barycenter = np.array(el[i].getCoords())/nat + barycenter
            barycenters.append(barycenter)
        return barycenters
    


    def scatter1(self, b):
            """
            Affiche les points contenus dans la liste passee en argument
            """
            x = []
            y = []
            mark = [".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "P", "*",
                "h", "H", "+", "x", "X", "D", "d", "|", "_", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
            k=0
            for el in self.neighbour:  
                k=k+1  
                for i in range(len(el)):
                    if el[i].getlabel() == 'H':
                        coords = el[i].getCoords()
                        plt.scatter(coords[0], coords[1],c = 'r' )
                    else:
                        coords = el[i].getCoords()
                        plt.scatter(coords[0], coords[1],c = 'b' )
                    

                
            plt.scatter(x,y)
            for i in range(len(b)):
                plt.scatter(b[i][0],b[i][1], c = 'black', marker='*')
            plt.title("DETECTION")
            plt.show()

   








                                                                                                                
                                                                                    
    

        





       
                                              
    

