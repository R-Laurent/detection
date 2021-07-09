#          for i in range(len(L)):
#  #            for j in range(i+1, len(L)):
#  # ici il y a un problème.
#  #prenons le cas où 2 est voision de 1 et 3 cet algorithme teste le voisinage de 2 avec 3 mais pas de 3 avec 2 car 3>2
#  # on peut donc faire:
#              for j in range(len(L)):
#                  if (i!=j) :
#                      v1=[self.coordinates[i][0]-self.coordinates[j][0],self.coordinates[i][1],self.coordinates[j][1]]
#                      d1=np.linalg.norm(v1)
#                      if d1<=2 and d1!=0:
#                          M=[self.coordinates[i],self.coordinates[j]]
#                          self.neighbor.append(M)
#          for i in range(len(self.neighbor)):
#              print(self.neighbor[i])