class Atome:
    def __init__(self, label, x, y, z):
        self.valeur = {}
        self.setLabel(label)
        self.setCoord(x,y,z)

    def setLabel(self, label):
        self.valeur['label'] = label

    def setCoord(self, x,y,z):
        self.valeur['x'] = x
        self.valeur['y'] = y
        self.valeur['z'] = z

    def getlabel(self):
        return self.valeur['label']

    def getX(self):
        return self.valeur['x']

    def getY(self):
        return self.valeur['y']

    def getZ(self):
        return self.valeur['z']

    def getCoords(self):
        return [self.getX(), self.getY(), self.getZ()]

    def __str__(self):
        return "{} {} {} {}".format(self.valeur['label'], self.valeur['x'], self.valeur['y'], self.valeur['z'])

    def print(self):
        print(self.valeur)

# Si le code est appele  en tant que programme principal
# Sinon (cas d'un import), cettepartie du code n'est pas utilisée
if __name__=='__main__':
    atome = Atome('C',1,2,3)
    atome.print()
