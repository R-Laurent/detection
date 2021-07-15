import molecule as mol

def barycenters(filename):
    molecule = mol.Molecule(filename)
    v = {
        'K_region' : molecule.getBarycenter(molecule.preBarycenterK_Regions(molecule.kRegion())),

        'Zig_Zag_Region' : molecule.getBarycenter(molecule.preBarycenterZig_Zag(molecule.zig_zag())),

        'Bay_Region_And_Gulf_Region' : molecule.getBarycenter(molecule.bay_dist_Hydrogene())

    }
    return v

#print(barycenters("geoms/acene5.xyz"))
