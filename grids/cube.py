import numpy as np
from data.constants import Bohr2Angstrom

def generate_cubefile(geom, grid_values, dx, dy, dz, nptx, npty, nptz):
    grid = None
    return generate_cubefile(geom, grid, grid_values, dx, dy, dz, nptx, npty, nptz)

def generate_cubefile(geom, grid, grid_values, dx, dy, dz, nptx, npty, nptz):
    """
    Generate a cube file for the given geometry and grid
    """
    cubefile = "ims.cube"
    fio = open(cubefile, "w+")
    nat = len(geom)
    fio.write("head 1\n".format())
    fio.write("head 2\n".format())
    fio.write("{0:5d} {1[0]:12.6f} {1[1]:12.6f} {1[2]:12.6f}\n".format(
        nat, grid[0] / Bohr2Angstrom))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        nptx, dx / Bohr2Angstrom, 0, 0))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        npty, 0, dy / Bohr2Angstrom, 0))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        nptz, 0, 0, dz / Bohr2Angstrom))
    for el in geom:
        if (el['label'].strip().lower() == "c"):
            nuclear_charge = 6
        elif (el['label'].strip().lower() == "h"):
            nuclear_charge = 1
        charge = 0.0
        fio.write(
            "{:5d} {:12.6f} {:12.6f} {:12.6f}  {:12.6f}\n".format(
                nuclear_charge,
                charge,
                el['x'] / Bohr2Angstrom,
                el['y'] / Bohr2Angstrom,
                el['z'] / Bohr2Angstrom))
    idx = 0
    for x in range(0, nptx + 1):
        for y in range(0, npty + 1):
            for z in range(0, nptz + 1):
                k = str(idx)
                val = 0
                if (k in grid_values.keys()):
                    val = grid_values[k]
                fio.write("{:12.6f}".format(val))
                if (z % 6 == 5):
                    fio.write("\n".format())
                idx = idx + 1
        fio.write("\n".format())
    return cubefile


def closest_node(node, nodes):
    # from
    # https://codereview.stackexchange.com/questions/28207/finding-the-closest-point-to-a-list-of-points
    return nodes[cdist([node], nodes).argmin()]


def generate_values_on_grid(geom, ims_grid, npts):
    """
    Put the ims_values on a cubic grid based on the proximity of a measured point and a grid point
    """
    xmin = min(a['x'] for a in ims_grid)
    xmax = max(a['x'] for a in ims_grid)
    ymin = min(a['y'] for a in ims_grid)
    ymax = max(a['y'] for a in ims_grid)
    zmin = min(a['z'] for a in ims_grid)
    zmax = max(a['z'] for a in ims_grid)
    #
    # Add 5% margin to handle volume plotting
    #
    xmin = xmin - xmin * .05
    ymin = ymin - ymin * .05
    zmin = zmin - zmin * .05
    xmax = xmax + xmax * .05
    ymax = ymax + ymax * .05
    zmax = zmax + zmax * .05
    logger.info(
        "xmin={} xmax={} ymin={} ymax={} zmin={} zmax={}".format(
            xmin,
            xmax,
            ymin,
            ymax,
            zmin,
            zmax))
    nptx = npty = nptz = npts
    tmp_lst = []
    # Generate grid : there is probably more efficient
    logger.info("Generating grid")
    for x in np.linspace(xmin, xmax, nptx):         # see dz
        for y in np.linspace(ymin, ymax, npty):     # see dz
            # dz = ( zmax - zmin ) / ( npts - 1 )
            for z in np.linspace(zmin, zmax, nptz):
                tmp_lst.append([x, y, z])
    grid = np.array(tmp_lst)
    grid_values = {}
    # remember dz = ( zmax - zmin ) / ( npts - 1 )
    dx = (xmax - xmin) / (nptx - 1)
    dy = (ymax - ymin) / (npty - 1)  # idem
    dz = (zmax - zmin) / (nptz - 1)  # idem
    for el in ims_grid:
        el_prox = closest_node([el['x'], el['y'], el['z']], grid)
        idz = (el_prox[2] - zmin) / dz
        idy = (el_prox[1] - ymin) / dy
        idx = (el_prox[0] - xmin) / dx
        index = idx * npty * nptz + idy * nptz + idz
        grid_values[str(int(index))] = el['ims']
        # store the value as a pair index (in str format), ims value: it is much better than storing a large number of
        # zero ims values on the grid
    return grid, grid_values, dx, dy, dz, nptx, npty, nptz

def generate_values_on_grid(geom, xmin, xmax, nptx, ymin, ymax, npty, zmin, zmax, nptz):
    """
    Put the ims_values on a cubic grid by computing the indices
    """
    #
    # Add 5% margin to handle volume plotting
    #
    xmin = xmin - xmin * .05
    ymin = ymin - ymin * .05
    zmin = zmin - zmin * .05
    xmax = xmax + xmax * .05
    ymax = ymax + ymax * .05
    zmax = zmax + zmax * .05
    logger.info(
        "xmin={} xmax={} ymin={} ymax={} zmin={} zmax={}".format(
            xmin,
            xmax,
            ymin,
            ymax,
            zmin,
            zmax))
    tmp_lst = []
    # Generate grid : there is probably more efficient
    logger.info("Generating grid")
    for x in np.linspace(xmin, xmax, nptx):         # see dz
        for y in np.linspace(ymin, ymax, npty):     # see dz
            # dz = ( zmax - zmin ) / ( npts - 1 )
            for z in np.linspace(zmin, zmax, nptz):
                tmp_lst.append([x, y, z])
    grid = np.array(tmp_lst)
    grid_values = {}
    # remember dz = ( zmax - zmin ) / ( npts - 1 )
    dx = (xmax - xmin) / (nptx - 1)
    dy = (ymax - ymin) / (npty - 1)  # idem
    dz = (zmax - zmin) / (nptz - 1)  # idem
    for el in ims_grid:
        el_prox = closest_node([el['x'], el['y'], el['z']], grid)
        idz = (el_prox[2] - zmin) / dz
        idy = (el_prox[1] - ymin) / dy
        idx = (el_prox[0] - xmin) / dx
        index = idx * npty * nptz + idy * nptz + idz
        grid_values[str(int(index))] = el['ims']
        # store the value as a pair index (in str format), ims value: it is much better than storing a large number of
        # zero ims values on the grid
    return grid, grid_values, dx, dy, dz, nptx, npty, nptz


def getIdx(value, val_min, val_max=None, npts=None, dval=None):
    """
    Return the index of a 1D value on a linear grid.
    Either val_max and npts
    or
    dval
    must be present
    """
    if (dval==None):
        dval = (val_max-val_min)/(npts - 1)
    return np.round((value-val_min)/dval)

def getIdx_cube(x, y, z, xmin, xmax, nptx, ymin, ymax, npty, zmin, zmax, nptz):
    """
    Returns the index of a point on a 3D rectangular grid
    The cube order is assumed (from http://paulbourke.net/dataformats/cube/)
       for (ix=0;ix<NX;ix++) {
           for (iy=0;iy<NY;iy++) {
               for (iz=0;iz<NZ;iz++) {
                   printf("%g ",data[ix][iy][iz]);
                   if (iz % 6 == 5)
                       printf("\n");
                   }
               printf("\n");
           }
       }
    """
    idx = getIdx(x, xmin, xmax, nptx)
    idy = getIdx(y, ymin, ymax, npty)
    idz = getIdx(z, zmin, zmax, nptz)
    index = idx * npty * nptz + idy * nptz + idz
    return index

def generate_cube_volume(geom, points, xmin, xmax, nptx, ymin, ymax, npty, zmin, zmax, nptz):
    index_values = []
    for pt in points:
        x = pt[0]
        y = pt[1]
        z = pt[2]
        index_values.append(getIdx_cube(x, y, z, xmin, xmax, nptx, ymin, ymax, npty, zmin, zmax, nptz))
    return index_values

def generate_cubefile_new(filename, geom, grid, xmin, xmax, nptx, ymin, ymax, npty, zmin, zmax, nptz):
    """
    Generate a cube file for the given geometry and grid
    """
    dx = ( xmax - xmin ) / ( nptx - 1 )
    dy = ( ymax - ymin ) / ( npty - 1 )
    dz = ( zmax - zmin ) / ( nptz - 1 )
    cubefile = filename
    nat = len(geom.atoms)
    fio = open(cubefile, "w+")
    nat = len(geom.atoms)
    fio.write("head 1\n".format())
    fio.write("head 2\n".format())
    fio.write("{0:5d} {1[0]:12.6f} {1[1]:12.6f} {1[2]:12.6f}\n".format(
        nat, [xmin/Bohr2Angstrom, ymin/Bohr2Angstrom, zmin/Bohr2Angstrom]))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        nptx, dx / Bohr2Angstrom, 0, 0))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        npty, 0, dy / Bohr2Angstrom, 0))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        nptz, 0, 0, dz / Bohr2Angstrom))
    for el in geom.atoms:
        if (el['label'].strip().lower() == "c"):
            nuclear_charge = 6
        elif (el['label'].strip().lower() == "h"):
            nuclear_charge = 1
        charge = 0.0
        fio.write(
            "{:5d} {:12.6f} {:12.6f} {:12.6f}  {:12.6f}\n".format(
                nuclear_charge,
                charge,
                el['x'] / Bohr2Angstrom,
                el['y'] / Bohr2Angstrom,
                el['z'] / Bohr2Angstrom))
    idx = 0
    for ix in range(0, nptx):
        for iy in range(0, npty):
            for iz in range(0, nptz):
                val=0
                if idx in grid:
                    val=1
                fio.write("{:12.6f}".format(val))
                if (iz % 6 == 5):
                    fio.write("\n".format())
                idx = idx + 1
        fio.write("\n".format())
    return cubefile

