import numpy as np

def generate_jmolfile(jmol_filename, cube_filename, plane, png_filename):
    #WARNING : the plane is the average plane of the cycle
    # its equation is of the form : ax + by + c = z
    # To generate the equation of the plane parallel to this plane that goes through a specific point
    #  one must do this :
    #   ax + by -z + c = 0
    # V_normal = [a, b, -1]
    # equation of the plane that goes through the points we want:
    #  ax + by - z = d
    # Thus
    #  d = -1 * ( a * x_point + b * x_point - z_point )
    # As a point, we take the origin to which we add offset * the normalized normal vector
    origin = plane["origin"]
    offset = plane["offset"]
    vect = [-plane['a'], -plane['b'], 1]
    vect = vect / np.linalg.norm(vect)
    vect = offset * vect
    point = origin + vect
    d = -1 * ( plane['a'] * point[0] + plane['b'] * point[1] - point[2] )
    fio = open(jmol_filename, "w+")
    # set background
    fio.write('set backgroundColor "lightgray"\n'.format())
    # load cube file
    fio.write('load "{}"\n'.format(cube_filename))
    # settings
    fio.write('spacefill 5%\n'.format())
    fio.write('wireframe 5\n'.format())
    fio.write('color atoms black\n'.format())
    # define isosurface
    fio.write('isosurface iso_surface'.format())
    #   with bb
#    fio.write(
#        ' boundbox {{ {} {} {} }} {{ {} {} {} }}'.format(
#            plane["bbox"][0][0],
#            plane["bbox"][0][1],
#            plane["bbox"][0][2],
#            plane["bbox"][1][0],
#            plane["bbox"][1][1],
#            plane["bbox"][1][2]))
    fio.write(' color range -75 75 colorscheme "rwb"'.format())
    fio.write(
        ' "{0}" plane {{ {1} {2} {3} {4} }} map "{0}"\n'.format(
            cube_filename,
            plane['a'],
            plane['b'],
            -1,
            d))
    # define isocontours
    fio.write('isosurface iso_contours'.format())
    #   with bb
#    fio.write(
#        ' boundbox {{ {} {} {} }} {{ {} {} {} }}'.format(
#            plane["bbox"][0][0],
#            plane["bbox"][0][1],
#            plane["bbox"][0][2],
#            plane["bbox"][1][0],
#            plane["bbox"][1][1],
#            plane["bbox"][1][2]))
    fio.write(' color range -75 75 colorscheme "rwb"'.format())
    fio.write(' contour increment {{ -45,75,5 }}'.format())
    fio.write(
        ' "{0}" plane {{ {1} {2} {3} {4} }} map "{0}"\n'.format(
            cube_filename,
            plane['a'],
            plane['b'],
            -1,
            d))
    fio.write(' color isosurface phase red blue\n'.format())
    #
    fio.write('draw diameter 20 points [{{ {0[0]} {0[1]} {0[2]} }},{{ {1[0]} {1[1]} {1[2]} }}]\n'.format(
        [i for i in plane["bbox"][0]], [i for i in plane["bbox"][1]]))
    # write png file
    fio.write('write "{}" as PNGJ\n'.format(png_filename))
    fio.close()


def run_jmol(jmol_filename):
    proc = subprocess.Popen(["jmol",
                             jmol_filename],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            universal_newlines=True)
    stdout, stderr = proc.communicate()
    if (logger):
        # stderr is None because of stderr=subprocess.STDOUT in the Popen call
        logger.info(stdout)
