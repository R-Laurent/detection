import numpy as np
import json

"""
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

#
#from https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
#
a = np.array([[1, 2, 3], [4, 5, 6]])
print(a.shape)
json_dump = json.dumps({'a': a, 'aa': [2, (2, 3, 4), a], 'bb': [2]}, cls=NumpyEncoder)
print(json_dump)

to restore
json_load = json.loads(json_dump)
a_restored = np.asarray(json_load["a"])
print(a_restored)
print(a_restored.shape)
"""


class NumpyEncoder(json.JSONEncoder):
    # from
    # https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def dump_state(state, filename, indent=None):
    #
    # Serializes numpy lists in json
    #
    fio = open(filename, "w+")
    json.dump(state, fio, cls=NumpyEncoder, indent=indent)
    fio.close()

def load_state(filename):
    fio = open(filename, "r")
    state_loaded = json.load(fio)
    fio.close()
    #
    # Unserializes numpy lists in json
    #
    grids_loaded = state_loaded["grids"]
    grids = []
    for g in grids_loaded:
        grid = np.asarray(g)
        grids.append(grid)
    state = {}
    state["cycles"] = state_loaded["cycles"]
    state["planes"] = state_loaded["planes"]
    state["grids"] = grids
    return state

def dump_colormap(colormap, filename, indent=None):
    fio = open(filename, "w+")
    json.dump(colormap, fio, indent=indent)
    fio.close()

def load_colormap(filename):
    fio = open(filename, "r")
    colormap_loaded = json.load(fio)
    fio.close()
    return colormap_loaded
