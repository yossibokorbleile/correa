import _correa

def create_polygon(path_file : str):
    return _correa.PyPolygon(path_file)

def print_polygon(p) :
    _correa.print_polygon(p)