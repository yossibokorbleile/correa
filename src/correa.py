import _correa

def create_polygon(path_file : str):
    return _correa.PyPolygon(path_file)

def print_attributes(polygon : PyPolygon):
    