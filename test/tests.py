import correa

p1 = correa.create_polygon("../examples/contour1.csv")
p2 = correa.create_polygon("../examples/contour2.csv")

comp = correa.compare_polygons(p1, p2, q=100, verbose=True)
print(comp)