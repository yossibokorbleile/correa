import correa

p1 = correa.create_polygon("/Users/yossi/correa/examples/contour1.csv")
p2 = correa.create_polygon("/Users/yossi/correa/examples/contour2.csv")
correa.compare_polygons(p1,p2)
correa.print_polygon(p1)
correa.print_polygon(p2)
