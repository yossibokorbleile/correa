import correa
p = correa.PyPolygon('/Users/yossi/AAU/ToMaCo/contours/contours/01kPa/24h_Trypsin24h_2-FITC_001_png_contour.csv')
p.vertices()
p.ellipse_min()
p.ellipse_max()
p.ellipse_lsq()
p.willmore()

