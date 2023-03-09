import sys
sys.path.append("/Users/yossi/correa_private/src/")

import correa
p = correa.create_polygon('/Users/yossi/AAU/ToMaCo/contours/contours/01kPa/24h_Trypsin24h_2-FITC_001_png_contour.csv')
correa.print_polygon(p)
