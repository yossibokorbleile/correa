#!/bin/sh

#  CompareAll.sh
#  Created by Yossi Bokor on February 2, 2022.

list=$1
name=$2
s="Sphericity-$name"
f="Frechet-$name"
e1="Inscribed-ellipse-$name"
e2="Inscribing-ellipse-$name"
e3="Least-square-ellipse-$name"
wm="Willmore-$name"
wc="Wasserstein-curvature-$name"
wd="Wasserstein-distance-$name"
rm $s $f $e1 $e2 $e3 $wm $wc $wd
touch $s $f $e1 $e2 $e3 $wm $wc $wd

while read f1; do
	file1="Contours/contours_24h/01kPa/contour_$f1.csv"
	file1s=""
	file1f=""
	file1e1=""
	file1e2=""
	file1e3=""
	file1wm=""
	file1wc=""
	file1wd=""

	while read f2; do
		file2="Contours/contours_24h/01kPa/contour_$f2.csv"
		echo $file1" and "$file2" starting comparison"
		bin/Comp2DShapes -i1 $file1 -i2 $file2 -d 4 &> CellComp.log

		sdist=$(grep "Sphericity" CellComp.log | awk '{print $6}')
		file1s+=" "$sdist","
		#echo $file1" "$file2" "$sdist >> $s

		fdist=$(grep "Frechet" CellComp.log | awk '{print $4}')
		file1f+=" "$fdist","
		#echo $file1" "$file2" "$fdist >> $f

		e1dist=$(grep "(inscribed ellipse" CellComp.log | awk '{print $5}')
		file1e1+=" "$e1dist","
		#echo $file1" "$file2" "$e1dist >> $e1

		e2dist=$(grep "(inscribing ellipse" CellComp.log | awk '{print $5}')
		file1e2+=" "$e2dist","
		#echo $file1" "$file2" "$e2dist >> $e2

		e3dist=$(grep "(least square ellipse" CellComp.log | awk '{print $6}')
		file1e3+=" "$e3dist","
		#echo $file1" "$file2" "$e3dist >> $e3

		wmdist=$(grep "(Willmore" CellComp.log | awk '{print $4}')
		file1wm+=" "$wmdist","
		#echo $file1" "$file2" "$wmdist >> $wm

		wcdist=$(grep "Wasserstein-curvature" CellComp.log | awk '{print $4}')
		file1wc+=" "$wcdist","
		#echo $file1" "$file2" "$wcdist >> $wc

		wddist=$(grep "2-Wasserstein" CellComp.log | awk '{print $6}')
		file1wd+=" "$wddist","
		#echo $file1" "$file2" "$wddist >> $wd
	done < $list

	echo ${file1s} >> $s
	echo ${file1f} >> $f
	echo ${file1e1} >> $e1
	echo ${file1e2} >> $e2
	echo ${file1e3} >> $e3
	echo ${file1wm} >> $wm
	echo ${file1wc} >> $wc
	echo ${file1wd} >> $wd

done < $list
echo "Finished comparing them all"
