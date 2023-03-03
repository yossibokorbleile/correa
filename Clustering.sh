#!/bin/sh

#  Clustering.sh
#  Created by Yossi Bokor on February 2, 2022.

matrix=$1
npts=$2
k=$3


while read file1; do
	file1s=""
	file1f=""
	file1e1=""
	file1e2=""
	file1e3=""
	file1wm=""
	file1wc=""
	file1wd=""

	while read file2; do
		echo $file1" and "$file2" starting comparison"
		bin/Comp2DShapes -i1 $file1 -i2 $file2 -d 4 &> Comp2DShapes.log

		sdist=$(grep "Sphericity" Comp2DShapes.log | awk '{print $6}')
		file1s+=" "$sdist","
		#echo $file1" "$file2" "$sdist >> $s

		fdist=$(grep "Frechet" Comp2DShapes.log | awk '{print $4}')
		file1f+=" "$fdist","
		#echo $file1" "$file2" "$fdist >> $f

		e1dist=$(grep "(inscribed ellipse" Comp2DShapes.log | awk '{print $5}')
		file1e1+=" "$e1dist","
		#echo $file1" "$file2" "$e1dist >> $e1

		e2dist=$(grep "(inscribing ellipse" Comp2DShapes.log | awk '{print $5}')
		file1e2+=" "$e2dist","
		#echo $file1" "$file2" "$e2dist >> $e2

		e3dist=$(grep "(least square ellipse" Comp2DShapes.log | awk '{print $6}')
		file1e3+=" "$e3dist","
		#echo $file1" "$file2" "$e3dist >> $e3

		wmdist=$(grep "(Willmore" Comp2DShapes.log | awk '{print $4}')
		file1wm+=" "$wmdist","
		#echo $file1" "$file2" "$wmdist >> $wm

		wcdist=$(grep "Wasserstein-curvature" Comp2DShapes.log | awk '{print $4}')
		file1wc+=" "$wcdist","
		#echo $file1" "$file2" "$wcdist >> $wc

		wddist=$(grep "2-Wasserstein" Comp2DShapes.log | awk '{print $6}')
		file1wd+=" "$wddist","
		#echo $file1" "$file2" "$wddist >> $wd
	done < $list

	echo ${file1s%?} >> $s
	echo ${file1f%?} >> $f
	echo ${file1e1%?} >> $e1
	echo ${file1e2%?} >> $e2
	echo $file1e3%?} >> $e3
	echo ${file1wm%?} >> $wm
	echo ${file1wc%?} >> $wc
	echo ${file1wd%?} >> $wd

done < $list
echo "Finished comparing them all"
