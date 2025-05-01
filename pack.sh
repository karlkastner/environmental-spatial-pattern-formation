#!/bin/bash
d=$(basename $(pwd))
cd ../

zip  -9 ~/$d.zip -r $d \
	--exclude \
		$d'/.git/*' \
		$d'/.gitignore' \
		$d'/img/*' \
		$d'/mat/*' \
		$d'/other/*' \
		$d'/pack.sh' \
		$d'/patterns/metastudy*' \
		$d'/todo/*'

cp ../lib/auxiliar/addpath_recursive.m $d
zip  -j -9 ~/$d.zip $d/addpath_recursive.m
# ../../lib/auxiliar/addpath_recursive.m
#zip  -j -9 ~/$d.zip ../../lib/auxiliar/addpath_recursive.m

