#!/bin/bash
set -e
echo $#
if [[ $# -ne 2 ]];  then
    echo "wrong number of arguments supplied"
	exit 1
fi

mkdir -v ~/Documents/Work/results/$1
cp -v backup/*$2.bkp ~/Documents/Work/results/$1
mkdir ~/Documents/Work/results/$1/scripts
cp -v scripts/calculationinformation.info ~/Documents/Work/results/$1/scripts
sh ~/Documents/Work/results/rename.sh ~/Documents/Work/results/$1
