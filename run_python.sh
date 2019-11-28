#!/bin/bash
#
#PBS -q condo
#PBS -l nodes=1:ppn=8
#PBS -N run_python
#PBS -j oe
#PBS -V
#PBS -o run_python.out
#PBS -m a
#PBS -M mrdavidson@ucsd.edu

datebeg=$(date +%s)

cat $PBS_NODEFILE
echo 'The list above shows the nodes this job has exclusive access to.'
cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh

## BEGIN Script

python3 $1

## END Script

dateend=$(date +%s)
diffdate=$(($dateend-$datebeg))

if [ "$diffdate" -lt "120" ]
	then echo "Only $diffdate seconds have elapsed!"
	exit 1
fi

exit 0


