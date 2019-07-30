#!/bin/sh

#  simulate_move.sh
#  
#
#  Created by tkym on 2019/07/30.
#  


##### define option (i:insertion) #####
#while getopts "i:" option
#do
#case ${option} in
#i)
#result_name="$OPTARG"
#;;
#\?)
#printf "\n     Usage : simulate_move.sh [-i result_name]   \n\n "  1>&2
#exit 1
#;;
#esac
#done

###### if there's not file to input #####
#if [ ! $result_name ]
#then
#printf "\n     Usage : simulate_move.sh [-i result_name]   \n\n "  1>&2
#exit 1
#fi
#
#if [ -e result/${result_name} ]
#then
#printf "\n      result/${result_name} already exists. Renew it? (\"yes\" or \"no\") : "
#read answer
#
#if [ $answer != yes ]
#then
#printf "\n\n    Please input not existence result_name.     \n\n"
#exit 0
#else
#rm -rf result/${result_name}
#fi
#fi

MAX=$1
for sample in `seq 1 ${MAX}`
do
printf "\t calculating ${sample}  "
if [ ! -e ${sample} ]
then mkdir ${sample}
fi
./octa_sampling ${sample} 100


done
