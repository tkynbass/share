#!/bin/sh

#  simulate_move.sh
#  
#
#  Created by tkym on 2018/06/06.
#  


##### define option (i:insertion) #####
while getopts "i:" option
do
case ${option} in
i)
result_name="$OPTARG"
;;
\?)
printf "\n     Usage : simulate_move.sh [-i result_name]   \n\n "  1>&2
exit 1
;;
esac
done

###### if there's not file to input #####
if [ ! $result_name ]
then
printf "\n     Usage : simulate_move.sh [-i result_name]   \n\n "  1>&2
exit 1
fi

if [ -e result/${result_name} ]
then
printf "\n      result/${result_name} already exists. Renew it? (\"yes\" or \"no\") : "
read answer

if [ $answer != yes ]
then
printf "\n\n    Please input not existence result_name.     \n\n"
exit 0
else
rm -rf result/${result_name}
fi
fi

###### result create #####
mkdir result/${result_name}
mkdir result/${result_name}/coordinates
cp src/common.h result/${result_name}/
cp src/func_move.h result/${result_name}/

echo ""

icc -O3 -xHost -DHAVE_SSE2=1 -DSFMT_MEXP=19937 -g -o simulate_move src/simulate_move.c src/dSFMT/dSFMT.c -qopenmp -fast

./src/simulate_move ${result_name}
