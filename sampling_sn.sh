#!/bin/sh

#  simulate_move.sh
#  
#
#  Created by tkym on 2019/07/30.
#  


MAX=$1

printf "\n"

for k in 0.5 0.6 0.7 0.8 0.9
do
printf "\t calculating ${k}  "
if [ ! -e ${sample} ]
then mkdir ${sample}
fi
./octa_sampling 100000 ${k}


done
