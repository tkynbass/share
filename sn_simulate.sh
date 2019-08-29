#!/bin/sh



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

NUM=$1
K_SN=$2

start=$((${NUM}*10000+1))
end=$((${start}+9999))

#if [ -e result.txt ]
#then
#rm result.txt
#fi

#if [ -e init.txt ]
#then
#rm init.txt
#fi

printf "\n"
for sample in `seq ${start} ${end}`
do
printf "\t calculating ${sample}  "
if [ ! -e ${sample} ]
then mkdir ${sample}
fi
./octa_sampling_sn ${sample} 100 ${K_SN}
printf "\r"
done
