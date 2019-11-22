#!/bin/sh


NUM=$1
K_MN=$2
K_SN=$3
K_SM=$4

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
./octa_sampling_sn ${sample} 300 ${K_MN} ${SN} ${K_SM}
printf "\r"
done
