K_SN=$1
MAX=$2

CLASS=$((${MAX}/10000))

if [ -e result.txt ]
then
rm result.txt
fi

if [ -e init.txt ]
then
rm init.txt
fi

if [ ! -e ${K_SN}/ ]
then
mkdir ${K_SN}
fi

parallel -j 5 "sh sn_simulate.sh {} ${K_SN}" ::: `seq 0 $((${CLASS}-1))`

echo "\n Complete! \n"
