MAX=$2
K_SN=$1

CLASS=$((${MAX}/10000))

if [ -e result.txt ]
then
rm result.txt
fi

if [ -e init.txt ]
then
rm init.txt
fi

seq 0 $((${CLASS}-1)) | parallel -j 5 'sh sn_simulate.sh {} ${K_SN}'

echo "\n Complete! \n"
