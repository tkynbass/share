MAX=$1

##### define option (i:insertion) #####
#while getopts "s:" option
#do
#case ${option} in
#-s)
#START="$OPTARG"
#;;
#\?)
#START=0
#;;
#esac
#done


#if [ ${START} -ne 0 ]
#then
#START=$((${START}/10000))
#fi

START=0

CLASS=$((${MAX}/1000))

for K_MN in `seq 0 7`
do
for K_SN in `seq 0 $((7 - ${K_MN}))`
do
K_SM=$((7 - ${K_MN} - ${K_SN}))

if [ ! -e $((${K_MN} + 1))_$((${K_SN} + 1))_$((${K_SM} + 1))/ ]
then
mkdir $((${K_MN} + 1))_$((${K_SN} + 1))_$((${K_SM} + 1))
else
if [ -e $((${K_MN} + 1))_$((${K_SN} + 1))_$((${K_SM} + 1))/result.txt ]
then
rm $((${K_MN} + 1))_$((${K_SN} + 1))_$((${K_SM} + 1))/result.txt
fi

if [ -e $((${K_MN} + 1))_$((${K_SN} + 1))_$((${K_SM} + 1))/init.txt ]
then
rm $((${K_MN} + 1))_$((${K_SN} + 1))_$((${K_SM} + 1))/init.txt
fi
fi
echo "\t Calculating... $((${K_MN} + 1)) $((${K_SN} + 1)) $((${K_SM} + 1))"
parallel --no-notice -j 7 "sh linear_simulate.sh {} $((${K_MN} + 1)) $((${K_SN} + 1)) $((${K_SM} + 1))" ::: `seq ${START}  $((${CLASS}-1))`
done
done


echo "\n Complete! \n"
