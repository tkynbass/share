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

CLASS=$((${MAX}/10000))

for K_MN in `seq 0 7`
do
for K_SN in `seq 0 $((7 - ${i}))`
do
K_SM=$((7 - ${K_MN} - ${K_SN}))

if [ ! -e ${K_MN}_${SN}_${K_SM}/ ]
then
mkdir ${K_SN}
else
if [ -e ${K_MN}_${SN}_${K_SM}/result.txt ]
then
rm result.txt
fi

if [ -e ${K_MN}_${SN}_${K_SM}/init.txt ]
then
rm init.txt
fi
fi

parallel -j 7 "sh exp_simulate.sh {} $((${K_MN} + 1)) $((${K_SN} + 1)) $((${K_SM} + 1))" ::: `seq ${START}  $((${CLASS}-1))`
done
done


echo "\n Complete! \n"
