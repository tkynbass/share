K_SN=$1
MAX=$2

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

parallel -j 5 "sh sn_simulate.sh {} ${K_SN}" ::: `seq ${START}  $((${CLASS}-1))`

echo "\n Complete! \n"
