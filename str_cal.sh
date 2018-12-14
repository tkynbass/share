#!bin/sh
NUMBER=$1
ALGORHITHM=$2

for i in 1long 1short 2long 2short 3long 3short;
do
echo "\n\t${i}"

mkdir -p opt_${i}
mv ${ALGORHITHM}_data_${i}.txt opt_${i}

cd opt_${i}
cp ../strcal .
./strcal 2000 ${ALGORHITHM}_data_${i}.txt
cd ..
done
