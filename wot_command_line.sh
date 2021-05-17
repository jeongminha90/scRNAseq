#!/bin/bash
set -e
set -u
set -o pipefail


wot optimal_transport --matrix ../DIC_fD5_iPSC_mat_all.txt --cell_days ../days_DIC_fD5_iPSC_all.txt



sample_list=('7' '6' '5')
num_sample=${#sample_list[@]}

mkdir iPSC
cd iPSC
for (( i=0;i<$num_sample;i++ ));
do
mkdir D${sample_list[${i}]}
cd D${sample_list[${i}]}
wot trajectory --tmap ../tmaps --cell_set ../../cell_cluster_DIC_fD5_iPSC_all.gmt --day ${sample_list[${i}]}
mv wot_trajectory.txt D${sample_list[${i}]}_wot_trajectory.txt
cd ../

done


mkdir ../iNSC
cd ../iNSC

wot optimal_transport --matrix ../DIC_fD5_iNSC_mat_all.txt --cell_days ../days_DIC_fD5_iNSC_all.txt

for (( i=0;i<$num_sample;i++ ));
do
mkdir D${sample_list[${i}]}
cd D${sample_list[${i}]}
wot trajectory --tmap ../tmaps --cell_set ../../cell_cluster_DIC_fD5_iNSC_all.gmt --day ${sample_list[${i}]}
mv wot_trajectory.txt D${sample_list[${i}]}_wot_trajectory.txt
cd ../

done
