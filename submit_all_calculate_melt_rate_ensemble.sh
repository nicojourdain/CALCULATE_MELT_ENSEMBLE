#!/bin/bash

NN_TUNE=`grep "nn_tuning  =" calculate_melt_rate_ensemble.f90 | awk '{print $3}'`
echo "nn_tuning = $NN_TUNE"

rm -f trpj_*

# SPLIT JOBS FOR ICE SHELVES FOR WHICH WE CANNOT CALCULATE THE ENSEMBLE IN 24h walltime :

for ID in 10 21 11 # ROSS, RONNE, FILCHNER
do
  rm -f trpj_ID${ID}_[0-9]*
  if [ $ID -eq 10 ]; then
    walltime=10 # for nn_tuning=19
  elif [ $ID -eq 21 ]; then
    walltime=06 # for nn_tuning=19
  else
    walltime=03
  fi
  for ii in $(seq 1 ${NN_TUNE})
  do
    sed -e "s/   DO kk_tuning=1,nn_tuning/   DO kk_tuning=${ii},${ii}/" \
        -e "s/'melt_ensemble_',i3.3,'.nc'/'melt_ensemble_',i3.3,'_${ii}.nc'/" \
        -e "s/DO kisf=2,mNisf/DO kisf=${ID},${ID}/" calculate_melt_rate_ensemble.f90 > trpj_ID${ID}_${ii}.f90
    ifort -qopenmp -c $NC_INC trpj_ID${ID}_${ii}.f90
    ifort -qopenmp -o trpj_ID${ID}_${ii} trpj_ID${ID}_${ii}.o $NC_LIB
    ./submit_OMP.sh trpj_ID${ID}_${ii} $walltime 64
  done
  rm -f trpj_ID${ID}_${ii}.f90 trpj_ID${ID}_${ii}.o
done

# THEN SPLIT JOBS FOR GROUPS OF ICE SHELVES THAT CAN BE CALCULATED IN 24h walltime :

for GROUP in $(seq 1 6)
do
  if [ $GROUP -eq 1 ]; then
    isf1=12
    isf2=20
  elif [ $GROUP -eq 2 ]; then
    isf1=22
    isf2=30
  elif [ $GROUP -eq 3 ]; then
    isf1=31
    isf2=42
  elif [ $GROUP -eq 4 ]; then
    isf1=43
    isf2=53
  elif [ $GROUP -eq 5 ]; then
    isf1=54
    isf2=65
  elif [ $GROUP -eq 6 ]; then
    isf1=66
    isf2=75
  fi
    sed -e "s/DO kisf=2,mNisf/DO kisf=${isf1},${isf2}/" calculate_melt_rate_ensemble.f90 > trpj_${GROUP}.f90
    ifort -qopenmp -c $NC_INC trpj_${GROUP}.f90
    ifort -qopenmp -o trpj_${GROUP} trpj_${GROUP}.o $NC_LIB
    ./submit_OMP.sh trpj_${GROUP} 12 64
    rm -f trpj_${GROUP}.f90 trpj_${GROUP}.o
done
