#!/bin/bash

NN_K=`grep "nn_K       =" calculate_melt_rate_ensemble_DeltaT_IDEALIZED.f90 | awk '{print $3}'`
echo "nn_K = $NN_K"

rm -f gxtm_*

# SPLIT JOBS FOR ICE SHELVES FOR WHICH WE CANNOT CALCULATE THE ENSEMBLE IN 24h walltime :

for ID in 10 21 11 # ROSS, RONNE, FILCHNER
do
  rm -f gxtm_ID${ID}_[0-9]*
  if [ $ID -eq 10 ]; then
    walltime=24 # for nn_K=19
  elif [ $ID -eq 21 ]; then
    walltime=17 # for nn_K=19
  else
    walltime=06
  fi
  for ii in $(seq 1 ${NN_K})
  do
    if [ $ID -eq 10 ] && [ $ii -ge 15 ]; then
      sed -e "s/   DO kk_K=1,nn_K/   DO kk_K=${ii},${ii}/" \
          -e "s/'melt_ensemble_DeltaT_IDEALIZED_',i3.3,'.nc'/'melt_ensemble_DeltaT_IDEALIZED_',i3.3,'_${ii}.nc'/" \
          -e "s/Ntun = 50/Ntun = 45/" \
          -e "s/DO kisf=2,mNisf/DO kisf=${ID},${ID}/" calculate_melt_rate_ensemble_DeltaT_IDEALIZED.f90 > gxtm_ID${ID}_${ii}.f90
    else
      sed -e "s/   DO kk_K=1,nn_K/   DO kk_K=${ii},${ii}/" \
          -e "s/'melt_ensemble_DeltaT_IDEALIZED_',i3.3,'.nc'/'melt_ensemble_DeltaT_IDEALIZED_',i3.3,'_${ii}.nc'/" \
          -e "s/DO kisf=2,mNisf/DO kisf=${ID},${ID}/" calculate_melt_rate_ensemble_DeltaT_IDEALIZED.f90 > gxtm_ID${ID}_${ii}.f90
    fi
    ifort -qopenmp -c $NC_INC gxtm_ID${ID}_${ii}.f90
    ifort -qopenmp -o gxtm_ID${ID}_${ii} gxtm_ID${ID}_${ii}.o $NC_LIB
    if [ $ID -eq 10 ]; then
      ./submit_OMP12.sh gxtm_ID${ID}_${ii} $walltime 64
    else
      ./submit_OMP.sh gxtm_ID${ID}_${ii} $walltime 64
    fi
  done
  rm -f gxtm_ID${ID}_${ii}.f90 gxtm_ID${ID}_${ii}.o
done

# THEN SPLIT JOBS FOR GROUPS OF ICE SHELVES THAT CAN BE CALCULATED IN 24h walltime :

for GROUP in $(seq 1 8)
do
  if [ $GROUP -eq 1 ]; then
    isf1=12
    isf2=20
  elif [ $GROUP -eq 2 ]; then
    isf1=22
    isf2=30
  elif [ $GROUP -eq 3 ]; then
    isf1=31
    isf2=38
  elif [ $GROUP -eq 4 ]; then
    isf1=39
    isf2=42
  elif [ $GROUP -eq 5 ]; then
    isf1=43
    isf2=51
  elif [ $GROUP -eq 6 ]; then
    isf1=52
    isf2=63
  elif [ $GROUP -eq 7 ]; then
    isf1=64
    isf2=70
  elif [ $GROUP -eq 8 ]; then
    isf1=71
    isf2=75
  fi
    sed -e "s/DO kisf=2,mNisf/DO kisf=${isf1},${isf2}/" calculate_melt_rate_ensemble_DeltaT_IDEALIZED.f90 > gxtm_${GROUP}.f90
    ifort -qopenmp -c $NC_INC gxtm_${GROUP}.f90
    ifort -qopenmp -o gxtm_${GROUP} gxtm_${GROUP}.o $NC_LIB
    ./submit_OMP.sh gxtm_${GROUP} 24 64
    rm -f gxtm_${GROUP}.f90 gxtm_${GROUP}.o
done
