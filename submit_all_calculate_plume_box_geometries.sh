#!/bin/bash

rm -f trpj_*

# THEN SPLIT JOBS FOR GROUPS OF ICE SHELVES THAT CAN BE CALCULATED IN 24h walltime :

for GROUP in $(seq 1 9)
do
  if [ $GROUP -eq 1 ]; then
    isf1=10
    isf2=10
  elif [ $GROUP -eq 2 ]; then
    isf1=11
    isf2=20
  elif [ $GROUP -eq 3 ]; then
    isf1=21
    isf2=21
  elif [ $GROUP -eq 4 ]; then
    isf1=22
    isf2=22
  elif [ $GROUP -eq 5 ]; then
    isf1=23
    isf2=30
  elif [ $GROUP -eq 6 ]; then
    isf1=31
    isf2=50
  elif [ $GROUP -eq 7 ]; then
    isf1=51
    isf2=70
  elif [ $GROUP -eq 8 ]; then
    isf1=71
    isf2=130
  elif [ $GROUP -eq 9 ]; then
    isf1=131
    isf2=216
  fi
  # reduce nb of experiments here, the melt netcdf files will just be created for a quick check
  sed -e "s/DO kisf=2,mNisf/DO kisf=${isf1},${isf2}/" \
      -e "s#^nn_TS_pres =#nn_TS_pres = 5 !#" \
      -e "s#^nn_tuning  =#nn_tuning  = 1 !#" calculate_melt_rate_ensemble.f90 > trpj_${GROUP}.f90
  ifort -qopenmp -c $NC_INC trpj_${GROUP}.f90
  ifort -qopenmp -o trpj_${GROUP} trpj_${GROUP}.o $NC_LIB
  ./submit_OMP.sh trpj_${GROUP} 24 64
  rm -f trpj_${GROUP}.f90 trpj_${GROUP}.o
done
