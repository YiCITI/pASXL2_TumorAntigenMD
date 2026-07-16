#!/usr/bin/env bash

traj=$1
prefix=${traj%.*}

ulimit -v $((64 * 1024 * 1024))

nohup MMPBSA.py -O \
  -i mmpbsa_gb.in \
  -cp complex.prmtop \
  -rp receptor.prmtop \
  -lp ligand.prmtop \
  -y "$traj" \
  -o "FINAL_RESULTS_MMPBSA_GB_${prefix}.dat" \
  -eo "FINAL_RESULTS_MMPBSA_GB_${prefix}.csv" \
  > "mmpbsa_gb_${prefix}.log" 2>&1 &
