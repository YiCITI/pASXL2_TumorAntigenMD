#!/usr/bin/env bash

traj=$1
prefix=$(basename "$traj")
prefix=${prefix%.*}

ulimit -v $((64 * 1024 * 1024))

nohup MMPBSA.py -O \
  -i mmpbsa_pb.in \
  -cp complex.prmtop \
  -rp receptor.prmtop \
  -lp ligand.prmtop \
  -y "$traj" \
  -o "FINAL_RESULTS_MMPBSA_PB_${prefix}.dat" \
  -eo "FINAL_RESULTS_MMPBSA_PB_${prefix}.csv" \
  > "mmpbsa_pb_${prefix}.log" 2>&1 &
