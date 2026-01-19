#!/usr/bin/env bash

ante-MMPBSA.py \
  -p system.prmtop \
  -c complex.prmtop \
  -r receptor.prmtop \
  -l ligand.prmtop \
  -s ":WAT,:Na+,:CL" \
  -n ":1-11" 

