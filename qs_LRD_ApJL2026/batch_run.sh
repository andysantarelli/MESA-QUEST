#!/bin/bash

# Batch run script for quasi-star models
for mass in 1d4 1d5 1d6; do
    for z in 0.0002 0.000002; do
        shmesa change inlist_project log_directory "'M${mass}_z${z}'" initial_mass "${mass}" Zbase "${z}" initial_z "${z}" && 
        ./rn &&
        echo "Completed run for mass = ${mass}, z = ${z}"
    done
done