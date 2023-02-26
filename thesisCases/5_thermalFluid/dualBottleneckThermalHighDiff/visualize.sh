#!/bin/bash

grep 'Objective Function (Power Dissipated) J:' log.txt > objectiveVal.txt
sed -i s/'Objective Function (Power Dissipated) J: '/''/g objectiveVal.txt

python visualize.py

xdg-open Objective.png
