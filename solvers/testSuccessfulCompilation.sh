#!/bin/sh
bold=$(tput bold)
normal=$(tput sgr0)
echo "${bold}Testing successful Compilation of: \n"

echo "${bold}**************************** TEST **************************** \n\nadjointShapeOptimizationFoamPower\n\n${normal}"

adjointShapeOptimizationFoamPower

echo "${bold}**************************** TEST **************************** \n\nVCadjointShapeOptimizationFoam\n\n${normal}"

VCadjointShapeOptimizationFoam

echo "${bold}**************************** TEST **************************** \n\nVCadjointShapeOptimizationFoamPower\n\n${normal}"

VCadjointShapeOptimizationFoamPower

echo "${bold}**************************** TEST **************************** \n\nVCTopOpThermalFluid\n\n${normal}"

VCTopOpThermalFluid


echo "${bold}If you see the FOAM header after a solver is run, then the solver was compiled successfully.${normal}"
