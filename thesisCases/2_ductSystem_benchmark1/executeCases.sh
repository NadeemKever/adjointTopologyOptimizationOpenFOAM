#!/bin/bash


cd ductSystemPowerVC
VCadjointShapeOptimizationFoamPower > log.txt &
paraFoam -touch

cd ../ductSystemPressureVC
VCadjointShapeOptimizationFoam > log.txt &
paraFoam -touch
