#!/bin/bash


cd adjointShapeOptPower
adjointShapeOptimizationFoamPower > log.txt &
paraFoam -touch


cd ../VCadjointShapeOptPower70
VCadjointShapeOptimizationFoamPower > log.txt &
paraFoam -touch
cd ../VCadjointShapeOptPower80
VCadjointShapeOptimizationFoamPower > log.txt &
paraFoam -touch
cd ../VCadjointShapeOptPower90
VCadjointShapeOptimizationFoamPower > log.txt &
paraFoam -touch
