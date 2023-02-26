#!/bin/bash

cd adjointShapeOptimizationFoamPower

wclean; wmake>makeLog.txt 

cd ../VCadjointShapeOptimizationFoam

wclean; wmake>makeLog.txt &

cd ../VCadjointShapeOptimizationFoamPower

wclean; wmake>makeLog.txt &

cd ../VCTopOpThermalFluid

wclean; wmake>makeLog.txt 

cd ..


