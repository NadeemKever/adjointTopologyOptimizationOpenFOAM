#!/bin/bash


cd dualBottleneckThermalHighDiff
VCTopOpThermalFluid > log.txt &
paraFoam -touch

cd ../dualBottleneckThermalLowDiff
VCTopOpThermalFluid > log.txt &
paraFoam -touch
