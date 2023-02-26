#!/bin/bash


cd ductSystem2Power
VCadjointShapeOptimizationFoamPower > log.txt &
paraFoam -touch
