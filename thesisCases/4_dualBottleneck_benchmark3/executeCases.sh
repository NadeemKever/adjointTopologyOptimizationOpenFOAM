#!/bin/bash


cd dualBottleneckPower
VCadjointShapeOptimizationFoamPower > log.txt &
paraFoam -touch
