#!/bin/sh

cd 1_toyProblem
./executeCases.sh

cd ..
cd 2_ductSystem_benchmark1
./executeCases.sh

cd ..
cd 3_ductSystemAlt_benchmark2
./executeCases.sh

cd ..
cd 4_dualBottleneck_benchmark3
./executeCases.sh

cd ..
cd 5_thermalFluid
./executeCases.sh
