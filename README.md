# adjointTopologyOptimizationOpenFOAM
Implementation of a volume constrained topology optimization solver using pressure losses or power dissipation as an objective function, with example cases.

For a detailed tutorial, please see Appendix C of the provided thesis found [here](https://github.com/NadeemKever/adjointTopologyOptimizationOpenFOAM/blob/main/Master_Thesis_Kever_Nadeem.pdf). 

## Quick Guide:
To get up and running quickly, enter the `solvers` directory and execute `$./compileSolvers.sh` to compile all of the solvers in the thesis. You can quickly test that all solvers have been successfully compiled by running `$./testSuccessfulCompilation.sh`. If everything has compiled successfully, then you should see the name of the specific solver followed with the OpenFOAM header. After verifying this, navigate to the `thesisCases` directory and execute `$./runAllCases.sh`. This should solve all the cases discussed in the thesis and provide output files that can be visualized.


