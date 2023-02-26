// alpha or DT, thermal diffusivity coefficient = k/rho/Cp
// Cp, specific heat capacity, for example PVC: Cp = 1530 J/kg/K
// k, thermal conductivity, for example PVC: k = 0.166 W/K/m
// rho, density, for example PVC: 1380 kg/m^3

// DT = k/rho/Cp = 7.8621e-8


// DT and Cp are read from transportProperties, see createFields.H



fvScalarMatrix TEqn
(
        //fvm::ddt(T) // <--- steady state simulation, so ddt disabled...
    //+
        fvm::div(phi, T)
      - fvm::laplacian(DT, T)
);

    TEqn.relax();


    TEqn.solve().initialResidual();