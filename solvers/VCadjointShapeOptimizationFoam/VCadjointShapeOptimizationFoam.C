/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ajointShapeOptimizationFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with optimisation of duct shape by applying "blockage" in regions
    causing pressure loss as estimated using an adjoint formulation.

    References:
    \verbatim
        "Implementation of a continuous adjoint for topology optimization of
         ducted flows"
        C. Othmer,
        E. de Villiers,
        H.G. Weller
        AIAA-2007-3947
        http://pdf.aiaa.org/preview/CDReadyMCFD07_1379/PV2007_3947.pdf
    \endverbatim

    Note that this solver optimises for total pressure loss whereas the
    above paper describes the method for optimising power-loss.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flow"
        " of non-Newtonian fluids with duct shape optimisation"
        " by applying 'blockage' in regions causing pressure loss"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    int count = 0;
    while (simple.loop())
    {   
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "Count = " << count << nl << endl;

        //alpha +=
        //    mesh.relaxationFactor("alpha")
        //   *(lambda*max(Ua & U, zeroSensitivity) - alpha);
        // if (count < 1000){
        //     alpha += mesh.fieldRelaxationFactor("alpha")
        //         *(min(max(alpha + lambda*(Ua & U), zeroAlpha), alphaMax) - alpha);}
        // 
        alpha += mesh.fieldRelaxationFactor("alpha")
            *(min(max(alpha + lambda*(Ua & U), zeroAlpha), alphaMax) - alpha);
        zeroCells(alpha, inletCells);
        zeroCells(alpha, outletCells);

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              + fvm::Sp(alpha, U)
             ==
                fvOptions(U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvOptions.correct(U);

            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            tUEqn.clear();
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }

        // Adjoint Pressure-velocity SIMPLE corrector
        {
            // Adjoint Momentum predictor

            volVectorField adjointTransposeConvection((fvc::grad(Ua) & U));
            //volVectorField adjointTransposeConvection
            //(
            //    fvc::reconstruct
            //    (
            //        mesh.magSf()*fvc::dotInterpolate(fvc::snGrad(Ua), U)
            //    )
            //);

            zeroCells(adjointTransposeConvection, inletCells);

            tmp<fvVectorMatrix> tUaEqn
            (
                fvm::div(-phi, Ua)
              - adjointTransposeConvection
              + turbulence->divDevReff(Ua)
              + fvm::Sp(alpha, Ua)
             ==
                fvOptions(Ua)
            );
            fvVectorMatrix& UaEqn = tUaEqn.ref();

            UaEqn.relax();

            fvOptions.constrain(UaEqn);

            solve(UaEqn == -fvc::grad(pa));

            fvOptions.correct(Ua);

            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HbyAa("HbyAa", Ua);
            HbyAa = rAUa*UaEqn.H();
            tUaEqn.clear();
            surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
            adjustPhi(phiHbyAa, Ua, pa);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
                );

                paEqn.setReference(paRefCell, paRefValue);
                paEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phia = phiHbyAa - paEqn.flux();
                }
            }

            #include "adjointContinuityErrs.H"

            // Explicitly relax pressure for adjoint momentum corrector
            pa.relax();

            // Adjoint momentum corrector
            Ua = HbyAa - rAUa*fvc::grad(pa);
            Ua.correctBoundaryConditions();
            fvOptions.correct(Ua);
        }

        laminarTransport.correct();
        turbulence->correct();







        sens = Ua & U;
		

        double maxSens = -10000000.0;
        double minSens =  10000000.0;
        
        double cellSens = 0.0;
        double totalVol = 0.0;
        forAll(mesh.C(), ID) {
            totalVol += mesh.V()[ID];
            cellSens = sens[ID];
            //The negative is because ... look at the plot and look at the pt of high and low sens really neg --> no alpha
            if (cellSens > maxSens) {
                maxSens = cellSens;
            }
            if (cellSens < minSens) {
                minSens = cellSens;
            }

//			intGammaVol += gamma[ID]*mesh.V()[ID];

        }
        Info<<"Total Volume:  " <<totalVol<<endl;
        Info<<"Max Sensitivity: "<<maxSens<<endl;
        Info<<"Min Sensitivity: "<<minSens<<endl;
        if (count >= iterationToEnforceVolumeConstraint.value() ){
            int steps = 250;
            double x[steps+1];
            int indx = 0;
            double step;
            step = (maxSens-minSens)/double(steps);
            Info<<"Step Size: "<<step<<endl;
            for (double i = maxSens; i > minSens; i-= step){
                
                x[indx] = i;
                indx += 1;
                // Info<<"Step: "<<indx<<'\n'<<"Sensitivity array:  " <<i<<endl;
            }
            


            double targetVol = volumeFraction.value() * totalVol;		
            double vol, sensL;
		
            for (int j = 0; j < steps; j++){
                vol = 0.0;
                sensL = x[j];
                Info<<"Loop Iteration:  " <<j<<endl;
                forAll(mesh.C(), ID) {
                    cellSens = sens[ID];
    //				alpha[ID] = pos(cellSens - sensL)*alphaMax.value();
    //				alphaTmp = pos(cellSens - sensL)*alphaMax.value();
                    if (cellSens > sensL){ 
                        alpha[ID] = alphaMax.value();
                    }
                    // else if (strictVolumeConstraint.value() == 1.0 ){
                    //     // Info<<"strict constraint hit!"<<endl;
                    //     alpha[ID] = 0.0;
                    //     // Info<<"alpha: "<<alpha[ID]<<"Cell: "<<ID<<endl;
                    // }
                }
                forAll(mesh.C(), ID) {
                    vol += (alpha[ID]/alphaMax.value())*mesh.V()[ID];
                }
                Info<<"Volume:  " <<vol<<endl;
                if (vol >= targetVol){
                    Info<<"Volume:  " <<vol<<endl;
                    break;
                }
            }
        }
		zeroCells(alpha, inletCells);
        zeroCells(alpha, outletCells);
		count += 1;
		
				



        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
