Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhia.H"


label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.setFluxRequired(pa.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


dimensionedScalar zeroSensitivity(dimVelocity*dimVelocity, Zero);
dimensionedScalar zeroAlpha(dimless/dimTime, Zero);

dimensionedScalar lambda
(
    "lambda",
    dimTime/sqr(dimLength),
    laminarTransport
);

dimensionedScalar alphaMax
(
    "alphaMax",
    dimless/dimTime,
    laminarTransport
);
dimensionedScalar volumeFraction
(
    "volumeFraction",
    dimless,
    laminarTransport
);
dimensionedScalar iterationToEnforceVolumeConstraint
(
    "iterationToEnforceVolumeConstraint",
    dimless,
    laminarTransport
);

const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
//const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    lambda*max(Ua & U, zeroSensitivity)
);
zeroCells(alpha, inletCells);
//zeroCells(alpha, outletCells);

#include "createFvOptions.H"

dictionary optFunc = mesh.solutionDict().subDict("objectiveFunctionDict");
///////////////////////////// Cost Function ///////////////////////////////////

int nObjPatch = optFunc.lookupOrDefault<scalar>("numberObjectivePatches",0);

int objFunction = optFunc.lookupOrDefault<scalar>("objectiveFunction",0);

wordList objPatchNames = optFunc.lookup("objectivePatchesNames");

Info<<"Initializing objective function calculation:"<<endl ;
Info<<"The objective function chosen is:"<<objFunction<<endl ;
Info<<"Name of the patches for which the cost function will be calculated:"<<objPatchNames<<endl;
Info<<"Number of patches:"<<nObjPatch<<endl ;


label objPatchList[nObjPatch];
int iLoop;
for (iLoop=0; iLoop < nObjPatch; iLoop++)
	{ 
		objPatchList[iLoop] = mesh.boundaryMesh().findPatchID(objPatchNames[iLoop]);
	}


///////////////////////////////////////////////////////////////////////////////

vector Udlookup = optFunc.lookupOrDefault<vector>("Uduserdefnodim",vector(1,0,0));

volVectorField Ud
(
	IOobject
		(
			"Ud",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
	mesh,
	Udlookup
);

volScalarField sens
	(
	IOobject
		(
			"sensitivity",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
	Ua&U
	);
			
