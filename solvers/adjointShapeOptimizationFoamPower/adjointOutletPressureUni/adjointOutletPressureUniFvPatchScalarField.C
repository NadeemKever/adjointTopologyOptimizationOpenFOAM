/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "adjointOutletPressureUniFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletPressureUniFvPatchScalarField::
adjointOutletPressureUniFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::adjointOutletPressureUniFvPatchScalarField::
adjointOutletPressureUniFvPatchScalarField
(
    const adjointOutletPressureUniFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adjointOutletPressureUniFvPatchScalarField::
adjointOutletPressureUniFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::adjointOutletPressureUniFvPatchScalarField::
adjointOutletPressureUniFvPatchScalarField
(
    const adjointOutletPressureUniFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletPressureUniFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");

	//............Begin Mods............//

    scalarField Up_n = Up & patch().nf();   //Primal normal velocity

    scalarField Uap_n = Uap & patch().nf();  //Adjoint normal velocity

    const incompressible::RASModel& rasModel = db().lookupObject<incompressible::RASModel>("RASProperties");

    scalarField nueff = rasModel.nuEff()().boundaryField()[patch().index()];  //effective viscosity

    const scalarField& deltainv = patch().deltaCoeffs(); // 1/distance    dist = dist b/w patches

    scalarField Uaneigh_n = (Uap.patchInternalField() & patch().nf());

	const fvPatchField<vector>& Udp = patch().lookupPatchField<volVectorField, vector>("Ud");

	scalarField Udp_n = (Udp & patch().nf());

    operator== ((Uap&Up) + (Up_n * Uap_n) + nueff*deltainv* (Uap_n - Uaneigh_n) + (Up_n - Udp_n)) ;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointOutletPressureUniFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointOutletPressureUniFvPatchScalarField
    );
}

// ************************************************************************* //
