/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.
    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "logLawVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

logLawVelocityFvPatchVectorField::logLawVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    uStar_(0.0),
    z0_(0.0001),
    kappa_(0.40),
    n_(1, 0, 0),
    h_(0, 0, 1),
    refHeight_(0.0, 0.0, 0.0)
{}


logLawVelocityFvPatchVectorField::logLawVelocityFvPatchVectorField
(
    const logLawVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    uStar_(ptf.uStar_),
    z0_(ptf.z0_),
    kappa_(ptf.kappa_),
    n_(ptf.n_),
    h_(ptf.h_),
    refHeight_(ptf.refHeight_)
{}


logLawVelocityFvPatchVectorField::logLawVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    uStar_(readScalar(dict.lookup("uStar"))),
    z0_(readScalar(dict.lookup("z0"))),
    kappa_(readScalar(dict.lookup("kappa"))),
    n_(dict.lookup("n")),
    h_(dict.lookup("h")),
    refHeight_(dict.lookupOrDefault<vector>("refHeight", vector::zero))
{
    if (mag(n_) < SMALL || mag(h_) < SMALL)
    {
        FatalErrorIn("logLawVelocityFvPatchVectorField(dict)")
            << "n or h given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    h_ /= mag(h_);

    evaluate();
}


logLawVelocityFvPatchVectorField::logLawVelocityFvPatchVectorField
(
    const logLawVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    uStar_(fcvpvf.uStar_),
    z0_(fcvpvf.z0_),
    kappa_(fcvpvf.kappa_),
    n_(fcvpvf.n_),
    h_(fcvpvf.h_),
    refHeight_(fcvpvf.refHeight_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void logLawVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get range and orientation
    boundBox bb(patch().patch().localPoints(), true);

    // Does not need the patch center here
    // vector ctr = 0.5*(bb.max() + bb.min());

    // Face center vector field
    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate alog the vertical direction
    // scalarField coord = (c & h_);
    scalarField coord = (c - refHeight()) & h_;

    vectorField::operator=(n_*(uStar_/kappa_)*(Foam::log(coord/z0_)));
}


// Write
void logLawVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("uStar") << uStar_ << token::END_STATEMENT << nl;
    os.writeKeyword("z0")    << z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")     << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("h")     << h_ << token::END_STATEMENT << nl;
    os.writeKeyword("refHeight") << refHeight_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, logLawVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
