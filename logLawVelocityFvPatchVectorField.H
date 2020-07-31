/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::velocityABLWallFunctionFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Boundary condition specifies a logarithmic velocity inlet profile
    (fixed value), given mean friction velocity, aerodynamic roughness and von-
    Karman constant, patch normal velocity at wall-normal heights are eastimated
    using the log-law of the wall.

SourceFiles
    logLawVelocityFvPatchVectorField.C

Author
    Rajib Roy, University of Wyoming, roy.rajib@live.com
\*---------------------------------------------------------------------------*/

#ifndef logLawVelocityFvPatchVectorField_H
#define logLawVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class logLawVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class logLawVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Friction velocity
        scalar uStar_;

	//- Aerodynamic roughness height
	scalar z0_;

	//- Von-Karman Constant, kappa
	scalar kappa_;

        //- Flow direction
        vector n_;

        //- Velocity Profile Direction (here the y-coordinate)
        vector h_;

        // reference height of the terrain
        vector refHeight_;


public:

    //- Runtime type information
    TypeName("logLawVelocity");


    // Constructors

        //- Construct from patch and internal field
        logLawVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        logLawVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given logLawVelocityFvPatchVectorField
        //  onto a new patch
        logLawVelocityFvPatchVectorField
        (
            const logLawVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new logLawVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        logLawVelocityFvPatchVectorField
        (
            const logLawVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new logLawVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Return friction velocity
        scalar& uStar()
        {
            return uStar_;
        }

        //- Return aerodynamic roughness height
        scalar& z0()
        {
            return z0_;
        }

        //- Return kappa
        scalar& kappa()
        {
            return kappa_;
        }


        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return height (h) direction
        vector& h()
        {
            return h_;
        }

        //- Return height (h) direction
        vector& refHeight()
        {
            return refHeight_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
