/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Class
    Foam::radialEquilibriumPressureOutletFvPatchScalarField

Description
    

SourceFiles
    radialEquilibriumPressureOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef radialEquilibriumPressureOutletFvPatchScalarField_H
#define radialEquilibriumPressureOutletFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class radialEquilibriumPressureOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class radialEquilibriumPressureOutletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;
        
	scalar rRef_;

	scalar pRef_;

        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;



public:

    //- Runtime type information
    TypeName("radialEquilibriumPressureOutlet");


    // Constructors

        //- Construct from patch and internal field
        radialEquilibriumPressureOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        radialEquilibriumPressureOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given radialEquilibriumPressureOutletFvPatchScalarField
        //  onto a new patch
        radialEquilibriumPressureOutletFvPatchScalarField
        (
            const radialEquilibriumPressureOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        radialEquilibriumPressureOutletFvPatchScalarField
        (
            const radialEquilibriumPressureOutletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new radialEquilibriumPressureOutletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        radialEquilibriumPressureOutletFvPatchScalarField
        (
            const radialEquilibriumPressureOutletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new radialEquilibriumPressureOutletFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

        const word& UName() const
        {
            return UName_;
        }
        
	scalar& rRef()
        {
            return rRef_;
        }
        
        scalar& pRef()
        {
            return pRef_;
        }

        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return y direction
        vector& y()
        {
            return y_;
        }
        
        // Mapping functions
/*
            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );
            */
        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
