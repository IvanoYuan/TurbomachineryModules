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

\*---------------------------------------------------------------------------*/

#include "pressureCCSDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cylindricalCS.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureCCSDirectedInletVelocityFvPatchVectorField::
pressureCCSDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    inletDir_(p.size()),
    cylindricalCCS_(0),
    n_(0, 0, 1),
    y_(1, 0, 0)
{}


pressureCCSDirectedInletVelocityFvPatchVectorField::
pressureCCSDirectedInletVelocityFvPatchVectorField
(
    const pressureCCSDirectedInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletDir_(ptf.inletDir_, mapper),
    cylindricalCCS_(ptf.cylindricalCCS_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


pressureCCSDirectedInletVelocityFvPatchVectorField::
pressureCCSDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    inletDir_("inletDirection", dict, p.size()),
    cylindricalCCS_(dict.lookup("cylindricalCCS")),
    n_(dict.lookup("n")),
    y_(dict.lookup("y"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    /*if (!inletDir_.empty())
    {
        if (min(mag(inletDir_)) < SMALL)
        {
            FatalErrorIn
            (
                "pressureCCSDirectedInletVelocityFvPatchVectorField::\n"
                "pressureCCSDirectedInletVelocityFvPatchVectorField\n"
                "(\n"
                "    const fvPatch& p,\n"
                "    const DimensionedField<vector, volMesh>& iF,\n"
                "    const dictionary& dict\n"
                ")"
            )    << "Badly defined inlet direction for field "
                 << this->dimensionedInternalField().name()
                 << " and patch " << this->patch().name()
                 << abort(FatalError);
        }
    }     */  

    // Normalise to obtain the flow direction
    //inletDir_ /= (mag(inletDir_) + SMALL);
}


pressureCCSDirectedInletVelocityFvPatchVectorField::
pressureCCSDirectedInletVelocityFvPatchVectorField
(
    const pressureCCSDirectedInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_),
    cylindricalCCS_(pivpvf.cylindricalCCS_),
    n_(pivpvf.n_),
    y_(pivpvf.y_)
{}


pressureCCSDirectedInletVelocityFvPatchVectorField::
pressureCCSDirectedInletVelocityFvPatchVectorField
(
    const pressureCCSDirectedInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_),
    cylindricalCCS_(pivpvf.cylindricalCCS_),
    n_(pivpvf.n_),
    y_(pivpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureCCSDirectedInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    inletDir_.autoMap(m);
}


void pressureCCSDirectedInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const pressureCCSDirectedInletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureCCSDirectedInletVelocityFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void pressureCCSDirectedInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!this->db().objectRegistry::found(phiName_))
    {
        // Flux not available, do not update
        InfoIn
        (
            "void pressureCCSDirectedInletVelocityFvPatchVectorField::"
            "updateCoeffs()"
        )   << "Flux field " << phiName_ << " not found.  "
            << "Performing fixed value update" << endl;

        fixedValueFvPatchVectorField::updateCoeffs();

        return;
    }
    
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);
    
    vectorField inletDir = inletDir_ /mag(inletDir_);// Normalise to obtain the flow direction
    vectorField n = patch().nf();
    scalarField ndmagS = (n & inletDir)*patch().magSf();
    
    

    if (cylindricalCCS_)
    {
        const vectorField& carc = patch().Cf(); //get centre point of every surface element
        
        cylindricalCS cyl( "cylindricalCS", point((carc[0]&n_),0,0), n_, y_, false ); //Construct cylidricalCS from origin and 2 axes
        
        vectorField cylc = cyl.localVector(carc); //location transformed into cylindrical coordinates
        
        //put radius and theta into scalarFields
        scalarField radius = cylc.component(vector::X);
        scalarField theta = cylc.component(vector::Y);
        
        //calculate cylindrical velocity
        vectorField inletDirCCS = inletDir;
        vector z_(0,0,0);
        z_ = n_ ^ y_;
        inletDirCCS.replace(vector::X, ((inletDir & y_)*cos(theta))+((inletDir & z_)*sin(theta)));//Vr
        inletDirCCS.replace(vector::Y, (-1*(inletDir & y_)*sin(theta))+((inletDir & z_)*cos(theta)));//Vteta
        inletDirCCS.replace(vector::Z, inletDir & n_); //Vz
        
        inletDir = inletDirCCS;
        
    }

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(inletDir*phip/ndmagS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        if (!this->db().objectRegistry::found(rhoName_))
        {
            // Rho not available, do not update
            fixedValueFvPatchVectorField::updateCoeffs();

            return;
        }

        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(inletDir*phip/(rhop*ndmagS));
    }
    else
    {
        FatalErrorIn
        (
            "pressureCCSDirectedInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void pressureCCSDirectedInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }
    inletDir_.writeEntry("inletDirection", os);
    os.writeKeyword("cylindricalCCS") << cylindricalCCS_ << token::END_STATEMENT << nl;
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y") << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void pressureCCSDirectedInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(inletDir_*(inletDir_ & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pressureCCSDirectedInletVelocityFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
