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

Author
    Aeron Yuan

\*---------------------------------------------------------------------------*/

#include "compressibleTotalTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleTotalTemperature, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        compressibleTotalTemperature,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleTotalTemperature::compressibleTotalTemperature
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    UName_(dict.lookup("UName"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating compressibleTotalTemperature from field "
        << UName_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::compressibleTotalTemperature::start()
{
    return true;
}


bool Foam::compressibleTotalTemperature::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    if
    (
        mesh.foundObject<volVectorField>(UName_)
     && mesh.foundObject<basicThermo>("thermophysicalProperties")
    )
    {
        const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);
        
        const basicThermo& thermo =
            mesh.lookupObject<basicThermo>("thermophysicalProperties");

        volScalarField Cp = thermo.Cp();
        volScalarField Cv = thermo.Cv();
        volScalarField T = thermo.T();
        
        volScalarField gamma = Cp/Cv;
        volScalarField R = Cp - Cv;
	volScalarField Ma = mag(U)/sqrt(gamma*R*T);
	
	
        volScalarField Ttot
        (
            IOobject
            (
                UName_ + "Ttot",
                U.instance(),
                mesh
            ),
            T*(1+0.5*(gamma-1)*sqr(Ma))
        );

        Info<< "compressibleTotalTemperature min = " << Foam::min(Ttot.internalField())
            << " max = " << Foam::max(Ttot.internalField()) << endl;

        if (mesh.time().outputTime())
        {
            Info << "Writing compressibleTotalTemperature field" << endl;
            Ttot.write();
        }

        return true;
    }
    else
    {
        Info<< "Field "  << UName_ << " or thermo not found.  Skipping."
            << endl;

        return false;
    }
}


bool Foam::compressibleTotalTemperature::read(const dictionary& dict)
{
    UName_ = word(dict.lookup("UName"));

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return false;
}

// ************************************************************************* //
