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

Application
    UPtot

Description
    Calculates and optionally writes the local UPtot number from the velocity
    field U at each time.

    The -nowrite option just outputs the max value without writing the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "basicPsiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject Theader
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
    IOobject pheader
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check U and T exists
    if (Uheader.headerOk() && Theader.headerOk() && pheader.headerOk())
    {
        autoPtr<volScalarField> UPtotPtr;
        
        autoPtr<volScalarField> UTtotPtr;

        volVectorField U(Uheader, mesh);
        
        volScalarField p(pheader, mesh);

        if (isFile(runTime.constantPath()/"thermophysicalProperties"))
        {
            // thermophysical UPtot
            autoPtr<basicPsiThermo> thermo
            (
                basicPsiThermo::New(mesh)
            );

            volScalarField Cp = thermo->Cp();
            volScalarField Cv = thermo->Cv();
            
            volScalarField Ma = mag(U)/(sqrt((Cp/Cv)*(Cp - Cv)*thermo->T()));
            
            volScalarField a = 1+0.5*(Cp/Cv-1)*sqr(Ma);

            UPtotPtr.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "UPtot",
                        runTime.timeName(),
                        mesh
                    ),
                    p*pow(a, Cp/(Cp-Cv))
                )
            );
            
            UTtotPtr.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "UTtot",
                        runTime.timeName(),
                        mesh
                    ),
                    a*thermo->T()
                )
            );
        }
        /*
        else
        {
            // thermodynamic UPtot
            IOdictionary thermoProps
            (
                IOobject
                (
                    "thermodynamicProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            dimensionedScalar R(thermoProps.lookup("R"));
            dimensionedScalar Cv(thermoProps.lookup("Cv"));

            volScalarField T(Theader, mesh);

            UPtotPtr.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "Ma",
                        runTime.timeName(),
                        mesh
                    ),
                    mag(U)/(sqrt(((Cv + R)/Cv)*R*T))
                )
            );
        }
        */

        Info<< "UPtot max : " << max(UPtotPtr()).value() << endl
        << "UPtot min : " << min(UPtotPtr()).value() << endl;

        if (writeResults)
        {
            UPtotPtr().write();
            UTtotPtr().write();
        }
    }
    else
    {
        Info<< "    Missing U or T or p" << endl;
    }
}


// ************************************************************************* //
