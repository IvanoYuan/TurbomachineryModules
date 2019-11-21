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
    

Description
    

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "basicPsiThermo.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    
    bool writeResults = !args.optionFound("noWrite");

    IOobject UPtotheader
    (
        "UPtot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject UTtotheader
    (
        "UTtot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
    IOobject phiheader
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
    if (UPtotheader.headerOk() && UTtotheader.headerOk() && phiheader.headerOk())
    {
    	label inletID = mesh.boundaryMesh().findPatchID("inlet");
    	label outletID = mesh.boundaryMesh().findPatchID("outlet");
    	
    	volScalarField UPtot(UPtotheader, mesh);
        volScalarField UTtot(UTtotheader, mesh);
        surfaceScalarField phi(phiheader, mesh);
        
        if (isFile(runTime.constantPath()/"thermophysicalProperties"))
        {
            // thermophysical UPtot
            autoPtr<basicPsiThermo> thermo
            (
                basicPsiThermo::New(mesh)
            );

            volScalarField Cp = thermo->Cp();
            volScalarField Cv = thermo->Cv();
            scalar gamma = 1.4;
            
            scalar Pin = 0;
            scalar Pout = 0;
            scalar Tin = 0;
            scalar Tout = 0;
            
            scalar totalPressureRatio = 0;
            scalar totalTemperatureRatio = 0;
            scalar isentropicEfficiency = 0;
            scalar polytropicEfficiency = 0;
            
            Pin = gSum(UPtot.boundaryField()[inletID]*mag(phi.boundaryField()[inletID]))/gSum(mag(phi.boundaryField()[inletID]));
            Pout = gSum(UPtot.boundaryField()[outletID]*mag(phi.boundaryField()[outletID]))/gSum(mag(phi.boundaryField()[inletID]));
            Tin = gSum(UTtot.boundaryField()[inletID]*mag(phi.boundaryField()[inletID]))/gSum(mag(phi.boundaryField()[inletID]));
            Tout = gSum(UTtot.boundaryField()[outletID]*mag(phi.boundaryField()[outletID]))/gSum(mag(phi.boundaryField()[inletID]));
        
            totalPressureRatio = Pout/Pin;
            totalTemperatureRatio = Tout/Tin;
            isentropicEfficiency = (pow(totalPressureRatio,(gamma-1)/gamma)-1)/(totalTemperatureRatio-1);
            polytropicEfficiency = ((gamma-1)/gamma)*(log(totalPressureRatio))/(log(totalTemperatureRatio));
            
            Info << endl
                 << "turboPerformance"<< endl
                 << "totalPressureRatio = " << totalPressureRatio << endl
                 << "inletMassFlow = " << -36*sum(phi.boundaryField()[inletID]) << endl
                 << "outletMassFlow = " << 36*sum(phi.boundaryField()[outletID]) << endl
                 << "totalTemperatureRatio = " << totalTemperatureRatio << endl
                 << "isentropicEfficiency = " << isentropicEfficiency << endl
                 << "polytropicEfficiency = " << polytropicEfficiency << endl;
                 
           if (writeResults)
           {
            OFstream os(runTime.timeName()/"turboPerformance");
            os << "inletMassFlow = " << -36*sum(phi.boundaryField()[inletID]) << endl;
            os << "outletMassFlow = " << 36*sum(phi.boundaryField()[outletID]) << endl;
            os << "totalPressureRatio = " << totalPressureRatio << endl;
            os << "totalTemperatureRatio = " << totalTemperatureRatio << endl;
            os << "isentropicEfficiency = " << isentropicEfficiency << endl;
            os << "polytropicEfficiency = " << polytropicEfficiency << endl;
           }   
        }
    }
    else
    {
        Info<< "Fuck openfoam: Missing UPtot or UTtot or phi" << endl;
    }
}

// ************************************************************************* //
