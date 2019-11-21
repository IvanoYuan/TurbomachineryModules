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

#include "fvCFD.H"

#include "cylindricalCS.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    word fieldName(args.additionalArgs()[0]);
    word patchName(args.additionalArgs()[1]);
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        
        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        mesh.readUpdate();
        
    	label patchID = mesh.boundaryMesh().findPatchID(patchName);
    	const fvPatch& cPatch = mesh.boundary()[patchID];
    	
    	//Info << cPatch.nf() << endl;
    	//get centre point of every surface element
    	const vectorField& carc = cPatch.Cf();
    	
    	//volVectorField field(fieldHeader, mesh);
    	
	
	//Construct cylidricalCS from origin and 2 axes
     	cylindricalCS cyl( "cylindricalCS", point(0.1067,0,0), vector(0,0,1), vector(1,0,0), false );
	
    	//location transformed into cylindrical coordinates
	vectorField cylc = cyl.localVector(carc);
	
	//put radius and theta into scalarFields
	scalarField radius = cylc.component(vector::X);
	scalarField theta = cylc.component(vector::Y);
	
    	//general way to get thetaN and radiusN
    	//1-initialize a structure
    	struct location
    	{
    		int num;
   	 	scalar radiusX;
    		scalar thetaY;
    	}loc[theta.size()],temp;
    	//2-assign the struct location
    	for(int i=0; i<theta.size(); i++)
    	{
    		loc[i].num = i;
    		loc[i].radiusX = radius[i];
    		loc[i].thetaY = theta[i];
    	}
    	//3-sort based on radius and theta
	for(int i=0;i<(theta.size()-1);i++)//Bubble sort
	{                     
		for(int j=0;j<(theta.size()-i-1);j++)
		{
			if(loc[j].radiusX > loc[j+1].radiusX)
			{
				temp=loc[j];
				loc[j]=loc[j+1];
				loc[j+1]=temp;
			}
		}
	}
	
    	//calculate thetaN
    	int thetaN = 1;
    	for(int i = 0; i<(theta.size()); i++)
    	{
    		//Info << loc[i].radiusX << endl; 
    		if(mag(loc[i+1].radiusX-loc[i].radiusX)<=1e-6)
    		{
    			++thetaN;
   		}
   		else
   		{
   			break;
   		}
    	}
    	
    	//finally, output thetaN and new point order
    	//fileName outputFile("pointOrdering");
    	OFstream os(runTime.constant()/"pointOrdering");
    	os << "(";
    	os << endl;
    	for(int i=0; i<(theta.size()+1); i++)
    	{
    		if(i==0)
    		{
    			os << thetaN;
    			os << endl;
    		}
    		else
    		{
    			os << loc[i-1].num;// << "=" << loc[i-1].radiusX;
    			os << endl;
    		}
    	}
    	os << ")";
    	os << endl;
    }
}



// ************************************************************************* //
