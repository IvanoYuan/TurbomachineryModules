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

#include "radialEquilibriumPressureOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "cylindricalCS.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
radialEquilibriumPressureOutletFvPatchScalarField::
radialEquilibriumPressureOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    rRef_(0),
    pRef_(101325),
    n_(0, 0, 1),
    y_(1, 0, 0)
    //thetaN(1)
{}

Foam::
radialEquilibriumPressureOutletFvPatchScalarField::
radialEquilibriumPressureOutletFvPatchScalarField
(
    const radialEquilibriumPressureOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rRef_(ptf.rRef_),
    pRef_(ptf.pRef_),
    n_(ptf.n_),
    y_(ptf.y_)
    //thetaN(ptf.thetaN)
{}

Foam::
radialEquilibriumPressureOutletFvPatchScalarField::
radialEquilibriumPressureOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rRef_(readScalar(dict.lookup("rRef"))),
    pRef_(readScalar(dict.lookup("pRef"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y"))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::updateCoeffs();
	//operator==(scalarField("pRef", dict, p.size()));
    }
    // get span bounds and check whether rRef is out of range
    /*const pointField& pp = patch().patch().localPoints();
    scalarField ppRadius = mag(pp - (pp & n_)*n_);
    scalar ppRmin = min(ppRadius);
    scalar ppRmax = max(ppRadius);
    //Info << ppRmin << endl<< ppRmax << endl;
    if( rRef_<ppRmin || rRef_>ppRmax )
    {
    	FatalErrorIn
        (
            "void radialEquilibriumPressureOutletFvPatchScalarField::updateCoeffs()"
        )   << " rRef= " << rRef_ << " is out of  "
            << patch().name() << " range "
            << abort(FatalError);
    }
    //operator==(scalarField("value", dict, p.size()));
    //fvPatchField<scalar>::updateCoeffs();*/
}

Foam::
radialEquilibriumPressureOutletFvPatchScalarField::
radialEquilibriumPressureOutletFvPatchScalarField
(
    const radialEquilibriumPressureOutletFvPatchScalarField& repof
)
:
    fixedValueFvPatchScalarField(repof),
    UName_(repof.UName_),
    rRef_(repof.rRef_),
    pRef_(repof.pRef_),
    n_(repof.n_),
    y_(repof.y_)
{}

Foam::
radialEquilibriumPressureOutletFvPatchScalarField::
radialEquilibriumPressureOutletFvPatchScalarField
(
    const radialEquilibriumPressureOutletFvPatchScalarField& repof,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(repof, iF),
    UName_(repof.UName_),
    rRef_(repof.rRef_),
    pRef_(repof.pRef_),
    n_(repof.n_),
    y_(repof.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void Foam::radialEquilibriumPressureOutletFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    pRef_.autoMap(m);
}


void Foam::totalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const radialEquilibriumPressureOutletFvPatchScalarField& tiptf =
        refCast<const radialEquilibriumPressureOutletFvPatchScalarField>(ptf);

    pRef_.rmap(tiptf.pRef_, addr);
}
*/

void Foam::radialEquilibriumPressureOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Ucar = lookupPatchField<volVectorField, vector>(UName_);
    const fvPatchScalarField& rho = lookupPatchField<volScalarField, scalar>("rho");
	
   /* // get span bounds and check whether rRef is out of range
    const pointField& pp = patch().patch().localPoints();
    scalarField ppRadius = mag(pp - (pp & n_)*n_);
    scalar ppRmin = min(ppRadius);
    scalar ppRmax = max(ppRadius);
    //Info << ppRmin << endl<< ppRmax << endl;
    if( rRef_<ppRmin || rRef_>ppRmax )
    {
    	FatalErrorIn
        (
            "void radialEquilibriumPressureOutletFvPatchScalarField::updateCoeffs()"
        )   << " rRef= " << rRef_ << " is out of  "
            << patch().name() << " range "
            << abort(FatalError);
    }*/
    
    //get centre point of every surface element
    const vectorField& carc = patch().Cf();
	
    //Construct cylidricalCS from origin and 2 axes
    cylindricalCS cyl( "cylindricalCS", point((carc[0]&n_),0,0), n_, y_, false ); 
    
    //location transformed into cylindrical coordinates
    vectorField cylc = cyl.localVector(carc);
    
    //put radius and theta into scalarFields
    scalarField radius = cylc.component(vector::X);
    scalarField theta = cylc.component(vector::Y);
    
    //calculate cylindrical velocity
    fvPatchVectorField Ucyl = Ucar;
    vector z_(0,0,0);
    z_ = n_ ^ y_;
    Ucyl.replace(vector::X, ((Ucar & y_)*cos(theta))+((Ucar & z_)*sin(theta)));//Ur
    Ucyl.replace(vector::Y, (-1*(Ucar & y_)*sin(theta))+((Ucar & z_)*cos(theta)));//Uteta
    Ucyl.replace(vector::Z, Ucar & n_); //Uz
    scalarField Uradius = Ucyl.component(vector::X);
    scalarField Utheta = Ucyl.component(vector::Y);
    
    //read preprocessed point order
    labelList UO(theta.size(), 0);
    IFstream is(db().time().constant()/"pointOrdering");
    is >> UO; 
    
    //start to get unstructure grid dimension
    //general way to get thetaN and radiusN
    //1-initialize a structure
    /*struct location
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
    //3-sort based on radius
	for(int i=1;i<theta.size();i++)//冒泡法排序
	{                     
		for(int j=1;j<(theta.size()-i);j++)
		{
			if(loc[j].radiusX > loc[j+1].radiusX)
			{
				temp=loc[j];
				loc[j]=loc[j+1];
				loc[j+1]=temp;
			}
		}
	}
    //4-calculate thetaN
    
    int thetaN = 0;
    for(int i = 0; i<(theta.size()); i++)
    {
    	if(mag(radius[0]-radius[i])<=1e-6)
    	{
    		++thetaN;
    	}
    	else
    		break;
    }
    
    //simple but noob way to get thetaN and radiusN
    int thetaN = 1;
    for(label i = 0; i<(theta.size()); i++)
    {
    	if((theta[i+1]-theta[i])>0)
    	{
    		++thetaN;
    	}
    	else
    		break;
    }
    int radiusN = theta.size()/thetaN;
    
    //find faces at the same radius
    */
    int thetaN = UO[0];//radiusN was calculated
    int radiusN = theta.size()/thetaN;
    
    //initialize 1D scalarField rhoAvg, UradiusAvg, UthetaAvg, rRE and pRE
    scalarList rhoAvg(radiusN, 0.0);
    scalarList UradiusAvg(radiusN, 0.0);
    scalarList UthetaAvg(radiusN, 0.0);
    scalarList rRE(radiusN, 0.0);
    scalarList pRE(radiusN, 0.0);
    
    //calculate radial-averged rho, Uradius, Utheta and 1D-radius
    for(int i = 0; i<radiusN; i++)
    {
    	for(int j = i*thetaN; j<((i+1)*thetaN); j++)
    	{
    		UradiusAvg[i] += Uradius[UO[j+1]]/thetaN;
    		UthetaAvg[i] += Utheta[UO[j+1]]/thetaN;
    		rhoAvg[i] += rho[UO[j+1]]/thetaN;
    		rRE[i] += radius[UO[j+1]]/thetaN;
    	}
    }
    //Info << UradiusAvg << endl;
    //Info << UthetaAvg << endl;
    //Info << rhoAvg << endl;
    //Info << rRE << endl;
    //calculate 1D pressure distribution
    //step1-find where rRef located
    int rRefLocation = 0;
    if(max(rRE)<rRef_)
    {
    	rRefLocation = radiusN;
    }
    else
    {
    	for(int i = 0; i<radiusN; i++)
    	{
    		if(rRE[i]>=rRef_)
    		{
    			rRefLocation = i;//rRef is between radius[rRefLocation-1] and radius[rRefLocation]
    			break;
    		}
    		else
    			continue;
    	}
    }
    
    //step2-calculate pressure--simplified radial equilibirium equation
    for(int i=rRefLocation; i<radiusN; i++)
    {
    	if(i==rRefLocation)
    	{
    		pRE[i] = rhoAvg[i]*pow(UthetaAvg[i],2)*(rRE[i]-rRef_)/rRE[i]
    			+pRef_;
    	}
    	else
    	{		
    		pRE[i] = rhoAvg[i]*pow(UthetaAvg[i],2)*(rRE[i]-rRE[i-1])/rRE[i]
    			+pRE[i-1];
    	}
    	
    }
    for(int i=(rRefLocation-1); i>=0; i--)
    {
    	if(i==(rRefLocation-1))
    	{
    		pRE[i] = rhoAvg[i]*pow(UthetaAvg[i],2)*(rRE[i]-rRef_)/rRE[i]
    			+pRef_;
    		
    	}
    	else
    	{		
    		pRE[i] = rhoAvg[i]*pow(UthetaAvg[i],2)*(rRE[i]-rRE[i+1])/rRE[i]
    			+pRE[i+1];
    	}
    	
    }
    /*
    //step2-calculate pressure--full radial equilibirium equation
    scalar rhoRef = rhoAvg[rRefLocation]
    			-(rRE[rRefLocation]-rRef_)*(rhoAvg[rRefLocation]-rhoAvg[rRefLocation-1])
    			/(rRE[rRefLocation]-rRE[rRefLocation-1]);
    scalar UradiusRef = UradiusAvg[rRefLocation]
    			-(rRE[rRefLocation]-rRef_)*(UradiusAvg[rRefLocation]-UradiusAvg[rRefLocation-1])
    			/(rRE[rRefLocation]-rRE[rRefLocation-1]);
    
    for(int i=rRefLocation; i<radiusN; i++)
    {
    	if(i==rRefLocation)
    	{
    		pRE[i] = rhoAvg[i]*(pow(UthetaAvg[i],2)-pow(UradiusAvg[i],2))*(rRE[i]-rRef_)/rRE[i]
    			-(rhoAvg[i]*pow(UradiusAvg[i],2)-(rhoRef*pow(UradiusRef,2)))
    			+pRef_;
    	}
    	else
    	{		
    		pRE[i] = rhoAvg[i]*(pow(UthetaAvg[i],2)-pow(UradiusAvg[i],2))*(rRE[i]-rRE[i-1])/rRE[i]
    			-(rhoAvg[i]*pow(UradiusAvg[i],2)-(rhoAvg[i-1]*pow(UradiusAvg[i-1],2)))
    			+pRE[i-1];
    	}
    	
    }
    //Info << "checkPoint3" << endl;
    for(int i=(rRefLocation-1); i>=0; i--)
    {
    	if(i==(rRefLocation-1))
    	{
    		pRE[i] = rhoAvg[i]*(pow(UthetaAvg[i],2)-pow(UradiusAvg[i],2))*(rRE[i]-rRef_)/rRE[i]
    			-(rhoAvg[i]*pow(UradiusAvg[i],2)-(rhoRef*pow(UradiusRef,2)))
    			+pRef_;
    		
    	}
    	else
    	{		
    		pRE[i] = rhoAvg[i]*(pow(UthetaAvg[i],2)-pow(UradiusAvg[i],2))*(rRE[i]-rRE[i+1])/rRE[i]
    			-(rhoAvg[i]*pow(UradiusAvg[i],2)-(rhoAvg[i+1]*pow(UradiusAvg[i+1],2)))
    			+pRE[i+1];
    	}
    	
    }
*/
    scalarField pREreal(theta.size(),pRef_);
    //scalarField pREreal = rho-rho;
    for(int i = 0; i<radiusN; i++)
    {
    	for(int j = i*thetaN; j<((i+1)*thetaN); j++)
    	{
    		pREreal[j] = pRE[i];
    		
    	}
    }
    
    operator==
    (
    	pREreal
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}



void Foam::radialEquilibriumPressureOutletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("UName")
        << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rRef")
        << rRef_ << token::END_STATEMENT << nl;
    os.writeKeyword("pRef")
        << pRef_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

    makePatchTypeField
    (
        fvPatchScalarField,
        radialEquilibriumPressureOutletFvPatchScalarField
    );
}

// ************************************************************************* //
