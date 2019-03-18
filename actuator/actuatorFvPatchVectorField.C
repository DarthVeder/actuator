/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "actuatorFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    f0_(1.0),
    ra_(1.0),
    nc_(1.0),
    xa_(Zero),
    fn_max_(0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const actuatorFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    f0_(tdpvf.f0_),
    ra_(tdpvf.ra_),
    nc_(tdpvf.nc_),
    xa_(tdpvf.xa_),
    fn_max_(tdpvf.fn_max_)
{}


actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    f0_(readScalar(dict.lookup("f0"))),
    ra_(readScalar(dict.lookup("ra"))),
    nc_(readScalar(dict.lookup("nc"))),
    xa_(vector(dict.lookup("xa"))),
    fn_max_(readScalar(dict.lookup("fn_max")))
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const actuatorFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    f0_(tdpvf.f0_),
    ra_(tdpvf.ra_),
    nc_(tdpvf.nc_),
    xa_(tdpvf.xa_),
    fn_max_(tdpvf.fn_max_)
{}


actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const actuatorFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    f0_(tdpvf.f0_),
    ra_(tdpvf.ra_),
    nc_(tdpvf.nc_),
    xa_(tdpvf.xa_),
    fn_max_(tdpvf.fn_max_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void actuatorFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
}


void actuatorFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);
}


void actuatorFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& mechanicalProperties =
        db().lookupObject<IOdictionary>("mechanicalProperties");

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    const fvPatchField<scalar>& rhoE =
        patch().lookupPatchField<volScalarField, scalar>("E");

    const fvPatchField<scalar>& nu =
        patch().lookupPatchField<volScalarField, scalar>("nu");

    scalarField E(rhoE/rho);
    scalarField mu(E/(2.0*(1.0 + nu)));
    scalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
    scalarField threeK(E/(1.0 - 2.0*nu));
    scalar t = db().time().value();
    const scalar pi = constant::mathematical::pi;


    if (mechanicalProperties.get<bool>("planeStress"))
    {
        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
        threeK = E/(1.0 - nu);
    }
    scalarField twoMuLambda(2*mu + lambda);
    vectorField n(patch().nf());
    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");
	vectorField snGrad = fvPatchField<vector>::snGrad();

    vectorField Cf_ = patch().Cf();   
    scalar t_burst = nc_/f0_; // s                     
    scalar As = pi*ra_*ra_;
    scalar pn_max = fn_max_/As;

    scalar fvalue_ = 0.0;
    scalar dist_Cf_xa = 0.0;
    forAll(Cf_, i)
    {
        //cout<<Cf_[i][0]<<" "<<Cf_[i][1]<<" "<<Cf_[i][2]<<"\n";
        dist_Cf_xa = 0.0;
  	for (int j=0; j<3; j++) 
   	{
			dist_Cf_xa += (Cf_[i][j]-xa_[j]) * (Cf_[i][j]-xa_[j]);
	   	}
    	dist_Cf_xa = sqrt(dist_Cf_xa);
		// Computing sine-burst, given frequency_ and ncycles_
    	fvalue_ = 0.0;
    	if (t<=t_burst && dist_Cf_xa <= ra_)
    	{ // NEUMANN
            cout<<"..........\n";	
        	fvalue_ = pn_max * sin(2.0*pi*f0_*t) * pow(sin(pi*t/t_burst),2.0);
			gradient()[i] = ( - fvalue_*n[i]/rho[i] 
						      + twoMuLambda[i]*snGrad[i] - (n[i] & sigmaD[i]) 
					        )/twoMuLambda[i];

    	}	
		else
		{ // NEUMANN
			gradient()[i] = ( twoMuLambda[i]*snGrad[i] - (n[i] & sigmaD[i]) 
					        )/twoMuLambda[i];

		}
    }
    
    fixedGradientFvPatchVectorField::updateCoeffs();
}


void actuatorFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("f0") << f0_ << token::END_STATEMENT << nl;
    os.writeKeyword("ra") << ra_ << token::END_STATEMENT << nl;
    os.writeKeyword("nc") << nc_ << token::END_STATEMENT << nl;
    os.writeKeyword("xa") << xa_ << token::END_STATEMENT << nl;
    os.writeKeyword("fn_max") << fn_max_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    actuatorFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
