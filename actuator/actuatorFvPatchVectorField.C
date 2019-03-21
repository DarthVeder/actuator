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
    pressureFile_(),
    ra_(1.0),
    xa_(Zero),
    fn_tot_(0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}

actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    pressureFile_(Function1<scalar>::New("pressureValue", dict)),
    ra_(readScalar(dict.lookup("ra"))),
    xa_(vector(dict.lookup("xa"))),
    fn_tot_(readScalar(dict.lookup("fn_tot")))
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
    pressureFile_(tdpvf.pressureFile_.clone()),
    ra_(tdpvf.ra_),
    xa_(tdpvf.xa_),
    fn_tot_(tdpvf.fn_tot_)
{}


actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const actuatorFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    pressureFile_(tdpvf.pressureFile_.clone()),
    ra_(tdpvf.ra_),
    xa_(tdpvf.xa_),
    fn_tot_(tdpvf.fn_tot_)
{}


actuatorFvPatchVectorField::
actuatorFvPatchVectorField
(
    const actuatorFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    pressureFile_(tdpvf.pressureFile_.clone()),
    ra_(tdpvf.ra_),
    xa_(tdpvf.xa_),
    fn_tot_(tdpvf.fn_tot_)
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

    if (mechanicalProperties.get<bool>("planeStress"))
    {
        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
        threeK = E/(1.0 - nu);
    }
    scalarField twoMuLambda(2*mu + lambda);
    vectorField n(patch().nf());
    scalarField Sf_(patch().magSf());
    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");
	vectorField snGrad = fvPatchField<vector>::snGrad();
    vectorField Cf_ = patch().Cf();

    // Computing total area (m^2) on which the actuator force is used.
    // NOTE To save computation time, maybe could it be passed by user?
    scalar As = 0.0; // m^2
    scalarField dist_Cf_xa(Cf_.size()); // m
    forAll(Sf_, i)
    {
        dist_Cf_xa[i] = 0.0;
  	    for (int j=0; j<3; j++)
   	    {
			dist_Cf_xa[i] += (Cf_[i][j]-xa_[j]) * (Cf_[i][j]-xa_[j]);
	   	}
        dist_Cf_xa[i] = sqrt(dist_Cf_xa[i]);
        if (dist_Cf_xa[i] <= ra_)
        {
            As += Sf_[i];
        }
    }
    if (As == 0.0) As = GREAT;

    const scalar pvalue_ = fn_tot_ * pressureFile_->value(t);
    forAll(Cf_, i)
    {
    	if (dist_Cf_xa[i] <= ra_)
    	{ // NEUMANN, inside actuator disk
			gradient()[i] = (   pvalue_*n[i]/rho[i]
						      + twoMuLambda[i]*snGrad[i] - (n[i] & sigmaD[i])
					        )/twoMuLambda[i];

    	}
		else
		{ // NEUMANN, outside actuator disk
			gradient()[i] = ( twoMuLambda[i]*snGrad[i] - (n[i] & sigmaD[i])
					        )/twoMuLambda[i];

		}
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void actuatorFvPatchVectorField::write(Ostream& os) const
{

    fvPatchVectorField::write(os);
    os.writeKeyword("fn_tot") << fn_tot_ << token::END_STATEMENT << nl;
    pressureFile_->writeData(os);
    os.writeKeyword("ra") << ra_ << token::END_STATEMENT << nl;
    os.writeKeyword("xa") << xa_ << token::END_STATEMENT << nl;
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
