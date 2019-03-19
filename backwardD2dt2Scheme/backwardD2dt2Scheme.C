/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "backwardD2dt2Scheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT_() const
{
    return mesh().time().deltaTValue();
}


template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0Value();
}


template<class Type>
template<class GeoField>
scalar backwardD2dt2Scheme<Type>::deltaT0_(const GeoField& vf) const
{
    if (mesh().time().timeIndex() < 3)
    {
        return GREAT;
    }
    else
    {
        return deltaT_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
backwardD2dt2Scheme<Type>::fvcD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT2 = 1.0/sqr(deltaT_());

    IOobject d2dt2IOobject
    (
        "d2dt2("+vf.name()+')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar coefft    = 1 + deltaT_()/deltaT0_();
    scalar coefft0   = 2 + 3*deltaT_()/deltaT0_();
    scalar coefft00  = 1 + 3*deltaT_()/deltaT0_();
    scalar coefft000 = 1 * deltaT_()/deltaT0_();

    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            d2dt2IOobject,
            rDeltaT2*
            (
                coefft*vf
              - coefft0*vf.oldTime()
              + coefft00*vf.oldTime().oldTime()
              - coefft000*vf.oldTime().oldTime().oldTime()
            )
       )
   );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
backwardD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT2 = 1.0/sqr(deltaT_());

    IOobject d2dt2IOobject
    (
        "d2dt2("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar coefft    = 1 + deltaT_()/deltaT0_();
    scalar coefft0   = 2 + 3*deltaT_()/deltaT0_();
    scalar coefft00  = 1 + 3*deltaT_()/deltaT0_();
    scalar coefft000 = 1 * deltaT_()/deltaT0_();

    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            d2dt2IOobject,
            rDeltaT2*
            (
                coefft*rho*vf
              - coefft0*rho.oldTime()*vf.oldTime()
              + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
              - coefft000*rho.oldTime().oldTime().oldTime()*vf.oldTime().oldTime().oldTime()
            )
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>>
backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT2 = 1.0/sqr(deltaT_());

    scalar coefft    = 1 + deltaT_()/deltaT0_();
    scalar coefft0   = 2 + 3*deltaT_()/deltaT0_();
    scalar coefft00  = 1 + 3*deltaT_()/deltaT0_();
    scalar coefft000 = 1 * deltaT_()/deltaT0_();

    fvm.diag() = (coefft*rDeltaT2)*mesh().V();

    fvm.source() = rDeltaT2*mesh().V()*
    (
        coefft0*vf.oldTime().primitiveField()
      - coefft00*vf.oldTime().oldTime().primitiveField()
      + coefft000*vf.oldTime().oldTime().oldTime().primitiveField()
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT2 = 1.0/sqr(deltaT_());

    scalar coefft    = 1 + deltaT_()/deltaT0_();
    scalar coefft0   = 2 + 3*deltaT_()/deltaT0_();
    scalar coefft00  = 1 + 3*deltaT_()/deltaT0_();
    scalar coefft000 = 1 * deltaT_()/deltaT0_();

    fvm.diag() = (coefft*rDeltaT2*rho.value())*mesh().V();

    fvm.source() = rDeltaT2*mesh().V()*rho.value()*
    (
        coefft0*vf.oldTime().primitiveField()
      - coefft00*vf.oldTime().oldTime().primitiveField()
      + coefft000*vf.oldTime().oldTime().oldTime().primitiveField()
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT2 = 1.0/sqr(deltaT_());
    
    scalar coefft    = 1 + deltaT_()/deltaT0_();
    scalar coefft0   = 2 + 3*deltaT_()/deltaT0_();
    scalar coefft00  = 1 + 3*deltaT_()/deltaT0_();
    scalar coefft000 = 1 * deltaT_()/deltaT0_();

    fvm.diag() = (coefft*rDeltaT2)*rho.primitiveField()*mesh().V();

    fvm.source() = rDeltaT2*mesh().V()*
    (
        coefft0*rho.oldTime().primitiveField()
       *vf.oldTime().primitiveField()
      - coefft00*rho.oldTime().oldTime().primitiveField()
       *vf.oldTime().oldTime().primitiveField()
      + coefft000*rho.oldTime().oldTime().oldTime().primitiveField()
       *vf.oldTime().oldTime().oldTime().primitiveField()
    );

    return tfvm;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
