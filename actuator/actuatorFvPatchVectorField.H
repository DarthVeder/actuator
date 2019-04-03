/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::actuatorFvPatchVectorField

Description
    Actuator force acting normal (p*n) to the patch with p read from a file. 
    It is modelled as a circular disk (pi*ra_^2) on the patch with position xa_
    and radius ra_. 

    Boundary format for file:
    patch_name
    {
        type            actuator;
        ra              0.019478; // disk radius, (m)
        xa              ( 0 0.0031 0.0015 ); // disk center, (m)
        pressureValue
        {
            type csvFile; // file ha n rows with  "time value"
            nHeaderLine 0;
            refColumn 0;
            componentColumns (1);
            separator " ";
            mergeSeparators yes;
            file "<constant>/force.txt";

        }
        value           uniform (0 0 0);
    }

SourceFiles
    actuatorFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorFvPatchVectorField_H
#define actuatorFvPatchVectorField_H
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class actuatorFvPatch Declaration
\*---------------------------------------------------------------------------*/

class actuatorFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data
        //- reader for pressure data in external file
        autoPtr<Function1<scalar>> pressureFile_;
        //- Radius of actuator (m)
        scalar ra_;
        //- Position of the actuator center
        vector xa_;
        //- Maximum pressure (Pa)
        scalar pn_max_;

public:

    //- Runtime type information
    TypeName("actuator");


    // Constructors

        //- Construct from patch and internal field
        actuatorFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        actuatorFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  actuatorFvPatchVectorField onto a new patch
        actuatorFvPatchVectorField
        (
            const actuatorFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        actuatorFvPatchVectorField
        (
            const actuatorFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new actuatorFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        actuatorFvPatchVectorField
        (
            const actuatorFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new actuatorFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            // - Return actuator position (m)
            virtual const vector& xa() const
            {
                return xa_;
            }

            virtual vector& xa()
            {
                return xa_;
            }

            // - Return actuator radius (m)
            virtual const scalar& ra() const
            {
                return ra_;
            }

            virtual  scalar& ra()
            {
                return ra_;
            }
            // - Return maximum normal force (N)
            virtual const scalar& fnMax() const
            {
                return pn_max_;
            }

            virtual  scalar& fnMax()
            {
                return pn_max_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


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
