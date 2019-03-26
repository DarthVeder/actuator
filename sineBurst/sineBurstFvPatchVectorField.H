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
    Foam::sineBurstFvPatchVectorField

Description
    Actuator force acting normal (p*n) to the patch and defined  as a pressure 
    sine burst.
    It is modelled as a circular disk (pi*ra_^2) on the patch with position xa_
    and radius ra_. The sine burst requires f0 and nc. The actual area is computed
    at run time, to allow for smaller sections of the actuator being present on the 
    patch. 

    Boundary format for sine burst:
    
    patch_name
    {
        type            sineBurst;
        f0              20000; // (Hz)
        nc              4; // (1)
        fn_tot          -1; // (N), total force. (+) is parallel with the normal() 
        ra              0.019478; // disk radius, (m)
        xa              ( 0 0.0031 0.0015 ); // disk center, (m)
        value           uniform (0 0 0);
    }
 
SourceFiles
    sineBurstFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef sineBurstFvPatchVectorField_H
#define sineBurstFvPatchVectorField_H
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sineBurstFvPatch Declaration
\*---------------------------------------------------------------------------*/

class sineBurstFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data
        //- Frequency of sine burst (Hz)
        scalar f0_;
        //- Radius of actuator (m)
        scalar ra_;
        //- Number of cycles (1)
        scalar nc_;
        //- Position of the actuator center
        vector xa_;
        //- Maximum force (N)
        scalar fn_tot_;

public:

    //- Runtime type information
    TypeName("sineBurst");


    // Constructors

        //- Construct from patch and internal field
        sineBurstFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sineBurstFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  sineBurstFvPatchVectorField onto a new patch
        sineBurstFvPatchVectorField
        (
            const sineBurstFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sineBurstFvPatchVectorField
        (
            const sineBurstFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new sineBurstFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sineBurstFvPatchVectorField
        (
            const sineBurstFvPatchVectorField&,
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
                new sineBurstFvPatchVectorField(*this, iF)
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

            // - Return sine burst frequency (Hz)
            virtual const scalar& f0() const
            {
                return f0_;
            }

            virtual  scalar& f0()
            {
                return f0_;
            }
            // - Return sine burst number of cycles (1)
            virtual const scalar& nc() const
            {
                return nc_;
            }

            virtual  scalar& nc()
            {
                return nc_;
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
                return fn_tot_;
            }

            virtual  scalar& fnMax()
            {
                return fn_tot_;
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
