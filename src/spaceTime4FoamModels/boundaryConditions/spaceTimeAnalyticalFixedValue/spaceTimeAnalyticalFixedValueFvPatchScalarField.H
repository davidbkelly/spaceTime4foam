/*---------------------------------------------------------------------------*\
License
    This file is part of spaceTime4foam.

    spaceTime4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    spaceTime4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with spaceTime4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField

Description
    Fixed value boundary condition that sets the patch face values to the
    analytical solution evaluated at the patch face centres (x, t).

    The analytical solution and the space-time velocity are taken from
    constant/spaceTimeProperties (keyword a and the analyticalSolution
    sub-dictionary). The dictionary registered on the Time database by the
    spaceTime4Foam solver is used if it exists; otherwise the file is read,
    so the condition also works in utilities such as foamToVTK and
    spaceTimeErrors.

    Any value entry is ignored: the value is always recomputed.

Usage
    Example of the boundary condition specification:
    \verbatim
    tStart
    {
        type            spaceTimeAnalyticalFixedValue;
    }
    \endverbatim

SourceFiles
    spaceTimeAnalyticalFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef spaceTimeAnalyticalFixedValueFvPatchScalarField_H
#define spaceTimeAnalyticalFixedValueFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class spaceTimeAnalyticalFixedValueFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class spaceTimeAnalyticalFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Analytical solution, constructed on first use
        mutable autoPtr<analyticalSolution> analyticalSolutionPtr_;


    // Private Member Functions

        //- Return the analytical solution, constructing it if needed
        const analyticalSolution& solution() const;

        //- Return the analytical solution at the patch face centres
        tmp<scalarField> analyticalValues() const;


public:

    //- Runtime type information
    TypeName("spaceTimeAnalyticalFixedValue");


    // Constructors

        //- Construct from patch and internal field
        spaceTimeAnalyticalFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        spaceTimeAnalyticalFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  spaceTimeAnalyticalFixedValueFvPatchScalarField onto a new patch
        spaceTimeAnalyticalFixedValueFvPatchScalarField
        (
            const spaceTimeAnalyticalFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        spaceTimeAnalyticalFixedValueFvPatchScalarField
        (
            const spaceTimeAnalyticalFixedValueFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        spaceTimeAnalyticalFixedValueFvPatchScalarField
        (
            const spaceTimeAnalyticalFixedValueFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new spaceTimeAnalyticalFixedValueFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new spaceTimeAnalyticalFixedValueFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
