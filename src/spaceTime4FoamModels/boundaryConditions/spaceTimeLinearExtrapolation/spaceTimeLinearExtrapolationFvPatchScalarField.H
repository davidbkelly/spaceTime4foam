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
    Foam::spaceTimeLinearExtrapolationFvPatchScalarField

Description
    Outflow boundary condition for space-time CCFV that sets the face value
    by linear extrapolation from the owner cell P with the cell gradient,
    using the full offset vector from the cell centre to the face centre:

        u_f = u_P + grad(u)_P & (C_f - C_P)

    With zeroGradient, u_f = u_P, which is O(h) wrong on outflow faces
    (OpenFOAM's linearUpwind adds no correction on non-coupled patches).

    Implementation: a fixedGradient condition whose gradient is

        g = deltaCoeffs*(grad(u)_P & (C_f - C_P))

    fixedGradient evaluates the face value as u_P + g/deltaCoeffs (in the
    convection matrix: valueInternalCoeffs = 1 and valueBoundaryCoeffs =
    g/deltaCoeffs), so deltaCoeffs cancels and the face value is exactly
    the expression above, including the tangential part of C_f - C_P, for
    any positive deltaCoeffs. (Setting g = grad(u)_P & nf instead would
    keep only the patch-normal part of C_f - C_P.)

    Why fixedGradient: u_P stays implicit in the matrix, as with
    zeroGradient, so the owner's diagonal coefficient is unchanged and only
    the correction grad(u)_P & (C_f - C_P) is explicit. This matches the
    deferred-correction outer iteration of the cellCentred model.

    When grad(u) is updated: the gradient is recomputed in updateCoeffs()
    with fvc::grad(u, gradSchemeName), i.e. with the gradSchemes entry of
    that name in system/fvSchemes, from the current field (internal and
    boundary values). In the spaceTime4Foam iteration updateCoeffs() is
    called once per iteration, when the matrix of div(phiST,u) is built
    from the result of the previous iteration; the linearUpwind correction
    is evaluated from the same field. After the linear solve,
    correctBoundaryConditions() re-evaluates the face value with the new
    u_P but the same g, i.e. with the gradient lagged by one iteration. At
    convergence the lag is of the order of the change of u in one
    iteration.

    Each patch with this condition computes the gradient in its own
    updateCoeffs(). To compute it only once per iteration, cache it in
    system/fvSolution (cache { grad(u); }): all the patches then reuse the
    cached field while u is unchanged.

    The written gradient entry reproduces the face value when the field is
    read again, e.g. by spaceTimeErrors. If the gradient entry is missing
    (e.g. in 0/u), the initial face value is the cell value, as for
    zeroGradient.

    Known limitation: in a cell with two faces with this condition (the
    top-right corner cell of the benchmark meshes) the cell balance does
    not see one direction of grad(u)_P, while the face values feed back
    into grad(u)_P. On the left-diagonal meshes with a = 1 this feedback
    has a gain of 1.5 per iteration with leastSquares, so the corner face
    values diverge while the residual converges; with Gauss linear the
    gain is 1 and that component of the corner gradient keeps its initial
    value (zero). Do not use it with leastSquares until this is resolved.

Usage
    \table
        Property       | Description                   | Required | Default
        gradSchemeName | gradSchemes entry for grad(u) | no       | grad(u)
    \endtable

    The default gradSchemeName is grad(<field name>).

    Example of the boundary condition specification:
    \verbatim
    tEnd
    {
        type            spaceTimeLinearExtrapolation;
    }
    \endverbatim

SourceFiles
    spaceTimeLinearExtrapolationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef spaceTimeLinearExtrapolationFvPatchScalarField_H
#define spaceTimeLinearExtrapolationFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class spaceTimeLinearExtrapolationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class spaceTimeLinearExtrapolationFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the gradSchemes entry used for grad(u)
        word gradSchemeName_;


    // Private Member Functions

        //- Return the default gradSchemeName: grad(<field name>)
        word defaultGradSchemeName() const;


public:

    //- Runtime type information
    TypeName("spaceTimeLinearExtrapolation");


    // Constructors

        //- Construct from patch and internal field
        spaceTimeLinearExtrapolationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        spaceTimeLinearExtrapolationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  spaceTimeLinearExtrapolationFvPatchScalarField onto a new patch
        spaceTimeLinearExtrapolationFvPatchScalarField
        (
            const spaceTimeLinearExtrapolationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        spaceTimeLinearExtrapolationFvPatchScalarField
        (
            const spaceTimeLinearExtrapolationFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        spaceTimeLinearExtrapolationFvPatchScalarField
        (
            const spaceTimeLinearExtrapolationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new spaceTimeLinearExtrapolationFvPatchScalarField(*this)
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
                new spaceTimeLinearExtrapolationFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the name of the gradSchemes entry used for grad(u)
            const word& gradSchemeName() const
            {
                return gradSchemeName_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field:
            //  set the gradient from the current grad(u)
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
