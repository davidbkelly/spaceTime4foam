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

    Corner fallback: in a cell with two or more faces with this condition
    (counted over all patches; on the benchmark meshes only the top-right
    outflow corner cell, with one tEnd and one xRight face) the
    extrapolation is unstable or non-unique. In such a cell the gradient
    component v with d_D . v = 0, where D is the remaining (third) face and
    d_D = C_D - C_P, is fed back into itself through the extrapolated face
    values:
    - Gauss linear: the gain is exactly 1 (sum_f S_f d_f^T = V I), a
      neutral mode. The discrete solution is not unique: the corner face
      values (and for a != 1 the corner cell value and cells downstream of
      it) depend on the starting field.
    - leastSquares: the gain is 1.5 on the left-diagonal meshes, so the
      mode grows from round-off. For a = 1 the corner cell balance does
      not see this direction, so the residual converges while the corner
      face values diverge.
    With cornerFallback on (the default), the faces of such cells use the
    zeroGradient value u_f = u_P (gradient 0). The number of faces using
    the fallback is printed once per patch, on the first update. Cost: for
    a != 1 the corner cell is O(h) wrong and then dominates Linf; for
    a = 1 the corner cell value is unaffected (its two outflow corrections
    cancel). Upgrade path: take the corner cell's gradient from its upwind
    neighbour.

    With cornerFallback off, a cell with two or more such faces gives a
    warning (neutral mode, non-unique solution) if the gradient scheme is
    Gauss-type, and is a fatal error otherwise. Gauss-type (an allow-list)
    means that the first word of the gradSchemes entry is Gauss or
    iterativeGauss, or that the first word is one of the limited schemes
    cellLimited, cellLimited<limiter>, cellMDLimited, faceLimited or
    faceMDLimited and the second word is Gauss or iterativeGauss (e.g.
    cellLimited Gauss linear 1). Any other scheme, including leastSquares
    (gain 1.5) and schemes from user libraries, is stopped.

Usage
    \table
        Property       | Description                   | Required | Default
        gradSchemeName | gradSchemes entry for grad(u) | no       | grad(u)
        cornerFallback | u_f = u_P in multi-EX cells   | no       | on
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
#include "Switch.H"

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

        //- Use the zeroGradient value in cells with two or more faces with
        //  this condition
        Switch cornerFallback_;

        //- Faces of this patch whose cell has two or more faces with this
        //  condition, found on the first update
        mutable autoPtr<boolList> multiFacePtr_;


    // Private Member Functions

        //- Return the default gradSchemeName: grad(<field name>)
        word defaultGradSchemeName() const;

        //- Return true if the gradSchemes entry is a Gauss-type scheme
        //  (see the allow-list in the class description)
        bool gaussTypeGradient() const;

        //- Return the faces whose cell has two or more faces with this
        //  condition (over all patches of the field), finding them, and
        //  checking and reporting the corner fallback, on the first call
        const boolList& multiFaces() const;


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

            //- Return true if the corner fallback is on
            bool cornerFallback() const
            {
                return cornerFallback_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchScalarField&, const labelList&);


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
