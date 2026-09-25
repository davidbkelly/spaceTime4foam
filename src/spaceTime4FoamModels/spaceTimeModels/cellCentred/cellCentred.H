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
    Foam::spaceTimeModels::cellCentred

Description
    Cell-centred finite volume (CCFV) space-time discretisation of the
    steady space-time advection problem

        div_st(A u) = 0,    A = (a, 1, 0)

    solved as fvm::div(phiST, u) == 0 with phiST = fvc::flux(Ust) and
    Ust = A uniform. The face interpolation is chosen in system/fvSchemes
    by div(phiST,u), e.g. Gauss upwind (CC-1) or Gauss linearUpwind grad(u)
    (CC-2). The matrix always holds the upwind part; for a scheme with an
    explicit correction (corrected() true, e.g. linearUpwind) OpenFOAM adds
    the correction from the current u to the source, so each iteration is
    one step of deferred correction.

    One linear solve per iteration (use a nonsymmetric solver, e.g.
    PBiCGStab with DILU). Converged when the initial residual of the solve
    falls below convergenceTolerance.

SourceFiles
    cellCentred.C

\*---------------------------------------------------------------------------*/

#ifndef cellCentred_H
#define cellCentred_H

#include "spaceTimeModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace spaceTimeModels
{

/*---------------------------------------------------------------------------*\
                         Class cellCentred Declaration
\*---------------------------------------------------------------------------*/

class cellCentred
:
    public spaceTimeModel
{
    // Private data

        //- Solution field u, read from 0/u
        volScalarField u_;

        //- Space-time velocity field, uniform A
        volVectorField Ust_;

        //- Space-time face flux, fvc::flux(Ust)
        surfaceScalarField phiST_;

        //- Initial residual of the latest iteration
        scalar residual_;


    // Private Member Functions

        //- No copy construct
        cellCentred(const cellCentred&) = delete;

        //- No copy assignment
        void operator=(const cellCentred&) = delete;


public:

    //- Runtime type information
    TypeName("cellCentred");


    // Constructors

        //- Construct from Time
        cellCentred(Time& runTime);


    //- Destructor
    virtual ~cellCentred();


    // Member Functions

        // Access

            //- Return the solution field
            const volScalarField& u() const
            {
                return u_;
            }


        // Evolution

            //- Perform one iteration; return true once converged
            virtual bool evolve();

            //- Return the initial residual of the latest iteration
            virtual scalar residual() const
            {
                return residual_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace spaceTimeModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
