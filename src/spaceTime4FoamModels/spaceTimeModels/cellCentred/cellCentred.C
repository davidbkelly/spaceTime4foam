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

\*---------------------------------------------------------------------------*/

#include "cellCentred.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace spaceTimeModels
{
    defineTypeNameAndDebug(cellCentred, 0);
    addToRunTimeSelectionTable(spaceTimeModel, cellCentred, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaceTimeModels::cellCentred::cellCentred(Time& runTime)
:
    spaceTimeModel(typeName, runTime),
    u_
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    Ust_
    (
        IOobject
        (
            "Ust",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("Ust", dimless, A())
    ),
    phiST_
    (
        IOobject
        (
            "phiST",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(Ust_)
    ),
    residual_(GREAT)
{
    // No source term in the CCFV discretisation yet: stop rather than
    // solve the wrong problem
    if (analytical().hasSource())
    {
        FatalErrorInFunction
            << "The analyticalSolution "
            << analytical().type() << " has a non-zero source term, but"
            << " the source term is not implemented in the cellCentred"
            << " spaceTimeModel yet" << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spaceTimeModels::cellCentred::~cellCentred()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::spaceTimeModels::cellCentred::evolve()
{
    // div_st(A u) = 0. The explicit correction of a corrected scheme such as
    // linearUpwind is evaluated from the current u (deferred correction).
    fvScalarMatrix uEqn
    (
        fvm::div(phiST_, u_)
    );

    const SolverPerformance<scalar> solverPerf = uEqn.solve();

    residual_ = solverPerf.initialResidual();

    return residual_ < convergenceTolerance();
}


// ************************************************************************* //
