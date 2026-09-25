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

#include "spaceTimeModel.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(spaceTimeModel, 0);
    defineRunTimeSelectionTable(spaceTimeModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaceTimeModel::spaceTimeModel
(
    const word& type,
    Time& runTime
)
:
    runTime_(runTime),
    spaceTimeProperties_
    (
        IOobject
        (
            "spaceTimeProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    ),
    A_(analyticalSolution::spaceTimeVelocity(spaceTimeProperties_)),
    convergenceTolerance_
    (
        spaceTimeProperties_.get<scalar>("convergenceTolerance")
    ),
    analyticalSolutionPtr_
    (
        analyticalSolution::New
        (
            spaceTimeProperties_.subDict("analyticalSolution"),
            A_
        )
    )
{
    Info<< "Space-time velocity A = " << A_ << nl
        << "Convergence tolerance on the initial residual = "
        << convergenceTolerance_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spaceTimeModel::~spaceTimeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::spaceTimeModel::writeFields()
{
    runTime_.writeNow();
}


// ************************************************************************* //
