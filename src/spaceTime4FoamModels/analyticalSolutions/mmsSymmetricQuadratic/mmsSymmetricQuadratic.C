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

#include "mmsSymmetricQuadratic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace analyticalSolutions
{
    defineTypeNameAndDebug(mmsSymmetricQuadratic, 0);
    addToRunTimeSelectionTable
    (
        analyticalSolution,
        mmsSymmetricQuadratic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolutions::mmsSymmetricQuadratic::mmsSymmetricQuadratic
(
    const dictionary& dict,
    const vector& A
)
:
    analyticalSolution(dict, A)
{
    Info<< "    mmsSymmetricQuadratic: u = x^2 + x t + t^2" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolutions::mmsSymmetricQuadratic::~mmsSymmetricQuadratic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::analyticalSolutions::mmsSymmetricQuadratic::value
(
    const point& p
) const
{
    const scalar x = p.x();
    const scalar t = p.y();

    return x*x + x*t + t*t;
}


Foam::vector Foam::analyticalSolutions::mmsSymmetricQuadratic::gradient
(
    const point& p
) const
{
    const scalar x = p.x();
    const scalar t = p.y();

    return vector(2*x + t, x + 2*t, 0);
}


// ************************************************************************* //
