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

#include "mmsQuadratic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace analyticalSolutions
{
    defineTypeNameAndDebug(mmsQuadratic, 0);
    addToRunTimeSelectionTable
    (
        analyticalSolution,
        mmsQuadratic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolutions::mmsQuadratic::mmsQuadratic
(
    const dictionary& dict,
    const vector& A
)
:
    analyticalSolution(dict, A)
{
    Info<< "    mmsQuadratic: u = 3 x^2 + 5 t^2" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolutions::mmsQuadratic::~mmsQuadratic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::analyticalSolutions::mmsQuadratic::value
(
    const point& p
) const
{
    const scalar x = p.x();
    const scalar t = p.y();

    return 3*x*x + 5*t*t;
}


Foam::vector Foam::analyticalSolutions::mmsQuadratic::gradient
(
    const point& p
) const
{
    const scalar x = p.x();
    const scalar t = p.y();

    return vector(6*x, 10*t, 0);
}


// ************************************************************************* //
