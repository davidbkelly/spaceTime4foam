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

#include "linearSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace analyticalSolutions
{
    defineTypeNameAndDebug(linearSolution, 0);
    addToRunTimeSelectionTable
    (
        analyticalSolution,
        linearSolution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolutions::linearSolution::linearSolution
(
    const dictionary& dict,
    const vector& A
)
:
    analyticalSolution(dict, A),
    c0_(dict.get<scalar>("c0")),
    c_(Zero)
{
    // c is (cx ct) or (cx ct 0)
    const scalarList c(dict.get<scalarList>("c"));

    if (c.size() != 2 && c.size() != 3)
    {
        FatalIOErrorInFunction(dict)
            << "c must have 2 or 3 components, not " << c
            << exit(FatalIOError);
    }

    if (c.size() == 3 && c[2] != 0)
    {
        FatalIOErrorInFunction(dict)
            << "The third (z) component of c must be zero: c = " << c
            << exit(FatalIOError);
    }

    c_ = vector(c[0], c[1], 0);

    Info<< "    linear: u = c0 + c . (x, t) with c0 = " << c0_ << " and c = "
        << c_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolutions::linearSolution::~linearSolution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::analyticalSolutions::linearSolution::value
(
    const point& p
) const
{
    return c0_ + c_.x()*p.x() + c_.y()*p.y();
}


Foam::vector Foam::analyticalSolutions::linearSolution::gradient
(
    const point& p
) const
{
    return c_;
}


// ************************************************************************* //
