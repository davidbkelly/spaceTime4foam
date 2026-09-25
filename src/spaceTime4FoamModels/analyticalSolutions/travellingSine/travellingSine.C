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

#include "travellingSine.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace analyticalSolutions
{
    defineTypeNameAndDebug(travellingSine, 0);
    addToRunTimeSelectionTable
    (
        analyticalSolution,
        travellingSine,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolutions::travellingSine::travellingSine
(
    const dictionary& dict,
    const vector& A
)
:
    analyticalSolution(dict, A),
    k_(dict.get<scalar>("k"))
{
    Info<< "    travellingSine: u = sin(2 pi k (x - a t)) with k = " << k_
        << " and a = " << A.x() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolutions::travellingSine::~travellingSine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::analyticalSolutions::travellingSine::value
(
    const point& p
) const
{
    const scalar a = A().x();
    const scalar x = p.x();
    const scalar t = p.y();

    return Foam::sin(constant::mathematical::twoPi*k_*(x - a*t));
}


Foam::vector Foam::analyticalSolutions::travellingSine::gradient
(
    const point& p
) const
{
    const scalar a = A().x();
    const scalar x = p.x();
    const scalar t = p.y();

    const scalar omega = constant::mathematical::twoPi*k_;
    const scalar dudx = omega*Foam::cos(omega*(x - a*t));

    return vector(dudx, -a*dudx, 0);
}


Foam::scalar Foam::analyticalSolutions::travellingSine::source
(
    const point& p
) const
{
    return 0;
}


// ************************************************************************* //
