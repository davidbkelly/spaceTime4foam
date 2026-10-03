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

#include "constantSolution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace analyticalSolutions
{
    defineTypeNameAndDebug(constantSolution, 0);
    addToRunTimeSelectionTable
    (
        analyticalSolution,
        constantSolution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolutions::constantSolution::constantSolution
(
    const dictionary& dict,
    const vector& A
)
:
    analyticalSolution(dict, A),
    c0_(dict.get<scalar>("c0"))
{
    Info<< "    constant: u = c0 with c0 = " << c0_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolutions::constantSolution::~constantSolution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::analyticalSolutions::constantSolution::value
(
    const point& p
) const
{
    return c0_;
}


Foam::vector Foam::analyticalSolutions::constantSolution::gradient
(
    const point& p
) const
{
    return Zero;
}


// ************************************************************************* //
