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

#include "mmsExponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace analyticalSolutions
{
    defineTypeNameAndDebug(mmsExponential, 0);
    addToRunTimeSelectionTable
    (
        analyticalSolution,
        mmsExponential,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolutions::mmsExponential::mmsExponential
(
    const dictionary& dict,
    const vector& A
)
:
    analyticalSolution(dict, A),
    k_(dict.getOrDefault<scalar>("k", 0.1))
{
    Info<< "    mmsExponential: u = exp(k (x + t)) with k = " << k_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolutions::mmsExponential::~mmsExponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::analyticalSolutions::mmsExponential::value
(
    const point& p
) const
{
    return Foam::exp(k_*(p.x() + p.y()));
}


Foam::vector Foam::analyticalSolutions::mmsExponential::gradient
(
    const point& p
) const
{
    const scalar dudx = k_*Foam::exp(k_*(p.x() + p.y()));

    return vector(dudx, dudx, 0);
}


// ************************************************************************* //
