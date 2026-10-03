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

#include "analyticalSolution.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(analyticalSolution, 0);
    defineRunTimeSelectionTable(analyticalSolution, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::analyticalSolution::analyticalSolution
(
    const dictionary& dict,
    const vector& A
)
:
    dict_(dict),
    A_(A)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::analyticalSolution::~analyticalSolution()
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::vector Foam::analyticalSolution::spaceTimeVelocity
(
    const dictionary& spaceTimeProperties
)
{
    const scalar a = spaceTimeProperties.get<scalar>("a");

    return vector(a, 1, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::analyticalSolution::source(const point& p) const
{
    // div_st(A u) = A & grad(u) for a constant A
    return A_ & gradient(p);
}


Foam::tmp<Foam::scalarField> Foam::analyticalSolution::value
(
    const vectorField& points
) const
{
    tmp<scalarField> tresult(new scalarField(points.size()));
    scalarField& result = tresult.ref();

    forAll(points, pointi)
    {
        result[pointi] = value(points[pointi]);
    }

    return tresult;
}


Foam::tmp<Foam::vectorField> Foam::analyticalSolution::gradient
(
    const vectorField& points
) const
{
    tmp<vectorField> tresult(new vectorField(points.size()));
    vectorField& result = tresult.ref();

    forAll(points, pointi)
    {
        result[pointi] = gradient(points[pointi]);
    }

    return tresult;
}


Foam::tmp<Foam::scalarField> Foam::analyticalSolution::source
(
    const vectorField& points
) const
{
    tmp<scalarField> tresult(new scalarField(points.size()));
    scalarField& result = tresult.ref();

    forAll(points, pointi)
    {
        result[pointi] = source(points[pointi]);
    }

    return tresult;
}


// ************************************************************************* //
