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
    Foam::analyticalSolutions::travellingSine

Description
    Travelling sine wave, the exact solution of linear advection
    du/dt + a du/dx = 0 (method of characteristics):

        u(x, t) = sin(2 pi k (x - a t))

    with a = A.x(), x = p.x() and t = p.y(). The space-time gradient is

        (du/dx, du/dt, 0) = 2 pi k cos(2 pi k (x - a t)) (1, -a, 0)

    and the source div_st(A u) = du/dt + a du/dx is zero.

    Example, in constant/spaceTimeProperties:
    \verbatim
    analyticalSolution
    {
        type            travellingSine;
        k               1;      // wave number: periods per unit length
    }
    \endverbatim

SourceFiles
    travellingSine.C

\*---------------------------------------------------------------------------*/

#ifndef travellingSine_H
#define travellingSine_H

#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace analyticalSolutions
{

/*---------------------------------------------------------------------------*\
                       Class travellingSine Declaration
\*---------------------------------------------------------------------------*/

class travellingSine
:
    public analyticalSolution
{
    // Private data

        //- Wave number k
        const scalar k_;


    // Private Member Functions

        //- No copy construct
        travellingSine(const travellingSine&) = delete;

        //- No copy assignment
        void operator=(const travellingSine&) = delete;


public:

    //- Runtime type information
    TypeName("travellingSine");


    // Constructors

        //- Construct from the coefficients dictionary and space-time velocity
        travellingSine(const dictionary& dict, const vector& A);


    //- Destructor
    virtual ~travellingSine();


    // Member Functions

        // Evaluation at one space-time point

            //- Return u at the point p = (x, t, z)
            virtual scalar value(const point& p) const;

            //- Return the space-time gradient (du/dx, du/dt, 0) at p
            virtual vector gradient(const point& p) const;

            //- Return the source div_st(A u), which is zero
            virtual scalar source(const point& p) const;


        // Evaluation at a list of space-time points

            using analyticalSolution::value;
            using analyticalSolution::gradient;
            using analyticalSolution::source;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace analyticalSolutions

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
