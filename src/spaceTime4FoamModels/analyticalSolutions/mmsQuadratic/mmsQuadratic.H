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
    Foam::analyticalSolutions::mmsQuadratic

Description
    Manufactured solution u2 of Tufillaro, Williams and Nishikawa (2026),
    Eq. 123, with y = t:

        u(x, t) = 3 x^2 + 5 t^2,   grad(u) = (6x, 10t, 0)

    The source is S = A & grad(u) (base class); with A = (1, 1) it is
    f2 = 6x + 10t, their Eq. 126.

    Example, in constant/spaceTimeProperties:
    \verbatim
    analyticalSolution
    {
        type            mmsQuadratic;
    }
    \endverbatim

SourceFiles
    mmsQuadratic.C

\*---------------------------------------------------------------------------*/

#ifndef mmsQuadratic_H
#define mmsQuadratic_H

#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace analyticalSolutions
{

/*---------------------------------------------------------------------------*\
                        Class mmsQuadratic Declaration
\*---------------------------------------------------------------------------*/

class mmsQuadratic
:
    public analyticalSolution
{
    // Private Member Functions

        //- No copy construct
        mmsQuadratic(const mmsQuadratic&) = delete;

        //- No copy assignment
        void operator=(const mmsQuadratic&) = delete;


public:

    //- Runtime type information
    TypeName("mmsQuadratic");


    // Constructors

        //- Construct from the coefficients dictionary and space-time velocity
        mmsQuadratic(const dictionary& dict, const vector& A);


    //- Destructor
    virtual ~mmsQuadratic();


    // Member Functions

        // Evaluation at one space-time point

            //- Return u at the point p = (x, t, z)
            virtual scalar value(const point& p) const;

            //- Return the space-time gradient (du/dx, du/dt, 0) at p
            virtual vector gradient(const point& p) const;

            //- The source A & grad(u) is not identically zero
            virtual bool hasSource() const
            {
                return true;
            }


        // Evaluation at a list of space-time points

            using analyticalSolution::value;
            using analyticalSolution::gradient;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace analyticalSolutions

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
