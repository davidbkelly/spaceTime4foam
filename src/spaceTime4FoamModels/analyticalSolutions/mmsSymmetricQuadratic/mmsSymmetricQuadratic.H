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
    Foam::analyticalSolutions::mmsSymmetricQuadratic

Description
    Manufactured solution u1 of Tufillaro, Williams and Nishikawa (2026),
    Eq. 122, with y = t:

        u(x, t) = x^2 + x t + t^2,   grad(u) = (2x + t, x + 2t, 0)

    The source is S = A & grad(u) (base class); with A = (1, 1) it is
    f1 = 3 (x + t), their Eq. 125.

    Example, in constant/spaceTimeProperties:
    \verbatim
    analyticalSolution
    {
        type            mmsSymmetricQuadratic;
    }
    \endverbatim

SourceFiles
    mmsSymmetricQuadratic.C

\*---------------------------------------------------------------------------*/

#ifndef mmsSymmetricQuadratic_H
#define mmsSymmetricQuadratic_H

#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace analyticalSolutions
{

/*---------------------------------------------------------------------------*\
                   Class mmsSymmetricQuadratic Declaration
\*---------------------------------------------------------------------------*/

class mmsSymmetricQuadratic
:
    public analyticalSolution
{
    // Private Member Functions

        //- No copy construct
        mmsSymmetricQuadratic(const mmsSymmetricQuadratic&) = delete;

        //- No copy assignment
        void operator=(const mmsSymmetricQuadratic&) = delete;


public:

    //- Runtime type information
    TypeName("mmsSymmetricQuadratic");


    // Constructors

        //- Construct from the coefficients dictionary and space-time velocity
        mmsSymmetricQuadratic(const dictionary& dict, const vector& A);


    //- Destructor
    virtual ~mmsSymmetricQuadratic();


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
