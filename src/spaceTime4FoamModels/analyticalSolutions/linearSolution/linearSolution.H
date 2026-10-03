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
    Foam::analyticalSolutions::linearSolution

Description
    Linear solution u(x, t) = c0 + c . (x, t), with gradient c = (cx, ct, 0)
    and source S = A & c. Used for the linear-exactness test of the
    vertex-centred residual (CLAUDE.md section 6, test 5).

    c is given as (cx ct) or (cx ct 0); a non-zero third component is an
    error. hasSource() is true unless A & c is exactly zero.

    The class is not called linear, to keep it apart from the Foam::linear
    interpolation scheme.

    Example, in constant/spaceTimeProperties:
    \verbatim
    analyticalSolution
    {
        type            linear;
        c0              0.3;
        c               (1.7 -2.3);
    }
    \endverbatim

SourceFiles
    linearSolution.C

\*---------------------------------------------------------------------------*/

#ifndef linearSolution_H
#define linearSolution_H

#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace analyticalSolutions
{

/*---------------------------------------------------------------------------*\
                       Class linearSolution Declaration
\*---------------------------------------------------------------------------*/

class linearSolution
:
    public analyticalSolution
{
    // Private data

        //- The value at the origin, c0
        const scalar c0_;

        //- The gradient c = (cx, ct, 0)
        vector c_;


    // Private Member Functions

        //- No copy construct
        linearSolution(const linearSolution&) = delete;

        //- No copy assignment
        void operator=(const linearSolution&) = delete;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from the coefficients dictionary and space-time velocity
        linearSolution(const dictionary& dict, const vector& A);


    //- Destructor
    virtual ~linearSolution();


    // Member Functions

        // Evaluation at one space-time point

            //- Return u = c0 + cx x + ct t
            virtual scalar value(const point& p) const;

            //- Return the gradient c
            virtual vector gradient(const point& p) const;

            //- True unless A & c is zero
            virtual bool hasSource() const
            {
                return (A() & c_) != 0;
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
