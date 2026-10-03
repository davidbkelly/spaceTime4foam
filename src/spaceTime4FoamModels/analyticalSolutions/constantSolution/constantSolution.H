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
    Foam::analyticalSolutions::constantSolution

Description
    Constant solution u(x, t) = c0, with zero gradient and zero source. Used
    for the freestream-preservation test of the vertex-centred residual
    (CLAUDE.md section 6, test 4).

    The class is not called constant, which would hide the namespace
    Foam::constant inside Foam::analyticalSolutions.

    Example, in constant/spaceTimeProperties:
    \verbatim
    analyticalSolution
    {
        type            constant;
        c0              0.7;
    }
    \endverbatim

SourceFiles
    constantSolution.C

\*---------------------------------------------------------------------------*/

#ifndef constantSolution_H
#define constantSolution_H

#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace analyticalSolutions
{

/*---------------------------------------------------------------------------*\
                      Class constantSolution Declaration
\*---------------------------------------------------------------------------*/

class constantSolution
:
    public analyticalSolution
{
    // Private data

        //- The constant value c0
        const scalar c0_;


    // Private Member Functions

        //- No copy construct
        constantSolution(const constantSolution&) = delete;

        //- No copy assignment
        void operator=(const constantSolution&) = delete;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from the coefficients dictionary and space-time velocity
        constantSolution(const dictionary& dict, const vector& A);


    //- Destructor
    virtual ~constantSolution();


    // Member Functions

        // Evaluation at one space-time point

            //- Return u = c0
            virtual scalar value(const point& p) const;

            //- Return the gradient, which is zero
            virtual vector gradient(const point& p) const;

            //- The source is identically zero
            virtual bool hasSource() const
            {
                return false;
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
