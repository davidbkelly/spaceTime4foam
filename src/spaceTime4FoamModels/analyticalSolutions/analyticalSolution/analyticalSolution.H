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
    Foam::analyticalSolution

Description
    Run-time selectable analytical solution u(x, t) of the space-time
    problem div_st(A u) = S, used for boundary data, initial data and error
    norms by both space-time methods.

    Points are space-time points: p.x() is x and p.y() is t (see CLAUDE.md
    section 5). The space-time velocity A = (a, 1, 0) is passed in by the
    caller, built from the keyword a of constant/spaceTimeProperties by
    spaceTimeVelocity(), so the analytical solution always uses the same
    speed as the discretisation.

    Selected from the analyticalSolution sub-dictionary of
    constant/spaceTimeProperties, e.g.
    \verbatim
    a                   1;
    analyticalSolution
    {
        type            travellingSine;
        k               1;
    }
    \endverbatim

SourceFiles
    analyticalSolution.C
    newAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalSolution_H
#define analyticalSolution_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "pointField.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class analyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class analyticalSolution
{
    // Private data

        //- Coefficients dictionary (the analyticalSolution sub-dictionary)
        const dictionary dict_;

        //- Space-time velocity A = (a, 1, 0)
        const vector A_;


    // Private Member Functions

        //- No copy construct
        analyticalSolution(const analyticalSolution&) = delete;

        //- No copy assignment
        void operator=(const analyticalSolution&) = delete;


public:

    //- Runtime type information
    TypeName("analyticalSolution");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            analyticalSolution,
            dictionary,
            (
                const dictionary& dict,
                const vector& A
            ),
            (dict, A)
        );


    // Constructors

        //- Construct from the coefficients dictionary and space-time velocity
        analyticalSolution(const dictionary& dict, const vector& A);


    // Selectors

        //- Select from the analyticalSolution sub-dictionary
        static autoPtr<analyticalSolution> New
        (
            const dictionary& dict,
            const vector& A
        );


    //- Destructor
    virtual ~analyticalSolution();


    // Static Member Functions

        //- Return the space-time velocity A = (a, 1, 0), with a read from
        //  the spaceTimeProperties dictionary
        static vector spaceTimeVelocity(const dictionary& spaceTimeProperties);


    // Member Functions

        // Access

            //- Return the coefficients dictionary
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return the space-time velocity A
            const vector& A() const
            {
                return A_;
            }


        // Evaluation at one space-time point

            //- Return u at the point p = (x, t, z)
            virtual scalar value(const point& p) const = 0;

            //- Return the space-time gradient (du/dx, du/dt, 0) at p
            virtual vector gradient(const point& p) const = 0;

            //- Return true if the source S = div_st(A u) is not identically
            //  zero. Pure virtual, so that every analytical solution must
            //  state it: a model that cannot handle a source (cellCentred)
            //  stops instead of silently dropping it
            virtual bool hasSource() const = 0;

            //- Return the source S = div_st(A u) = A & gradient at p
            virtual scalar source(const point& p) const = 0;


        // Evaluation at a list of space-time points

            //- Return u at the points
            tmp<scalarField> value(const vectorField& points) const;

            //- Return the space-time gradient at the points
            tmp<vectorField> gradient(const vectorField& points) const;

            //- Return the source at the points
            tmp<scalarField> source(const vectorField& points) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
