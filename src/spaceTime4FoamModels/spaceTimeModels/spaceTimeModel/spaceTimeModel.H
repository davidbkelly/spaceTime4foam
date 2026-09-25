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
    Foam::spaceTimeModel

Description
    Virtual base class for space-time discretisations of div_st(A u) = S,
    selected by the spaceTimeModel keyword of constant/spaceTimeProperties:
    \verbatim
    spaceTimeModel       cellCentred;
    a                    1;        // space-time velocity A = (a, 1, 0)
    convergenceTolerance 1e-10;    // stop when u initial residual < this
    analyticalSolution
    {
        type             travellingSine;
        k                1;
    }
    \endverbatim

    The space-time mesh has x as mesh x and t as mesh y. The problem is
    steady, so Time counts iterations: each call to evolve() performs one
    iteration.

    The spaceTimeProperties dictionary is registered on the Time database,
    so that boundary conditions such as spaceTimeAnalyticalFixedValue use
    the same dictionary.

SourceFiles
    spaceTimeModel.C
    newSpaceTimeModel.C

\*---------------------------------------------------------------------------*/

#ifndef spaceTimeModel_H
#define spaceTimeModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class spaceTimeModel Declaration
\*---------------------------------------------------------------------------*/

class spaceTimeModel
{
    // Private data

        //- Reference to the Time database (counts iterations)
        Time& runTime_;

        //- The spaceTimeProperties dictionary, registered on runTime
        IOdictionary spaceTimeProperties_;

        //- Space-time mesh
        fvMesh mesh_;

        //- Space-time velocity A = (a, 1, 0)
        const vector A_;

        //- Converged when the initial residual falls below this value
        const scalar convergenceTolerance_;

        //- Analytical solution
        autoPtr<analyticalSolution> analyticalSolutionPtr_;


    // Private Member Functions

        //- No copy construct
        spaceTimeModel(const spaceTimeModel&) = delete;

        //- No copy assignment
        void operator=(const spaceTimeModel&) = delete;


public:

    //- Runtime type information
    TypeName("spaceTimeModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            spaceTimeModel,
            dictionary,
            (
                Time& runTime
            ),
            (runTime)
        );


    // Constructors

        //- Construct from components
        spaceTimeModel(const word& type, Time& runTime);


    // Selectors

        //- Select from the spaceTimeModel keyword of spaceTimeProperties
        static autoPtr<spaceTimeModel> New(Time& runTime);


    //- Destructor
    virtual ~spaceTimeModel();


    // Member Functions

        // Access

            //- Return const reference to the time
            const Time& runTime() const
            {
                return runTime_;
            }

            //- Return the spaceTimeProperties dictionary
            const dictionary& spaceTimeProperties() const
            {
                return spaceTimeProperties_;
            }

            //- Return const reference to the mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return non-const reference to the mesh
            fvMesh& mesh()
            {
                return mesh_;
            }

            //- Return the space-time velocity A
            const vector& A() const
            {
                return A_;
            }

            //- Return the convergence tolerance on the initial residual
            scalar convergenceTolerance() const
            {
                return convergenceTolerance_;
            }

            //- Return the analytical solution
            const analyticalSolution& analytical() const
            {
                return analyticalSolutionPtr_();
            }


        // Evolution

            //- Perform one iteration. Return true once converged, i.e. when
            //  the initial residual is below convergenceTolerance
            virtual bool evolve() = 0;

            //- Return the initial residual of the latest iteration
            virtual scalar residual() const = 0;


        // Write

            //- Write the fields at the current time (iteration)
            virtual void writeFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
