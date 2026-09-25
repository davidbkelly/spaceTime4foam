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

Application
    spaceTime4Foam

Description
    Steady space-time solver: the space-time discretisation (spaceTimeModel)
    is chosen at run-time in constant/spaceTimeProperties.

    Time counts iterations, not physical time: each time step is one
    iteration of the selected spaceTimeModel.
    - When the model reports convergence, the fields are written at that
      iteration and the run ends with exit status 0.
    - If endTime is reached without convergence, the fields are written at
      endTime, a warning is printed and the exit status is 1, so that
      scripts and tests can detect the failure.

    In both cases postProcessing/spaceTimeSolverInfo.dat records the number
    of iterations, the converged flag (1 or 0), the final initial residual
    and the wall-clock time [s] of the iteration loop.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "clockTime.H"
#include "OFstream.H"
#include "spaceTimeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady space-time solver; the discretisation is selected in"
        " constant/spaceTimeProperties"
    );

    // Only serial runs have been tested
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    // Create the space-time model (reads the mesh and the fields)
    autoPtr<spaceTimeModel> model = spaceTimeModel::New(runTime);

    Info<< nl << "Starting the space-time iterations" << nl << endl;

    // Wall-clock time of the iteration loop only
    const clockTime loopClock;

    bool converged = false;
    label nIterations = 0;

    while (runTime.loop())
    {
        Info<< "Iteration = " << runTime.timeName() << nl << endl;

        converged = model().evolve();
        nIterations++;

        if (converged)
        {
            Info<< nl << "Converged in " << nIterations << " iterations:"
                << " initial residual " << model().residual()
                << " < " << model().convergenceTolerance() << nl << endl;

            model().writeFields();
            break;
        }
        else
        {
            runTime.write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    const scalar wallTime = loopClock.elapsedTime();

    if (!converged)
    {
        model().writeFields();
    }

#   include "writeSpaceTimeSolverInfo.H"

    if (!converged)
    {
        WarningInFunction
            << "NOT CONVERGED: endTime reached after " << nIterations
            << " iterations with initial residual " << model().residual()
            << " >= " << model().convergenceTolerance() << nl
            << "    Exiting with status 1" << nl << endl;

        return 1;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
