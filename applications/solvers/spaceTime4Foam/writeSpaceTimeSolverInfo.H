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

Description
    Write postProcessing/spaceTimeSolverInfo.dat: one data line with the
    number of iterations, the converged flag (1 or 0), the final initial
    residual and the wall-clock time [s] of the iteration loop.

Required Variables
    - runTime [Time]
    - model [autoPtr<spaceTimeModel>]
    - nIterations [label]
    - converged [bool]
    - wallTime [scalar]

\*---------------------------------------------------------------------------*/

{
    const fileName infoDir(runTime.globalPath()/"postProcessing");
    mkDir(infoDir);

    OFstream infoFile(infoDir/"spaceTimeSolverInfo.dat");
    infoFile.precision(10);

    infoFile
        << "# spaceTime4Foam solver information" << nl
        << "# spaceTimeModel "
        << model().spaceTimeProperties().get<word>("spaceTimeModel") << nl
        << "# nIterations converged finalInitialResidual wallTime[s]" << nl
        << nIterations << " " << label(converged) << " "
        << model().residual() << " " << wallTime << nl;

    Info<< "Written " << infoFile.name() << nl
        << "    nIterations = " << nIterations
        << ", converged = " << label(converged)
        << ", finalInitialResidual = " << model().residual()
        << ", wallTime = " << wallTime << " s" << endl;
}

// ************************************************************************* //
