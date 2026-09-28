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
    perturbSpaceTimeMesh

Description
    Random perturbation of the interior nodes of a one-layer space-time
    triangle mesh (CLAUDE.md section 6), as in Tufillaro et al. (2026):
    every interior node of the front triangulation (a node on no boundary
    edge, see medianDualMesh) is moved by (dx, dt), each drawn independently
    and uniformly from [-maxFraction h, +maxFraction h], with
    h = sqrt(A_domain/(nTriangles/2)). The same displacement is applied to
    the node's back point, so front and back stay paired. Boundary nodes
    are not moved.

    The numbers come from OpenFOAM's Random, seeded with seed; the nodes
    are visited in local (medianDualMesh) node order and dx is drawn
    before dt, so a mesh is reproduced from its seed and maxFraction.
    Portability: Random maps its Rand48 generator (portable) through
    std::uniform_real_distribution, whose algorithm is not fixed by the C++
    standard, so the same seed may give a different mesh with another C++
    standard library (e.g. libstdc++ instead of libc++). The mesh is
    reproducible on one toolchain; tests/perturbSpaceTimeMesh stores
    reference displacements and a checksum to detect any change.

    Before writing, the signed area of every triangle is recomputed: the
    utility stops with a fatal error if any sign changes or any area is
    near zero (below 1e-8 of the mean triangle area).

    Input: system/perturbSpaceTimeMeshDict (optional), e.g.
    \verbatim
    seed            12345;
    maxFraction     0.2;
    \endverbatim
    with the command-line overrides -seed and -maxFraction. seed must be
    given in one of them; maxFraction defaults to 0.2. The seed and
    maxFraction used, h, and the largest displacements are written to the
    log.

    Output: the points are written at full precision (17 digits) to a new
    time directory, or with -overwrite to the directory the points were
    read from (constant/polyMesh), as other OpenFOAM mesh utilities do.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
#include "pointIOField.H"
#include "medianDualMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Randomly move the interior nodes of a one-layer space-time triangle"
        " mesh by up to maxFraction h in x and in t"
    );

    argList::addOption
    (
        "seed",
        "label",
        "Random seed (overrides the dictionary)"
    );

    argList::addOption
    (
        "maxFraction",
        "scalar",
        "Maximum displacement in each direction as a fraction of h"
        " (overrides the dictionary; default 0.2)"
    );

    argList::addOption
    (
        "patch",
        "name",
        "Triangulated patch (default frontAndBack)"
    );

#   include "addOverwriteOption.H"

    // Only serial runs are supported by medianDualMesh
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const bool overwrite = args.found("overwrite");

    // Settings: the dictionary (if present), then the command line
    IOdictionary perturbDict
    (
        IOobject
        (
            "perturbSpaceTimeMeshDict",
            runTime.system(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    label seed = -1;
    bool haveSeed = perturbDict.readIfPresent("seed", seed);

    if (args.readIfPresent("seed", seed))
    {
        haveSeed = true;
    }

    if (!haveSeed)
    {
        FatalErrorInFunction
            << "No seed: give seed in system/perturbSpaceTimeMeshDict or"
            << " the -seed option" << exit(FatalError);
    }

    scalar maxFraction =
        perturbDict.getOrDefault<scalar>("maxFraction", 0.2);
    args.readIfPresent("maxFraction", maxFraction);

    if (maxFraction < 0)
    {
        FatalErrorInFunction
            << "maxFraction must not be negative: " << maxFraction
            << exit(FatalError);
    }

    const word patchName =
        args.getOrDefault<word>("patch", "frontAndBack");

    const medianDualMesh dual(mesh, patchName);

    dual.printSummary();

    const labelList& front = dual.frontMeshPoints();
    const labelList& back = dual.backMeshPoints();
    const boolList& isBoundaryNode = dual.isBoundaryNode();

    // Every mesh point must be a front or a back point
    if (2*dual.nNodes() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "The mesh has " << mesh.nPoints() << " points but the patch "
            << patchName << " has " << 2*dual.nNodes() << ": only one-layer"
            << " meshes are supported" << exit(FatalError);
    }

    const scalar h = dual.h();
    const scalar maxDisplacement = maxFraction*h;

    Info<< nl << "Perturbation settings:" << nl
        << "    seed        " << seed << nl
        << "    maxFraction " << maxFraction << nl
        << "    h           " << h << nl
        << "    max |dx|, max |dt| allowed: maxFraction h = "
        << maxDisplacement << endl;

    Random rnd(seed);

    pointField newPoints(mesh.points());
    pointField newNodePoints(dual.points());

    label nMoved = 0;
    scalar maxDx = 0;
    scalar maxDt = 0;

    forAll(newNodePoints, nodei)
    {
        if (isBoundaryNode[nodei])
        {
            continue;
        }

        // dx first, then dt
        const scalar dx =
            rnd.position<scalar>(-maxDisplacement, maxDisplacement);
        const scalar dt =
            rnd.position<scalar>(-maxDisplacement, maxDisplacement);
        const vector displacement(dx, dt, 0);

        newNodePoints[nodei] += displacement;
        newPoints[front[nodei]] += displacement;
        newPoints[back[nodei]] += displacement;

        maxDx = max(maxDx, mag(dx));
        maxDt = max(maxDt, mag(dt));
        nMoved++;
    }

    // Signed areas before and after: same sign, none near zero
    const scalarField oldAreas
    (
        medianDualMesh::signedTriangleAreas(dual.triangles(), dual.points())
    );
    const scalarField newAreas
    (
        medianDualMesh::signedTriangleAreas(dual.triangles(), newNodePoints)
    );

    const scalar areaTol = 1e-8*sum(mag(oldAreas))/oldAreas.size();
    scalar minRatio = GREAT;
    scalar maxRatio = 0;

    forAll(newAreas, trii)
    {
        const scalar ratio = newAreas[trii]/oldAreas[trii];

        if (ratio <= 0 || mag(newAreas[trii]) <= areaTol)
        {
            FatalErrorInFunction
                << "Triangle " << trii << " " << dual.triangles()[trii]
                << ": the signed area changes from " << oldAreas[trii]
                << " to " << newAreas[trii] << " (tolerance " << areaTol
                << "). Nothing is written; use a smaller maxFraction."
                << exit(FatalError);
        }

        minRatio = min(minRatio, ratio);
        maxRatio = max(maxRatio, ratio);
    }

    Info<< nl << "Perturbation:" << nl
        << "    interior nodes moved " << nMoved << ", boundary nodes kept "
        << dual.nNodes() - nMoved << nl
        << "    max |dx| = " << maxDx << " = " << maxDx/h << " h" << nl
        << "    max |dt| = " << maxDt << " = " << maxDt/h << " h" << nl
        << "    domain area before " << sum(mag(oldAreas))
        << ", after " << sum(mag(newAreas)) << nl
        << "    triangle area ratio new/old: min " << minRatio << ", max "
        << maxRatio << endl;

    // Write the points, as transformPoints does, at full precision
    if (!overwrite)
    {
        ++runTime;
    }

    fileName instance = runTime.timeName();

    if (overwrite)
    {
        instance = mesh.pointsInstance();
    }

    pointIOField pointsOut
    (
        IOobject
        (
            "points",
            instance,
            polyMesh::meshSubDir,
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        newPoints
    );

    IOstream::minPrecision(17);

    Info<< nl << "Writing points into directory "
        << runTime.relativePath(pointsOut.path()) << endl;

    if (!pointsOut.write())
    {
        FatalErrorInFunction
            << "Cannot write " << pointsOut.objectPath() << exit(FatalError);
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
