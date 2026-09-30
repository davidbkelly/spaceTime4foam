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
    and uniformly from [-maxFraction h, +maxFraction h), with
    h = sqrt(A_domain/(nTriangles/2)). The same displacement is applied to
    the node's back point, so front and back stay paired. Boundary nodes
    are not moved.

    The numbers come from OpenFOAM's Rand48 generator, seeded with seed
    (converted to a 32-bit unsigned integer, i.e. modulo 2^32); the nodes
    are visited in local (medianDualMesh) node order and dx is drawn
    before dt, so a mesh is reproduced from its seed and maxFraction.
    Each draw maps one raw Rand48 output r (an integer in [0, 2^31 - 1]:
    Rand48::operator() returns bits 47..17 of the 48-bit state, so 31 bits)
    to a displacement in [-m, m), m = maxFraction h, with
    \verbatim
        x = r/2^31                    (exact, in [0, 1))
        d = m (2 x - 1)               (= -m + (m - (-m)) x)
    \endverbatim
    2 x - 1 is exact in double precision (x has at most 31 significant
    bits), so d is one correctly rounded IEEE-754 product and cannot be
    changed by a fused multiply-add. Portability: Rand48 is the lrand48
    linear congruential generator (a = 0x5DEECE66D, c = 0xB, m = 2^48,
    initial state (seed << 16) | 0x330E) built on
    std::linear_congruential_engine, whose output sequence is fixed exactly
    by the C++ standard; no std distribution is used. The same seed
    therefore gives bit-identical displacements with any conforming
    compiler and standard library, for a given h. h is computed from the
    mesh points by ordinary double arithmetic (medianDualMesh), which a
    compiler may contract into fused multiply-adds; where the triangle
    areas are exact in double precision, as on the structured test meshes
    (h = 1/N exactly), this cannot change h.
    tests/perturbSpaceTimeMesh stores reference displacements and a
    checksum of the point list to detect any change.

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
#include "Rand48.H"
#include "pointIOField.H"
#include "medianDualMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Map the next raw Rand48 output r in [0, 2^31 - 1] to a displacement in
// [-m, m): x = r/2^31, d = m (2 x - 1). 2 x - 1 is exact, so d is a single
// correctly rounded product, the same on every platform (see Description)
static scalar uniformDisplacement(Rand48& generator, const scalar m)
{
    const double twoTo31 = 2147483648.0;

    const uint32_t r = generator();
    const double x = double(r)/twoTo31;
    const double s = 2.0*x - 1.0;

    return m*s;
}


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

    // The seed is converted to 32 bits (modulo 2^32), as Random does
    Rand48 generator(static_cast<uint32_t>(seed));

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
        const scalar dx = uniformDisplacement(generator, maxDisplacement);
        const scalar dt = uniformDisplacement(generator, maxDisplacement);
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
