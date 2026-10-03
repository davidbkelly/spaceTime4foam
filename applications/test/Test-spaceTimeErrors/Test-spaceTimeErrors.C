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
    Test-spaceTimeErrors

Description
    Known-answer data for the vertexCentred branch of spaceTimeErrors, on
    the mesh of a one-layer case whose front patch (frontAndBack faces at
    the minimum z) is triangulated.

    Writes the pointScalarField u = u_exact + delta at time 1, with the
    analytical solution of constant/spaceTimeProperties and the known
    perturbation
        delta(x, t) = 0.01 + 0.02 x + 0.05 t^2,
    a function of (x, t) only, so front and back points agree.

    Writes expected.dat with the values spaceTimeErrors must report,
    computed here without spaceTimeErrors and without medianDualMesh:
    - dual areas V_j = sum of |T|/3 over the front triangles T of node j;
    - boundary nodes: front points of a face of a non-empty patch other
      than frontAndBack;
    - t = T weights: for each tEnd face, half the distance between its two
      front points, added to each of them;
    - e_j = u_j - u_exact(x_j), and over each subset S
          L1 = sum |e| w / sum w, L2 = sqrt(sum e^2 w / sum w),
          Linf = max |e|.
    Columns of expected.dat (one data line):
        nNodes nInterior nBoundary nUnknownsWeak nUnknownsStrong
        L1 L2 Linf over all, interior, boundary and tEnd nodes (12 values)
    with nUnknownsWeak = nNodes and nUnknownsStrong = nInterior.

    Check: the t = T weights sum to the length of tEnd (the sum of its
    face lengths) and to the expected value given with -tEndLength.
    Exits with status 1 if a check fails.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "emptyPolyPatch.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "OFstream.H"
#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Weighted L1, L2 and Linf of the errors e over the nodes with w > 0
void norms
(
    const scalarField& e,
    const scalarField& w,
    const boolList& inSet,
    scalar& L1,
    scalar& L2,
    scalar& Linf
)
{
    scalar sumW = 0;
    L1 = 0;
    L2 = 0;
    Linf = 0;

    forAll(e, i)
    {
        if (inSet[i])
        {
            L1 += mag(e[i])*w[i];
            L2 += sqr(e[i])*w[i];
            Linf = max(Linf, mag(e[i]));
            sumW += w[i];
        }
    }

    L1 /= sumW;
    L2 = Foam::sqrt(L2/sumW);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write u = u_exact + delta and the expected vertexCentred error"
        " norms of spaceTimeErrors"
    );

    argList::addOption
    (
        "tEndLength",
        "L",
        "Expected length of the tEnd patch (default 1)"
    );

    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const IOdictionary props
    (
        IOobject
        (
            "spaceTimeProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<analyticalSolution> analytical = analyticalSolution::New
    (
        props.subDict("analyticalSolution"),
        analyticalSolution::spaceTimeVelocity(props)
    );

    const pointField& points = mesh.points();
    const scalar zMin = mesh.bounds().min().z();
    const scalar zTol = 1e-6*(mesh.bounds().max().z() - zMin);

    // Front points (z = zMin) and their node numbers
    labelList nodeOfPoint(mesh.nPoints(), -1);
    DynamicList<label> frontPoints;

    forAll(points, pointi)
    {
        if (mag(points[pointi].z() - zMin) < zTol)
        {
            nodeOfPoint[pointi] = frontPoints.size();
            frontPoints.append(pointi);
        }
    }

    const label nNodes = frontPoints.size();

    // u = u_exact + delta at every mesh point, and the nodal errors
    pointScalarField u
    (
        IOobject
        (
            "u",
            "1",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    );

    forAll(points, pointi)
    {
        const point p(points[pointi].x(), points[pointi].y(), 0);
        const scalar delta = 0.01 + 0.02*p.x() + 0.05*sqr(p.y());

        u[pointi] = analytical().value(p) + delta;
    }

    scalarField e(nNodes);

    forAll(frontPoints, nodei)
    {
        const point& q = points[frontPoints[nodei]];
        const point p(q.x(), q.y(), 0);

        e[nodei] = u[frontPoints[nodei]] - analytical().value(p);
    }

    // Dual areas from the front triangles of frontAndBack
    const label fbID = mesh.boundaryMesh().findPatchID("frontAndBack");
    const polyPatch& fb = mesh.boundaryMesh()[fbID];
    scalarField V(nNodes, 0);

    forAll(fb, facei)
    {
        const face& f = fb[facei];

        if (nodeOfPoint[f[0]] < 0)
        {
            continue;
        }

        const scalar area = f.mag(points);

        forAll(f, fp)
        {
            V[nodeOfPoint[f[fp]]] += area/f.size();
        }
    }

    // Boundary nodes and t = T weights from the side patches
    boolList isBoundary(nNodes, false);
    scalarField tEndWeight(nNodes, 0);
    scalar tEndLength = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if (patchi == fbID || isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        forAll(pp, facei)
        {
            const face& f = pp[facei];
            DynamicList<label> frontNodes(2);

            forAll(f, fp)
            {
                if (nodeOfPoint[f[fp]] >= 0)
                {
                    frontNodes.append(nodeOfPoint[f[fp]]);
                    isBoundary[nodeOfPoint[f[fp]]] = true;
                }
            }

            if (pp.name() == "tEnd" && frontNodes.size() == 2)
            {
                const scalar length = mag
                (
                    points[frontPoints[frontNodes[1]]]
                  - points[frontPoints[frontNodes[0]]]
                );

                tEndWeight[frontNodes[0]] += 0.5*length;
                tEndWeight[frontNodes[1]] += 0.5*length;
                tEndLength += length;
            }
        }
    }

    // Subsets
    boolList all(nNodes, true);
    boolList interior(nNodes);
    boolList onTEnd(nNodes);
    label nBoundary = 0;

    forAll(interior, nodei)
    {
        interior[nodei] = !isBoundary[nodei];
        onTEnd[nodei] = (tEndWeight[nodei] > 0);

        if (isBoundary[nodei])
        {
            nBoundary++;
        }
    }

    const label nInterior = nNodes - nBoundary;

    scalarField L1(4);
    scalarField L2(4);
    scalarField Linf(4);

    norms(e, V, all, L1[0], L2[0], Linf[0]);
    norms(e, V, interior, L1[1], L2[1], Linf[1]);
    norms(e, V, isBoundary, L1[2], L2[2], Linf[2]);
    norms(e, tEndWeight, onTEnd, L1[3], L2[3], Linf[3]);

    // Checks on the t = T weights
    label nFailed = 0;
    const scalar expectedLength = args.getOrDefault<scalar>("tEndLength", 1);
    const scalar sumWeights = sum(tEndWeight);

    Info<< "nNodes " << nNodes << ", nInterior " << nInterior
        << ", nBoundary " << nBoundary << ", sum V_j " << sum(V) << nl
        << "t = T weights: sum " << sumWeights << ", tEnd length "
        << tEndLength << ", expected " << expectedLength << endl;

    if (mag(sumWeights - tEndLength) > 1e-12*tEndLength)
    {
        Info<< "    FAIL: sum of t = T weights is not the tEnd length" << endl;
        nFailed++;
    }

    if (mag(sumWeights - expectedLength) > 1e-12*expectedLength)
    {
        Info<< "    FAIL: sum of t = T weights is not " << expectedLength
            << endl;
        nFailed++;
    }

    // Write u and the expected values
    runTime.setTime(1, 1);
    u.write();

    OFstream os("expected.dat");
    os.precision(17);

    os  << "# Test-spaceTimeErrors: expected vertexCentred values" << nl
        << "# nNodes nInterior nBoundary nUnknownsWeak nUnknownsStrong"
        << " L1 L2 Linf (all, interior, boundary, tEnd)" << nl
        << nNodes << " " << nInterior << " " << nBoundary << " " << nNodes
        << " " << nInterior;

    forAll(L1, seti)
    {
        os  << " " << L1[seti] << " " << L2[seti] << " " << Linf[seti];
    }

    os  << nl;

    Info<< "Written 1/u and " << os.name() << nl << endl;

    if (nFailed > 0)
    {
        Info<< "Test-spaceTimeErrors: " << nFailed << " check(s) FAILED" << nl
            << endl;
        return 1;
    }

    Info<< "Test-spaceTimeErrors: all checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
