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
    Test-medianDualMesh

Description
    Test of medianDualMesh on the mesh of a case (CLAUDE.md section 6,
    verification tests 2, 3 and the LSQ part of 5). Every check prints the
    value found and the tolerance used, and the application exits with
    status 1 if any check fails.

    eps is the machine epsilon. Unless stated otherwise the tolerance is
    1e3 eps times the local scale named in each check.

    Checks on any mesh:
     1. Front/back pairing: every node has a distinct back mesh point, not
        a front point, with the same (x, y) (scale: the (x, y) bounding-box
        diagonal).
     2. Sum of V_j equals the domain area sum |T| (scale A_domain), and the
        boundary integral 0.5 sum_B m_B . n_B (m_B the edge midpoint, exact
        for a polygon), and the value of -domainArea if given.
     3. Closure at every node j:
            sum_k n_jk (oriented out of j) + 0.5 sum_(B at j) n_B = 0
        (scale: the longest edge at j).
     4. V_j = 0.25 sum_k (p_k - p_j) . n_jk at every node (scale V_j). At
        interior nodes this is Eq. 20 of Tufillaro et al.; it also holds at
        boundary nodes, because n_B is normal to (p_k - p_j) there.
     5. Antisymmetry: n_kj is recomputed here independently from the
        higher-labelled node k, as (1/3) sum_T n_k^T + n_B/6 with normals
        and orientations computed here, and must equal -n_jk (scale: the
        edge length). This checks the boundary coefficient 1/6.
     6. n_jk . (p_k - p_j) > 0 on every edge (the minimum cosine is shown).
     7. n_B . S_f > 0 on every boundary edge, with S_f the area vector of
        its mesh face (outward by OpenFOAM convention), and that face
        contains both points of the edge.
     8. Boundary edges per patch equal the number of faces of that patch,
        for every non-empty patch other than the triangulated one.
     9. LSQ exactness: for u = c0 + c . x with c0 = 0.3, c = (1.7, -2.3),
        the gradient equals c at every node, reported for interior,
        boundary and corner nodes (boundary nodes on two patches). Scale:
        |c| (absolute tolerance 1e3 eps |c|, about 6e-13; the observed
        errors, from rounding u_k - u_j, are about eps max|u|/h).
    9b. LSQ weighting: for the nonlinear field
        u = sin(3x) cos(2t) + x^2 t, the gradient equals the unweighted
        least-squares gradient solved here from the 2 x 2 normal equations
        of each node's edge neighbours (Cramer's rule), independently of
        the precomputed coefficients. Linear data cannot detect a weighted
        fit; this can. Scale: the largest slope |u_k - u_j|/|p_k - p_j| at
        the node.
    The valence histogram (number of nodes per number of edge neighbours)
    is printed for interior and boundary nodes.

    With -rightTriangle (single triangle (0, 0), (1, 0), (0, 1)), all with
    tolerance 10 eps:
    10. 3 nodes, 3 edges (all boundary), 1 triangle.
    11. V = 1/6 at every node.
    12. n from (0, 0) to (1, 0) is (1/3, 1/6), from (0, 0) to (0, 1) is
        (1/6, 1/3) and from (1, 0) to (0, 1) is (-1/6, 1/6). The nodes are
        found by coordinates and each expected vector is flipped if the
        stored edge runs the other way (j < k).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "emptyPolyPatch.H"
#include "EdgeMap.H"
#include "Map.H"
#include "HashSet.H"
#include "medianDualMesh.H"

#include <limits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Print the result of a check "value <= tol" and return 1 on failure
label check(const string& description, const scalar value, const scalar tol)
{
    if (value <= tol)
    {
        Info<< "    PASS: " << description.c_str() << ": " << value
            << " <= " << tol << endl;
        return 0;
    }

    Info<< "    FAIL: " << description.c_str() << ": " << value
        << " > " << tol << endl;
    return 1;
}


// Print the result of a check "value > 0" and return 1 on failure
label checkPositive(const string& description, const scalar value)
{
    if (value > 0)
    {
        Info<< "    PASS: " << description.c_str() << ": " << value
            << " > 0" << endl;
        return 0;
    }

    Info<< "    FAIL: " << description.c_str() << ": " << value
        << " <= 0" << endl;
    return 1;
}


// Print the result of a check "value == expected" and return 1 on failure
label checkEqual
(
    const string& description,
    const label value,
    const label expected
)
{
    if (value == expected)
    {
        Info<< "    PASS: " << description.c_str() << ": " << value << endl;
        return 0;
    }

    Info<< "    FAIL: " << description.c_str() << ": " << value
        << ", expected " << expected << endl;
    return 1;
}


// Normal of the segment (a, b), length |b - a|, oriented away from "from".
// Written independently of medianDualMesh::normalAwayFrom: the other
// rotation is used and then oriented by the same kind of test.
vector independentNormal(const point& a, const point& b, const point& from)
{
    const vector d = b - a;
    vector normal(-d.y(), d.x(), 0);

    if (((0.5*(a + b) - from) & normal) < 0)
    {
        normal = -normal;
    }

    return normal;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test of the median-dual geometry (medianDualMesh) on the mesh of"
        " the case"
    );

    argList::addOption
    (
        "patch",
        "name",
        "Triangulated patch (default frontAndBack)"
    );

    argList::addOption
    (
        "domainArea",
        "A",
        "Expected domain area, checked against sum V_j"
    );

    argList::addBoolOption
    (
        "rightTriangle",
        "Also check the single right triangle (0, 0), (1, 0), (0, 1)"
    );

    // Only serial runs are supported by medianDualMesh
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const word patchName =
        args.getOrDefault<word>("patch", "frontAndBack");

    const medianDualMesh dual(mesh, patchName);

    dual.printSummary();

    const pointField& p = dual.points();
    const edgeList& edges = dual.edges();
    const scalarField& V = dual.V();
    const vectorField& n = dual.n();
    const vectorField& nB = dual.nB();
    const faceList& triangles = dual.triangles();
    const label nInternalEdges = dual.nInternalEdges();
    const boolList& isBoundaryNode = dual.isBoundaryNode();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    const scalar eps = std::numeric_limits<scalar>::epsilon();
    const scalar tolFactor = 1e3*eps;

    Info<< nl << "Machine epsilon eps = " << eps << ", tolerance factor"
        << " 1e3 eps = " << tolFactor << nl << endl;

    label nFailed = 0;

    // Longest and shortest edge at each node
    scalarField maxEdgeLength(dual.nNodes(), 0);
    scalarField minEdgeLength(dual.nNodes(), GREAT);

    forAll(edges, edgei)
    {
        const scalar length =
            mag(p[edges[edgei].end()] - p[edges[edgei].start()]);

        forAll(edges[edgei], ei)
        {
            const label nodei = edges[edgei][ei];
            maxEdgeLength[nodei] = max(maxEdgeLength[nodei], length);
            minEdgeLength[nodei] = min(minEdgeLength[nodei], length);
        }
    }


    // 1. Front/back pairing

    Info<< "1. Front/back pairing" << endl;
    {
        const labelList& front = dual.frontMeshPoints();
        const labelList& back = dual.backMeshPoints();
        const labelHashSet frontSet(front);
        labelHashSet backSet;
        label nBad = 0;
        scalar maxDiff = 0;

        forAll(back, nodei)
        {
            if (frontSet.found(back[nodei]) || !backSet.insert(back[nodei]))
            {
                nBad++;
            }

            const vector d =
                mesh.points()[back[nodei]] - mesh.points()[front[nodei]];
            maxDiff = max(maxDiff, max(mag(d.x()), mag(d.y())));
        }

        const scalar L = mag(boundBox(p).span());
        const label nPatchPoints = bMesh[bMesh.findPatchID(patchName)]
            .meshPoints().size();

        nFailed += checkEqual
        (
            "back points that are front points or used twice",
            nBad,
            0
        );
        nFailed += checkEqual
        (
            "patch points, 2 nNodes",
            nPatchPoints,
            2*dual.nNodes()
        );
        nFailed += check
        (
            "max |x_back - x_front|, |y_back - y_front|",
            maxDiff,
            tolFactor*L
        );
    }


    // 2. Sum of the dual areas

    Info<< nl << "2. Sum of the dual areas" << endl;
    {
        const scalar sumV = sum(V);
        const scalar sumT = sum(dual.triangleAreas());

        scalar boundaryArea = 0;

        forAll(nB, i)
        {
            const edge& e = edges[nInternalEdges + i];
            boundaryArea += 0.5*((0.5*(p[e.start()] + p[e.end()])) & nB[i]);
        }

        Info<< "    sum V_j = " << sumV << ", sum |T| = " << sumT
            << ", 0.5 sum_B m_B . n_B = " << boundaryArea << endl;

        nFailed += check
        (
            "|sum V_j - sum |T||",
            mag(sumV - sumT),
            tolFactor*sumT
        );
        nFailed += check
        (
            "|sum V_j - 0.5 sum_B m_B . n_B|",
            mag(sumV - boundaryArea),
            tolFactor*sumT
        );

        scalar expectedArea = 0;

        if (args.readIfPresent("domainArea", expectedArea))
        {
            nFailed += check
            (
                "|sum V_j - domainArea| (domainArea = "
              + Foam::name(expectedArea) + ")",
                mag(sumV - expectedArea),
                tolFactor*expectedArea
            );
        }
    }


    // 3. Closure and 4. volume identity

    Info<< nl << "3. Closure and 4. V_j = 0.25 sum_k (p_k - p_j) . n_jk"
        << endl;
    {
        vectorField closure(dual.nNodes(), Zero);
        scalarField volume(dual.nNodes(), 0);

        forAll(edges, edgei)
        {
            const label j = edges[edgei].start();
            const label k = edges[edgei].end();
            const vector d = p[k] - p[j];

            closure[j] += n[edgei];
            closure[k] -= n[edgei];

            // (p_k - p_j) . n_jk = (p_j - p_k) . n_kj
            volume[j] += 0.25*(d & n[edgei]);
            volume[k] += 0.25*(d & n[edgei]);
        }

        forAll(nB, i)
        {
            const edge& e = edges[nInternalEdges + i];

            closure[e.start()] += 0.5*nB[i];
            closure[e.end()] += 0.5*nB[i];
        }

        scalar maxClosureInterior = 0;
        scalar maxClosureBoundary = 0;
        scalar maxRelClosure = 0;
        scalar maxVolumeInterior = 0;
        scalar maxVolumeBoundary = 0;

        forAll(closure, nodei)
        {
            const scalar c = mag(closure[nodei]);
            const scalar relC = c/maxEdgeLength[nodei];
            const scalar relV = mag(volume[nodei] - V[nodei])/V[nodei];

            maxRelClosure = max(maxRelClosure, relC);

            if (isBoundaryNode[nodei])
            {
                maxClosureBoundary = max(maxClosureBoundary, c);
                maxVolumeBoundary = max(maxVolumeBoundary, relV);
            }
            else
            {
                maxClosureInterior = max(maxClosureInterior, c);
                maxVolumeInterior = max(maxVolumeInterior, relV);
            }
        }

        Info<< "    max |closure|: interior nodes " << maxClosureInterior
            << ", boundary nodes " << maxClosureBoundary << endl;

        nFailed += check
        (
            "3. max |closure_j|/(longest edge at j), all nodes",
            maxRelClosure,
            tolFactor
        );
        nFailed += check
        (
            "4. max |V_j - 0.25 sum (p_k - p_j) . n_jk|/V_j, interior nodes",
            maxVolumeInterior,
            tolFactor
        );
        nFailed += check
        (
            "4. max |V_j - 0.25 sum (p_k - p_j) . n_jk|/V_j, boundary nodes",
            maxVolumeBoundary,
            tolFactor
        );
    }


    // 5. Antisymmetry from the higher-labelled node

    Info<< nl << "5. Antisymmetry: n_kj recomputed from node k" << endl;
    {
        EdgeMap<label> edgeLookup(2*edges.size());

        forAll(edges, edgei)
        {
            edgeLookup.insert(edges[edgei], edgei);
        }

        vectorField nFromEnd(edges.size(), Zero);
        labelList nTrianglesOfEdge(edges.size(), Zero);
        labelList oppositeOfEdge(edges.size(), -1);
        label nMissing = 0;

        forAll(triangles, trii)
        {
            const face& f = triangles[trii];

            forAll(f, fp)
            {
                const label a = f[fp];
                const label b = f.nextLabel(fp);
                const label l = f.prevLabel(fp);
                const label j = min(a, b);
                const label k = max(a, b);

                const auto iter = edgeLookup.cfind(edge(j, k));

                if (!iter.good())
                {
                    nMissing++;
                    continue;
                }

                const label edgei = iter.val();

                // n_k^T: normal of the edge (j, l) opposite k, away from k
                nFromEnd[edgei] += independentNormal(p[j], p[l], p[k])/3.0;
                nTrianglesOfEdge[edgei]++;
                oppositeOfEdge[edgei] = l;
            }
        }

        // Boundary edges: one triangle. They must be the last edges.
        label nWrongBoundary = 0;

        forAll(edges, edgei)
        {
            const bool boundary = (nTrianglesOfEdge[edgei] == 1);

            if (boundary != (edgei >= nInternalEdges))
            {
                nWrongBoundary++;
            }

            if (boundary)
            {
                const label j = edges[edgei].start();
                const label k = edges[edgei].end();
                const vector nBIndependent =
                    independentNormal(p[j], p[k], p[oppositeOfEdge[edgei]]);

                nFromEnd[edgei] += nBIndependent/6.0;
            }
        }

        scalar maxRelDiff = 0;
        scalar maxDiff = 0;

        forAll(edges, edgei)
        {
            const scalar length =
                mag(p[edges[edgei].end()] - p[edges[edgei].start()]);
            const scalar diff = mag(nFromEnd[edgei] + n[edgei]);

            maxDiff = max(maxDiff, diff);
            maxRelDiff = max(maxRelDiff, diff/length);
        }

        Info<< "    max |n_kj + n_jk| = " << maxDiff << endl;

        nFailed += checkEqual("triangle edges not found", nMissing, 0);
        nFailed += checkEqual
        (
            "edges with one triangle that are not the last edges, or"
            " the other way round",
            nWrongBoundary,
            0
        );
        nFailed += check
        (
            "max |n_kj + n_jk|/|p_k - p_j|",
            maxRelDiff,
            tolFactor
        );
    }


    // 6. Orientation of n_jk

    Info<< nl << "6. n_jk . (p_k - p_j) > 0" << endl;
    {
        scalar minCos = GREAT;

        forAll(edges, edgei)
        {
            const vector d = p[edges[edgei].end()] - p[edges[edgei].start()];

            minCos = min(minCos, (n[edgei] & d)/(mag(n[edgei])*mag(d)));
        }

        nFailed += checkPositive
        (
            "min n_jk . (p_k - p_j)/(|n_jk| |p_k - p_j|)",
            minCos
        );
    }


    // 7. and 8. Boundary edges and patches

    Info<< nl << "7. n_B . S_f > 0 and 8. boundary edges per patch" << endl;
    {
        const labelList& faces = dual.boundaryEdgeFaces();
        const labelList& front = dual.frontMeshPoints();
        const vectorField& Sf = mesh.faceAreas();

        scalar minCos = GREAT;
        label nNotOnFace = 0;
        labelList nEdgesOfPatch(bMesh.size(), Zero);

        forAll(nB, i)
        {
            const edge& e = edges[nInternalEdges + i];
            const label facei = faces[i];
            const face& meshFace = mesh.faces()[facei];

            if
            (
                !meshFace.found(front[e.start()])
             || !meshFace.found(front[e.end()])
            )
            {
                nNotOnFace++;
            }

            minCos = min
            (
                minCos,
                (nB[i] & Sf[facei])/(mag(nB[i])*mag(Sf[facei]))
            );

            nEdgesOfPatch[dual.boundaryEdgePatches()[i]]++;
        }

        nFailed += checkEqual
        (
            "boundary edges whose face does not contain both points",
            nNotOnFace,
            0
        );
        nFailed += checkPositive
        (
            "7. min n_B . S_f/(|n_B| |S_f|)",
            minCos
        );

        forAll(bMesh, patchi)
        {
            if
            (
                bMesh[patchi].name() == patchName
             || isA<emptyPolyPatch>(bMesh[patchi])
            )
            {
                continue;
            }

            nFailed += checkEqual
            (
                "8. boundary edges on patch " + bMesh[patchi].name()
              + " (= number of patch faces " + Foam::name(bMesh[patchi].size())
              + ")",
                nEdgesOfPatch[patchi],
                bMesh[patchi].size()
            );
        }
    }


    // 9. LSQ exactness for linear data

    Info<< nl << "9. LSQ gradient of u = c0 + c . x" << endl;
    {
        const scalar c0 = 0.3;
        const vector c(1.7, -2.3, 0);

        scalarField u(dual.nNodes());

        forAll(p, nodei)
        {
            u[nodei] = c0 + (c & p[nodei]);
        }

        const vectorField grad(dual.gradient(u));
        const scalar maxU = max(mag(u));

        // Absolute scale |c|: the error comes from rounding u_k - u_j
        // (about eps max|u|/h, 1e-14 at h = 1/32)
        const scalar tolLinear = tolFactor*mag(c);

        // Corner nodes: boundary nodes with boundary edges on two patches
        labelList firstPatch(dual.nNodes(), -1);
        boolList isCorner(dual.nNodes(), false);

        forAll(nB, i)
        {
            const edge& e = edges[nInternalEdges + i];
            const label patchi = dual.boundaryEdgePatches()[i];

            forAll(e, ei)
            {
                const label nodei = e[ei];

                if (firstPatch[nodei] < 0)
                {
                    firstPatch[nodei] = patchi;
                }
                else if (firstPatch[nodei] != patchi)
                {
                    isCorner[nodei] = true;
                }
            }
        }

        scalar maxErrInterior = 0;
        scalar maxErrBoundary = 0;
        scalar maxErrCorner = 0;
        scalar maxErr = 0;
        scalar maxZ = 0;
        label nCorners = 0;

        forAll(grad, nodei)
        {
            const scalar err = mag(grad[nodei] - c);
            maxErr = max(maxErr, err);
            maxZ = max(maxZ, mag(grad[nodei].z()));

            if (isCorner[nodei])
            {
                maxErrCorner = max(maxErrCorner, err);
                nCorners++;
            }

            if (isBoundaryNode[nodei])
            {
                maxErrBoundary = max(maxErrBoundary, err);
            }
            else
            {
                maxErrInterior = max(maxErrInterior, err);
            }
        }

        Info<< "    c0 = " << c0 << ", c = " << c << ", max |u| = " << maxU
            << nl
            << "    max |grad(u) - c|: interior nodes " << maxErrInterior
            << ", boundary nodes " << maxErrBoundary << ", " << nCorners
            << " corner nodes " << maxErrCorner << endl;

        nFailed += check
        (
            "max |grad(u)_j - c| (tolerance 1e3 eps |c|)",
            maxErr,
            tolLinear
        );
        nFailed += check("max |grad(u)_j . z|", maxZ, 0);
    }


    // 9b. Nonlinear data against an independent unweighted LSQ

    Info<< nl << "9b. LSQ gradient of u = sin(3x) cos(2t) + x^2 t against an"
        << " unweighted LSQ solved here" << endl;
    {
        scalarField u(dual.nNodes());

        forAll(p, nodei)
        {
            const scalar x = p[nodei].x();
            const scalar t = p[nodei].y();

            u[nodei] = Foam::sin(3*x)*Foam::cos(2*t) + sqr(x)*t;
        }

        const vectorField grad(dual.gradient(u));

        // Edge neighbours of each node, from the edge list
        List<DynamicList<label>> neighbours(dual.nNodes());

        forAll(edges, edgei)
        {
            neighbours[edges[edgei].start()].append(edges[edgei].end());
            neighbours[edges[edgei].end()].append(edges[edgei].start());
        }

        scalar maxDiff = 0;
        scalar maxRelDiff = 0;
        scalar maxGrad = 0;

        forAll(p, nodei)
        {
            // Normal equations of the unweighted fit, solved by Cramer's
            // rule
            scalar sxx = 0;
            scalar sxy = 0;
            scalar syy = 0;
            scalar bx = 0;
            scalar by = 0;
            scalar maxSlope = 0;

            forAll(neighbours[nodei], i)
            {
                const label k = neighbours[nodei][i];
                const scalar dx = p[k].x() - p[nodei].x();
                const scalar dy = p[k].y() - p[nodei].y();
                const scalar du = u[k] - u[nodei];

                sxx += dx*dx;
                sxy += dx*dy;
                syy += dy*dy;
                bx += dx*du;
                by += dy*du;
                maxSlope =
                    max(maxSlope, mag(du)/Foam::sqrt(dx*dx + dy*dy));
            }

            const scalar det = sxx*syy - sqr(sxy);
            const scalar gx = (syy*bx - sxy*by)/det;
            const scalar gy = (sxx*by - sxy*bx)/det;

            const scalar diff = Foam::sqrt
            (
                sqr(grad[nodei].x() - gx)
              + sqr(grad[nodei].y() - gy)
              + sqr(grad[nodei].z())
            );

            maxDiff = max(maxDiff, diff);
            maxRelDiff = max(maxRelDiff, diff/(maxSlope + VSMALL));
            maxGrad = max(maxGrad, Foam::sqrt(sqr(gx) + sqr(gy)));
        }

        Info<< "    max |grad(u)| = " << maxGrad
            << ", max |grad(u) - grad_LSQ(u)| = " << maxDiff << endl;

        nFailed += check
        (
            "max |grad(u)_j - grad_LSQ(u)_j|/(max_k |u_k - u_j|/|p_k - p_j|)",
            maxRelDiff,
            tolFactor
        );
    }


    // Valence histogram

    Info<< nl << "Valence histogram (edge neighbours: nodes)" << endl;
    {
        Map<label> interiorValence;
        Map<label> boundaryValence;

        forAll(dual.nodeEdges(), nodei)
        {
            const label valence = dual.nodeEdges()[nodei].size();

            if (isBoundaryNode[nodei])
            {
                boundaryValence(valence)++;
            }
            else
            {
                interiorValence(valence)++;
            }
        }

        Info<< "    interior:";
        for (const label valence : interiorValence.sortedToc())
        {
            Info<< " " << valence << ": " << interiorValence[valence];
        }
        Info<< nl << "    boundary:";
        for (const label valence : boundaryValence.sortedToc())
        {
            Info<< " " << valence << ": " << boundaryValence[valence];
        }
        Info<< nl << "    distinct interior valences: "
            << interiorValence.size() << endl;
    }


    // Right-triangle checks

    if (args.found("rightTriangle"))
    {
        Info<< nl << "Right triangle (0, 0), (1, 0), (0, 1)" << endl;

        const scalar tol = 10*eps;

        nFailed += checkEqual("10. nodes", dual.nNodes(), 3);
        nFailed += checkEqual("10. edges", dual.nEdges(), 3);
        nFailed += checkEqual("10. boundary edges", dual.nBoundaryEdges(), 3);
        nFailed += checkEqual("10. triangles", dual.nTriangles(), 1);

        if (dual.nNodes() == 3 && dual.nEdges() == 3)
        {
            // Nodes by coordinates
            const pointField corners
            ({
                point(0, 0, 0),
                point(1, 0, 0),
                point(0, 1, 0)
            });

            labelList cornerNode(3, -1);

            forAll(corners, i)
            {
                forAll(p, nodei)
                {
                    if (mag(p[nodei] - corners[i]) < 1e-12)
                    {
                        cornerNode[i] = nodei;
                    }
                }
            }

            nFailed += checkEqual
            (
                "10. corner nodes not found",
                label(cornerNode.found(-1)),
                0
            );

            if (!cornerNode.found(-1))
            {
                forAll(cornerNode, i)
                {
                    nFailed += check
                    (
                        "11. |V - 1/6| at " + Foam::name(corners[i]),
                        mag(V[cornerNode[i]] - 1.0/6.0),
                        tol
                    );
                }

                // Expected n from corner a to corner b
                const labelPairList pairs
                ({
                    labelPair(0, 1),
                    labelPair(0, 2),
                    labelPair(1, 2)
                });
                const vectorField expected
                ({
                    vector(1.0/3.0, 1.0/6.0, 0),
                    vector(1.0/6.0, 1.0/3.0, 0),
                    vector(-1.0/6.0, 1.0/6.0, 0)
                });

                forAll(pairs, pairi)
                {
                    const label a = cornerNode[pairs[pairi].first()];
                    const label b = cornerNode[pairs[pairi].second()];

                    label edgei = -1;

                    forAll(edges, ei)
                    {
                        if (edges[ei] == edge(a, b))
                        {
                            edgei = ei;
                        }
                    }

                    if (edgei < 0)
                    {
                        Info<< "    FAIL: edge " << a << " " << b
                            << " not found" << endl;
                        nFailed++;
                        continue;
                    }

                    // Stored from the lower label: flip if that is b
                    vector nExpected = expected[pairi];

                    if (edges[edgei].start() != a)
                    {
                        nExpected = -nExpected;
                    }

                    Info<< "    edge " << corners[pairs[pairi].first()]
                        << " -> " << corners[pairs[pairi].second()]
                        << ": n = " << n[edgei] << " (stored from node "
                        << edges[edgei].start() << " to "
                        << edges[edgei].end() << "), expected "
                        << nExpected << endl;

                    nFailed += check
                    (
                        "12. |n - expected| for this edge",
                        mag(n[edgei] - nExpected),
                        tol
                    );
                }
            }
        }
    }

    Info<< nl;

    if (nFailed > 0)
    {
        Info<< "Test-medianDualMesh: " << nFailed << " check(s) FAILED" << nl
            << endl;
        return 1;
    }

    Info<< "Test-medianDualMesh: all checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
