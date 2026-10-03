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
    Test-vertexCentredResidual

Description
    Test of the residual of the vertexCentred spaceTimeModel
    (vertexCentred::calcResidual) on the mesh of a case, with the settings
    of its constant/spaceTimeProperties (CLAUDE.md section 6, verification
    tests 4 and 5). The nodal values are the analytical solution at the
    nodes. The application exits with status 1 if a check fails.

    -freestream (test 4): analyticalSolution constant (u = c0). Every node,
    interior and boundary, must have |Res_j| <= 1e-13. With weak boundary
    treatment this checks the closure of the dual cells including the
    boundary fluxes of area |n_B|/2.

    -linear (test 5): analyticalSolution linear (u = c0 + c . x) and
    reconstruction linear (VC-2). Res_j/V_j = div_h(A u)_j - f_j with
    f_j = A & c, so the discrete flux divergence div_h(A u)_j equals A & c
    exactly if and only if Res_j/V_j = 0. Checked at every interior node:
        |Res_j/V_j| <= 100 eps |A| max|u| / h
    (the round-off of the edge fluxes, about eps |A| max|u| h, divided by
    V_j, about h^2). At the boundary nodes the nodal weak closure is not
    linearly exact: max |Res_j/V_j| there is reported, not checked. The LSQ
    gradient must equal c to 100 eps max|u|/h at every node (as in
    Test-medianDualMesh, check 9).

    -quadratic (reconstruction consistency): a quadratic analytical
    solution (e.g. mmsQuadratic) with VC-2. At a node whose edge vectors
    come in pairs d, -d (a point-symmetric stencil, e.g. an interior node
    of the structured meshes not next to the boundary) the unweighted LSQ
    gradient of a quadratic is exact. If both ends of every edge of node j
    have such stencils, the reconstruction gives
        uL = u_j + 0.5 grad(u)_j . d = u_k - 0.5 grad(u)_k . d = uR
    (d = x_k - x_j), so the upwind term vanishes and
        Res_j = sum_k (n_jk . A) (u_j + 0.5 grad(u)_j . d_jk) - f_j V_j
    with n_jk oriented out of j. This is computed here from the exact
    gradient, without the model's reconstruction or flux, and compared at
    every such node j (at least one is required):
        |Res_j - Res_j,expected|/V_j <= 100 eps |A| max|u| / h.
    It checks the reconstruction against the exact gradient (a sign or
    factor error) and the source. It cannot detect a reconstruction that
    uses the other node's gradient: on these translation-invariant meshes
    with a constant Hessian the resulting extra edge flux is the same for
    every edge of a given direction, so its divergence is zero. -reference
    detects that.

    -reference (independent residual): any analytical solution and any
    settings. The residual is recomputed here, node by node instead of
    edge by edge, with its own unweighted LSQ gradients (2 x 2 normal
    equations of each node's edge neighbours, solved by Cramer's rule),
    the reconstruction of each side from its own node's gradient (VC-2) or
    none (VC-1), the upwind flux, the weak boundary fluxes (inflow from the
    sign of A . n_B, area |n_B|/2), the source and the strong treatment.
    Only the geometry n_jk, n_B and V_j of medianDualMesh (verified by
    Test-medianDualMesh) is shared. Checked at every node:
        |Res_j - Res_j,reference|/V_j <= 100 eps |A| max|u| / h.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "vertexCentred.H"

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


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Freestream (-freestream), linear-exactness (-linear), reconstruction"
        " (-quadratic) or independent-residual (-reference) test of the"
        " vertexCentred residual on the mesh of the case"
    );

    argList::addBoolOption("freestream", "Freestream preservation (test 4)");
    argList::addBoolOption("linear", "Linear exactness (test 5)");
    argList::addBoolOption
    (
        "quadratic",
        "Reconstruction consistency for a quadratic solution"
    );
    argList::addBoolOption
    (
        "reference",
        "Comparison with an independent residual computed here"
    );

    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    const bool freestream = args.found("freestream");
    const bool linear = args.found("linear");
    const bool quadratic = args.found("quadratic");
    const bool reference = args.found("reference");

    if
    (
        label(freestream) + label(linear) + label(quadratic)
      + label(reference) != 1
    )
    {
        FatalErrorInFunction
            << "Give exactly one of -freestream, -linear, -quadratic and"
            << " -reference"
            << exit(FatalError);
    }

    spaceTimeModels::vertexCentred model(runTime);

    const medianDualMesh& dual = model.dual();
    const scalarField& V = dual.V();
    const boolList& isBoundaryNode = dual.isBoundaryNode();
    const word solutionType(model.analytical().type());

    const scalar eps = std::numeric_limits<scalar>::epsilon();

    // Nodal values: the analytical solution at the nodes
    const scalarField u(model.analytical().value(dual.points()));

    scalarField res(dual.nNodes());
    model.calcResidual(u, res);

    label nFailed = 0;

    Info<< nl << "Mesh: " << dual.nNodes() << " nodes, h = " << dual.h()
        << ", reconstruction "
        << (model.linearReconstruction() ? "linear" : "none")
        << ", boundaryTreatment "
        << (model.strongBoundary() ? "strong" : "weak") << nl << endl;

    if (freestream)
    {
        if (solutionType != "constant")
        {
            FatalErrorInFunction
                << "-freestream needs the analyticalSolution constant, not "
                << solutionType << exit(FatalError);
        }

        scalar maxInterior = 0;
        scalar maxBoundary = 0;

        forAll(res, nodei)
        {
            if (isBoundaryNode[nodei])
            {
                maxBoundary = max(maxBoundary, mag(res[nodei]));
            }
            else
            {
                maxInterior = max(maxInterior, mag(res[nodei]));
            }
        }

        Info<< "Freestream u = " << u[0] << ":" << nl
            << "    max |Res_j|: interior nodes " << maxInterior
            << ", boundary nodes " << maxBoundary << endl;

        nFailed += check
        (
            "max |Res_j| over all nodes",
            max(maxInterior, maxBoundary),
            1e-13
        );
    }
    else if (linear)
    {
        if (solutionType != "linear")
        {
            FatalErrorInFunction
                << "-linear needs the analyticalSolution linear, not "
                << solutionType << exit(FatalError);
        }

        if (!model.linearReconstruction())
        {
            FatalErrorInFunction
                << "-linear needs reconstruction linear (VC-2)"
                << exit(FatalError);
        }

        // The exact gradient c and the exact divergence A & c
        const vector c(model.analytical().gradient(dual.points()[0]));
        const scalar divExact = model.A() & c;
        const scalar maxU = max(mag(u));

        if (divExact == 0)
        {
            FatalErrorInFunction
                << "-linear needs c not perpendicular to A: A & c = 0"
                << exit(FatalError);
        }

        scalar maxInterior = 0;
        scalar maxBoundary = 0;

        forAll(res, nodei)
        {
            const scalar r = mag(res[nodei]/V[nodei]);

            if (isBoundaryNode[nodei])
            {
                maxBoundary = max(maxBoundary, r);
            }
            else
            {
                maxInterior = max(maxInterior, r);
            }
        }

        const vectorField gradU(dual.gradient(u));
        const scalar maxGradError = max(mag(gradU - c));

        const scalar tolRes = 100*eps*mag(model.A())*maxU/dual.h();
        const scalar tolGrad = 100*eps*maxU/dual.h();

        Info<< "Linear u = c0 + c . x with c = " << c << ", A & c = "
            << divExact << ", max |u| = " << maxU << nl
            << "    max |grad(u)_j - c| = " << maxGradError << nl
            << "    max |Res_j/V_j| = max |div_h(A u)_j - A & c|:" << nl
            << "        interior nodes " << maxInterior << nl
            << "        boundary nodes " << maxBoundary
            << " (REPORTED ONLY: the nodal weak closure is not linearly"
            << " exact; relative to |A & c|: "
            << maxBoundary/mag(divExact) << ")"
            << endl;

        nFailed += check("max |grad(u)_j - c|", maxGradError, tolGrad);
        nFailed += check
        (
            "max |Res_j/V_j| over interior nodes",
            maxInterior,
            tolRes
        );
    }
    else if (quadratic)
    {
        if (!model.linearReconstruction())
        {
            FatalErrorInFunction
                << "-quadratic needs reconstruction linear (VC-2)"
                << exit(FatalError);
        }

        const pointField& p = dual.points();
        const edgeList& edges = dual.edges();
        const vectorField& n = dual.n();
        const labelListList& nodeEdges = dual.nodeEdges();
        const vector& A = model.A();

        // Point-symmetric stencils: every edge vector d has a partner -d
        boolList isSymmetric(dual.nNodes(), false);

        forAll(nodeEdges, nodei)
        {
            if (isBoundaryNode[nodei])
            {
                continue;
            }

            bool symmetric = true;

            forAll(nodeEdges[nodei], ei)
            {
                const edge& e = edges[nodeEdges[nodei][ei]];
                const vector d = p[e.otherVertex(nodei)] - p[nodei];
                bool found = false;

                forAll(nodeEdges[nodei], ej)
                {
                    const edge& f = edges[nodeEdges[nodei][ej]];
                    const vector dOther = p[f.otherVertex(nodei)] - p[nodei];

                    if (mag(d + dOther) < 1e-9*dual.h())
                    {
                        found = true;
                    }
                }

                if (!found)
                {
                    symmetric = false;
                }
            }

            isSymmetric[nodei] = symmetric;
        }

        // Exact gradients and sources at the nodes
        const vectorField gradExact(model.analytical().gradient(p));
        const scalarField f(model.analytical().source(p));
        const scalar maxU = max(mag(u));

        label nTested = 0;
        scalar maxDiff = 0;

        forAll(nodeEdges, nodei)
        {
            // Node j and all its neighbours must have symmetric stencils
            bool testable = isSymmetric[nodei];

            forAll(nodeEdges[nodei], ei)
            {
                const edge& e = edges[nodeEdges[nodei][ei]];

                if (!isSymmetric[e.otherVertex(nodei)])
                {
                    testable = false;
                }
            }

            if (!testable)
            {
                continue;
            }

            // Central flux with the face value from node j's exact gradient
            scalar expected = 0;

            forAll(nodeEdges[nodei], ei)
            {
                const label edgei = nodeEdges[nodei][ei];
                const edge& e = edges[edgei];
                const label k = e.otherVertex(nodei);
                const vector d = p[k] - p[nodei];

                // n_jk out of node j
                vector nOut = n[edgei];

                if (e.start() != nodei)
                {
                    nOut = -nOut;
                }

                const scalar uFace = u[nodei] + 0.5*(gradExact[nodei] & d);

                expected += (nOut & A)*uFace;
            }

            if (model.analytical().hasSource())
            {
                expected -= f[nodei]*V[nodei];
            }

            maxDiff = max(maxDiff, mag(res[nodei] - expected)/V[nodei]);
            nTested++;
        }

        const scalar tolRes = 100*eps*mag(A)*maxU/dual.h();

        Info<< "Quadratic u (" << solutionType << "), max |u| = " << maxU
            << nl
            << "    nodes with symmetric stencils at both ends of every"
            << " edge: " << nTested << " of " << dual.nNodes() << nl
            << "    max |Res_j - Res_j,expected|/V_j = " << maxDiff << endl;

        if (nTested == 0)
        {
            Info<< "    FAIL: no node with symmetric stencils" << endl;
            nFailed++;
        }

        nFailed += check
        (
            "max |Res_j - Res_j,expected|/V_j",
            maxDiff,
            tolRes
        );
    }
    else
    {
        const pointField& p = dual.points();
        const edgeList& edges = dual.edges();
        const vectorField& n = dual.n();
        const vectorField& nB = dual.nB();
        const labelListList& nodeEdges = dual.nodeEdges();
        const label nInternalEdges = dual.nInternalEdges();
        const vector& A = model.A();
        const scalarField uExact(model.analytical().value(p));
        const scalarField f(model.analytical().source(p));

        // Nodal values that differ from the exact solution, so that the
        // inflow data and the nodal values differ
        scalarField v(dual.nNodes());

        forAll(v, j)
        {
            v[j] = u[j] + 0.1*Foam::sin(37*p[j].x() + 23*p[j].y());
        }

        scalarField resV(dual.nNodes());
        model.calcResidual(v, resV);

        // Unweighted LSQ gradients, node by node, by Cramer's rule
        vectorField grad(dual.nNodes(), Zero);

        if (model.linearReconstruction())
        {
            forAll(nodeEdges, j)
            {
                scalar Mxx = 0;
                scalar Mxy = 0;
                scalar Myy = 0;
                scalar bx = 0;
                scalar by = 0;

                forAll(nodeEdges[j], ei)
                {
                    const label k = edges[nodeEdges[j][ei]].otherVertex(j);
                    const scalar dx = p[k].x() - p[j].x();
                    const scalar dy = p[k].y() - p[j].y();
                    const scalar du = v[k] - v[j];

                    Mxx += dx*dx;
                    Mxy += dx*dy;
                    Myy += dy*dy;
                    bx += dx*du;
                    by += dy*du;
                }

                const scalar det = Mxx*Myy - Mxy*Mxy;

                grad[j] = vector
                (
                    (bx*Myy - by*Mxy)/det,
                    (Mxx*by - Mxy*bx)/det,
                    0
                );
            }
        }

        // Upwind flux with the states uOwn (own side) and uOther, through
        // the area vector nOut pointing out of the own node
        auto flux = [&A](const vector& nOut, scalar uOwn, scalar uOther)
        {
            const scalar lambda = (nOut/mag(nOut)) & A;

            return
                (
                    0.5*lambda*(uOwn + uOther)
                  - 0.5*mag(lambda)*(uOther - uOwn)
                )*mag(nOut);
        };

        scalarField resRef(dual.nNodes(), 0);

        // Interior edge fluxes out of every node
        forAll(nodeEdges, j)
        {
            forAll(nodeEdges[j], ei)
            {
                const label edgei = nodeEdges[j][ei];
                const label k = edges[edgei].otherVertex(j);
                const vector d = p[k] - p[j];

                vector nOut = n[edgei];

                if (edges[edgei].start() != j)
                {
                    nOut = -nOut;
                }

                const scalar uOwn = v[j] + 0.5*(grad[j] & d);
                const scalar uOther = v[k] - 0.5*(grad[k] & d);

                resRef[j] += flux(nOut, uOwn, uOther);
            }
        }

        // Weak closure: half of each boundary edge to each of its nodes
        if (!model.strongBoundary())
        {
            forAll(nB, i)
            {
                const edge& e = edges[nInternalEdges + i];
                const bool inflow = ((nB[i] & A) < 0);

                forAll(e, ei)
                {
                    const label j = e[ei];
                    const scalar ub = inflow ? uExact[j] : v[j];

                    resRef[j] += 0.5*flux(nB[i], v[j], ub);
                }
            }
        }

        if (model.analytical().hasSource())
        {
            resRef -= f*V;
        }

        if (model.strongBoundary())
        {
            forAll(resRef, j)
            {
                if (isBoundaryNode[j])
                {
                    resRef[j] = 0;
                }
            }
        }

        const scalar maxU = max(mag(v));
        const scalar maxDiff = max(mag(resV - resRef)/V);
        const scalar tolRes = 100*eps*mag(A)*maxU/dual.h();

        Info<< "Independent residual (" << solutionType << "), max |u| = "
            << maxU << ", max |Res_j|/V_j = " << max(mag(resV)/V) << nl
            << "    max |Res_j - Res_j,reference|/V_j = " << maxDiff
            << endl;

        nFailed += check
        (
            "max |Res_j - Res_j,reference|/V_j over all nodes",
            maxDiff,
            tolRes
        );
    }

    Info<< nl;

    if (nFailed > 0)
    {
        Info<< "Test-vertexCentredResidual: " << nFailed << " check(s) FAILED"
            << nl << endl;
        return 1;
    }

    Info<< "Test-vertexCentredResidual: all checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
