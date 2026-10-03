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
        "Freestream (-freestream) or linear-exactness (-linear) test of the"
        " vertexCentred residual on the mesh of the case"
    );

    argList::addBoolOption("freestream", "Freestream preservation (test 4)");
    argList::addBoolOption("linear", "Linear exactness (test 5)");

    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    const bool freestream = args.found("freestream");
    const bool linear = args.found("linear");

    if (freestream == linear)
    {
        FatalErrorInFunction
            << "Give exactly one of -freestream and -linear"
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
    else
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
