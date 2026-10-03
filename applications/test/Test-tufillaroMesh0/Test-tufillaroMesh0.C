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
    Test-tufillaroMesh0

Description
    Closed-form check of the converged VC-2 solution with strong Dirichlet
    data on mesh 0 of Tufillaro, Williams and Nishikawa (2026), Table 1: the
    unit square split by both diagonals into 4 triangles, with the corners
    P00 = (0, 0), P10 = (1, 0), P11 = (1, 1), P01 = (0, 1) and the centre
    c = (1/2, 1/2), the only unknown. A = (1, 1) is required.

    The value u_c of the solution written by spaceTime4Foam (latest time) is
    compared with the closed form derived below, which is computed here
    from the analytical corner values and the source only: no part of the
    vertexCentred model or of medianDualMesh is used.

    Derivation. Each triangle has area 1/4, so V_c = 4 (1/4)/3 = 1/3. On the
    edge from c to a corner, n_ck = (1/3) sum_T n_c^T over its two
    triangles, where n_c^T is the outward unit-length normal of the square
    side opposite c:
        n_c,P00 = (-1, -1)/3,  n_c,P11 = (1, 1)/3,
        n_c,P10 = (1, -1)/3,   n_c,P01 = (-1, 1)/3,
    so n_ck . A = -2/3, 2/3, 0 and 0. The upwind flux Phi |n| is
    (n . A) uL if n . A > 0 and (n . A) uR if n . A < 0, so
        Res_c = (2/3) uL(c, P11) - (2/3) uR(c, P00) - f_c/3,
        uL(c, P11) = u_c + (1/4) g_c . (1, 1),
        uR(c, P00) = u00 + (1/4) g00 . (1, 1).
    LSQ at c: the neighbour offsets are (+-1/2, +-1/2), M_c = I and
    g_c . (1, 1) = u11 - u00 (the u_c, u10 and u01 terms cancel).
    LSQ at P00: offsets (1, 0), (0, 1) and (1/2, 1/2), so
    M = [5/4 1/4; 1/4 5/4], (1, 1) M^-1 = (2/3)(1, 1) and
        g00 . (1, 1) = (2/3) (u10 + u01 + u_c - 3 u00),
    which depends on the unknown. Res_c = 0 is therefore linear in u_c:
        (5/6) u_c = (3/4) u00 - (1/4) u11 + (1/6)(u10 + u01) + f_c/2,
        u_c = 0.9 u00 - 0.3 u11 + 0.2 (u10 + u01) + 0.6 f_c.
    Checks: u = 1 and f = 0 give u_c = 1; u = x + y and f = 2 give
    u_c = 1, the exact value.

    Check: |u_c - u_c,closed| <= 1e-10. Also checked: the mesh has 4 cells
    and 10 points, the front and back centre points agree, and the
    vertexCentredCoeffs are reconstruction linear and boundaryTreatment
    strong. Exits with status 1 on failure.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "analyticalSolution.H"

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
        "Closed-form check of VC-2 with strong Dirichlet data on mesh 0 of"
        " Tufillaro et al. (latest time)"
    );

    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    const instantList times = runTime.times();
    runTime.setTime(times.last(), times.size() - 1);

#   include "createMesh.H"

    Info<< "Time = " << runTime.timeName() << nl << endl;

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

    label nFailed = 0;

    const dictionary& coeffs = props.subDict("vertexCentredCoeffs");

    if
    (
        props.get<word>("spaceTimeModel") != "vertexCentred"
     || coeffs.get<word>("reconstruction") != "linear"
     || coeffs.get<word>("boundaryTreatment") != "strong"
    )
    {
        FatalErrorInFunction
            << "The case must use vertexCentred with reconstruction linear"
            << " and boundaryTreatment strong" << exit(FatalError);
    }

    const vector A(analyticalSolution::spaceTimeVelocity(props));

    if (A != vector(1, 1, 0))
    {
        FatalErrorInFunction
            << "The closed form needs A = (1, 1, 0), not " << A
            << exit(FatalError);
    }

    if (mesh.nCells() != 4 || mesh.nPoints() != 10)
    {
        FatalErrorInFunction
            << "Mesh 0 has 4 cells and 10 points, this mesh has "
            << mesh.nCells() << " and " << mesh.nPoints() << exit(FatalError);
    }

    autoPtr<analyticalSolution> analytical = analyticalSolution::New
    (
        props.subDict("analyticalSolution"),
        A
    );

    const pointScalarField u
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh)
    );

    // The two mesh points at (x, y) = (1/2, 1/2)
    const pointField& meshPoints = mesh.points();
    DynamicList<label> centrePoints;

    forAll(meshPoints, pointi)
    {
        if
        (
            mag(meshPoints[pointi].x() - 0.5) < 1e-12
         && mag(meshPoints[pointi].y() - 0.5) < 1e-12
        )
        {
            centrePoints.append(pointi);
        }
    }

    if (centrePoints.size() != 2)
    {
        FatalErrorInFunction
            << "Expected 2 mesh points at (0.5, 0.5), found "
            << centrePoints.size() << exit(FatalError);
    }

    const scalar uc = u[centrePoints[0]];

    nFailed += check
    (
        "|u_front - u_back| at the centre",
        mag(u[centrePoints[1]] - uc),
        0
    );

    // Closed form from the corner data and the source at the centre
    const scalar u00 = analytical().value(point(0, 0, 0));
    const scalar u10 = analytical().value(point(1, 0, 0));
    const scalar u11 = analytical().value(point(1, 1, 0));
    const scalar u01 = analytical().value(point(0, 1, 0));
    const point c(0.5, 0.5, 0);
    const scalar fc = analytical().source(c);

    const scalar ucClosed =
        0.9*u00 - 0.3*u11 + 0.2*(u10 + u01) + 0.6*fc;

    const scalar ucExact = analytical().value(c);

    Info<< "Corner values u00 = " << u00 << ", u10 = " << u10
        << ", u11 = " << u11 << ", u01 = " << u01 << ", f_c = " << fc << nl
        << "u_c (solver)      = " << uc << nl
        << "u_c (closed form) = " << ucClosed << nl
        << "u_c (exact)       = " << ucExact << nl
        << "closed-form nodal error |u_c - u_exact| = "
        << mag(ucClosed - ucExact) << endl;

    nFailed += check
    (
        "|u_c(solver) - u_c(closed form)|",
        mag(uc - ucClosed),
        1e-10
    );

    Info<< nl;

    if (nFailed > 0)
    {
        Info<< "Test-tufillaroMesh0: " << nFailed << " check(s) FAILED" << nl
            << endl;
        return 1;
    }

    Info<< "Test-tufillaroMesh0: all checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
