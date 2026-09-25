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
    Test-analyticalSolution

Description
    Unit test of the travellingSine analytical solution. No case is needed.

    For several (a, k) pairs and a set of space-time points p = (x, t, 0):
    - value(p) is compared with sin(2 pi k (x - a t)) written out here;
    - gradient(p) is compared with second-order central differences of
      value() in x and t (step d = 1e-5, truncation error about
      d^2 (2 pi k)^3 / 6 < 1e-7 for k <= 2);
    - source(p) is compared with the central-difference PDE residual
      du/dt + a du/dx (zero for this solution);
    - the field overloads are compared with the point versions.

    Prints the maximum error of each check and exits with status 1 if any
    check exceeds its tolerance.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "mathematicalConstants.H"
#include "analyticalSolution.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// check <description> <maximum error> <tolerance>
// Print the result and return 1 on failure, 0 on success
label check(const word& description, const scalar maxError, const scalar tol)
{
    if (maxError <= tol)
    {
        Info<< "    PASS: " << description << ": max error " << maxError
            << " <= " << tol << endl;
        return 0;
    }

    Info<< "    FAIL: " << description << ": max error " << maxError
        << " > " << tol << endl;
    return 1;
}


int main(int argc, char *argv[])
{
    argList::addNote("Unit test of the travellingSine analytical solution");
    argList::noParallel();
    argList::noCheckProcessorDirectories();

    argList args(argc, argv, false, false, false);

    // Central-difference step and tolerances
    const scalar d = 1e-5;
    const scalar valueTol = 1e-14;
    const scalar cdTol = 1e-6;
    const scalar fieldTol = 0;

    // Space-time test points (x, t, 0), including points outside [0, 1]^2
    List<point> points
    ({
        point(0, 0, 0),
        point(0.1, 0.2, 0),
        point(0.37, 0.81, 0),
        point(0.5, 0.5, 0),
        point(0.9, 0.05, 0),
        point(1, 1, 0),
        point(-0.3, 1.7, 0),
        point(2.2, -0.4, 0)
    });

    // (a, k) pairs: a != 1 and k != 1 also catch swapped x and t
    List<Pair<scalar>> cases
    ({
        Pair<scalar>(1, 1),
        Pair<scalar>(0.7, 1),
        Pair<scalar>(1.3, 2),
        Pair<scalar>(-0.5, 0.5)
    });

    label nFailed = 0;

    forAll(cases, casei)
    {
        const scalar a = cases[casei].first();
        const scalar k = cases[casei].second();
        const vector A(a, 1, 0);

        dictionary dict;
        dict.add("type", word("travellingSine"));
        dict.add("k", k);

        Info<< nl << "Case a = " << a << ", k = " << k << endl;

        autoPtr<analyticalSolution> solPtr = analyticalSolution::New(dict, A);
        const analyticalSolution& sol = solPtr();

        scalar maxValueError = 0;
        scalar maxGradXError = 0;
        scalar maxGradTError = 0;
        scalar maxGradZError = 0;
        scalar maxSourceError = 0;

        forAll(points, pointi)
        {
            const point& p = points[pointi];
            const scalar x = p.x();
            const scalar t = p.y();

            // Value against the formula
            const scalar uFormula =
                Foam::sin(constant::mathematical::twoPi*k*(x - a*t));

            maxValueError = max(maxValueError, mag(sol.value(p) - uFormula));

            // Gradient against central differences
            const vector dx(d, 0, 0);
            const vector dt(0, d, 0);

            const scalar dudx =
                (sol.value(p + dx) - sol.value(p - dx))/(2*d);
            const scalar dudt =
                (sol.value(p + dt) - sol.value(p - dt))/(2*d);

            const vector gradU = sol.gradient(p);

            maxGradXError = max(maxGradXError, mag(gradU.x() - dudx));
            maxGradTError = max(maxGradTError, mag(gradU.y() - dudt));
            maxGradZError = max(maxGradZError, mag(gradU.z()));

            // Source against the central-difference PDE residual
            const scalar residual = dudt + a*dudx;

            maxSourceError =
                max(maxSourceError, mag(sol.source(p) - residual));
        }

        nFailed += check("value against formula", maxValueError, valueTol);
        nFailed += check("du/dx against CD", maxGradXError, cdTol);
        nFailed += check("du/dt against CD", maxGradTError, cdTol);
        nFailed += check("du/dz is zero", maxGradZError, 0);
        nFailed += check("source against CD", maxSourceError, cdTol);

        // Field overloads against the point versions
        const pointField pts(points);
        const scalarField values(sol.value(pts));
        const vectorField gradients(sol.gradient(pts));
        const scalarField sources(sol.source(pts));

        scalar maxFieldError = 0;

        forAll(pts, pointi)
        {
            maxFieldError = max
            (
                maxFieldError,
                mag(values[pointi] - sol.value(pts[pointi]))
            );
            maxFieldError = max
            (
                maxFieldError,
                mag(gradients[pointi] - sol.gradient(pts[pointi]))
            );
            maxFieldError = max
            (
                maxFieldError,
                mag(sources[pointi] - sol.source(pts[pointi]))
            );
        }

        nFailed += check("field overloads", maxFieldError, fieldTol);
    }

    Info<< nl;

    if (nFailed > 0)
    {
        Info<< "Test-analyticalSolution: " << nFailed << " check(s) FAILED"
            << nl << endl;
        return 1;
    }

    Info<< "Test-analyticalSolution: all checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
