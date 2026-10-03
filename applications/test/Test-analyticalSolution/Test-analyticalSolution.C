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
    Unit test of every analyticalSolution type: travellingSine, constant,
    linear, mmsSymmetricQuadratic, mmsQuadratic and mmsExponential. No case
    is needed.

    Each type is tested with several coefficient sets and space-time
    velocities A = (a, 1, 0), at a set of space-time points p = (x, t, 0):
    - value(p) is compared with the formula written out here;
    - gradient(p) is compared with second-order central differences of
      value() in x and t (step d = 1e-5; the truncation error is at most
      d^2 (2 pi k)^3 / 6 < 1e-7 for travellingSine with k <= 2, and zero
      for the polynomials);
    - source(p) is compared with the central-difference PDE residual
      du/dt + a du/dx;
    - hasSource() and source() agree: if hasSource() is false, |source(p)|
      must not exceed zeroSourceTol at every point; if it is true, it must
      exceed zeroSourceTol at some point;
    - the field overloads are compared with the point versions.

    For the manufactured solutions with A = (1, 1, 0), source(p) is also
    compared with the forcing functions of Tufillaro, Williams and
    Nishikawa (2026), Eq. 125-127:
        f1 = 3 (x + t), f2 = 6x + 10t, f3 = 2 k exp(k (x + t)).

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
label check(const string& description, const scalar maxError, const scalar tol)
{
    if (maxError <= tol)
    {
        Info<< "    PASS: " << description.c_str() << ": max error "
            << maxError << " <= " << tol << endl;
        return 0;
    }

    Info<< "    FAIL: " << description.c_str() << ": max error " << maxError
        << " > " << tol << endl;
    return 1;
}


// The value of the analytical solution, written out independently of the
// library, for the coefficients in dict
scalar formula
(
    const dictionary& dict,
    const scalar a,
    const point& p
)
{
    const word type(dict.get<word>("type"));
    const scalar x = p.x();
    const scalar t = p.y();

    if (type == "travellingSine")
    {
        const scalar k = dict.get<scalar>("k");
        return Foam::sin(constant::mathematical::twoPi*k*(x - a*t));
    }
    else if (type == "constant")
    {
        return dict.get<scalar>("c0");
    }
    else if (type == "linear")
    {
        const scalarList c(dict.get<scalarList>("c"));
        return dict.get<scalar>("c0") + c[0]*x + c[1]*t;
    }
    else if (type == "mmsSymmetricQuadratic")
    {
        return sqr(x) + x*t + sqr(t);
    }
    else if (type == "mmsQuadratic")
    {
        return 3*sqr(x) + 5*sqr(t);
    }
    else if (type == "mmsExponential")
    {
        const scalar k = dict.getOrDefault<scalar>("k", 0.1);
        return Foam::exp(k*(x + t));
    }

    FatalErrorInFunction
        << "No formula for type " << type << exit(FatalError);

    return 0;
}


// Set f to the forcing function of Tufillaro et al. (Eq. 125-127) for
// A = (1, 1); return false if the type is not one of their solutions
bool paperSource(const dictionary& dict, const point& p, scalar& f)
{
    const word type(dict.get<word>("type"));
    const scalar x = p.x();
    const scalar t = p.y();

    if (type == "mmsSymmetricQuadratic")
    {
        f = 3*(x + t);
        return true;
    }
    else if (type == "mmsQuadratic")
    {
        f = 6*x + 10*t;
        return true;
    }
    else if (type == "mmsExponential")
    {
        const scalar k = dict.getOrDefault<scalar>("k", 0.1);
        f = 2*k*Foam::exp(k*(x + t));
        return true;
    }

    return false;
}


int main(int argc, char *argv[])
{
    argList::addNote("Unit test of the analyticalSolution types");
    argList::noParallel();
    argList::noCheckProcessorDirectories();

    argList args(argc, argv, false, false, false);

    // Central-difference step and tolerances
    const scalar d = 1e-5;
    const scalar valueTol = 1e-14;
    const scalar cdTol = 1e-6;
    const scalar fieldTol = 0;
    const scalar zeroSourceTol = 1e-12;
    const scalar paperTol = 1e-13;

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

    // Coefficient dictionaries, in dictionary syntax
    List<string> dicts
    ({
        "type travellingSine; k 1;",
        "type travellingSine; k 2;",
        "type travellingSine; k 0.5;",
        "type constant; c0 0.7;",
        "type linear; c0 0.3; c (1.7 -2.3);",
        "type linear; c0 -0.2; c (0.4 1.1 0);",
        "type linear; c0 0.5; c (1 -1);",        // A & c = 0 when a = 1
        "type mmsSymmetricQuadratic;",
        "type mmsQuadratic;",
        "type mmsExponential;",
        "type mmsExponential; k 0.3;"
    });

    // Speeds a: a != 1 also catches swapped x and t
    List<scalar> speeds({1, 0.7, 1.3, -0.5});

    label nFailed = 0;
    label nCases = 0;

    forAll(dicts, dicti)
    {
        IStringStream dictStream(dicts[dicti]);
        const dictionary dict(dictStream);

        forAll(speeds, speedi)
        {
            const scalar a = speeds[speedi];
            const vector A(a, 1, 0);

            Info<< nl << "Case {" << dicts[dicti].c_str() << "}, a = " << a
                << endl;

            autoPtr<analyticalSolution> solPtr =
                analyticalSolution::New(dict, A);
            const analyticalSolution& sol = solPtr();

            nCases++;

            scalar maxValueError = 0;
            scalar maxGradXError = 0;
            scalar maxGradTError = 0;
            scalar maxGradZError = 0;
            scalar maxSourceError = 0;
            scalar maxSource = 0;
            scalar maxPaperError = 0;
            bool hasPaperSource = false;

            forAll(points, pointi)
            {
                const point& p = points[pointi];

                // Value against the formula
                maxValueError = max
                (
                    maxValueError,
                    mag(sol.value(p) - formula(dict, a, p))
                );

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
                maxSource = max(maxSource, mag(sol.source(p)));

                // Source against the paper's forcing function, A = (1, 1)
                scalar f = 0;

                if (a == 1 && paperSource(dict, p, f))
                {
                    hasPaperSource = true;
                    maxPaperError = max
                    (
                        maxPaperError,
                        mag(sol.source(p) - f)/(1 + mag(f))
                    );
                }
            }

            nFailed += check("value against formula", maxValueError, valueTol);
            nFailed += check("du/dx against CD", maxGradXError, cdTol);
            nFailed += check("du/dt against CD", maxGradTError, cdTol);
            nFailed += check("du/dz is zero", maxGradZError, 0);
            nFailed += check("source against CD", maxSourceError, cdTol);

            if (hasPaperSource)
            {
                nFailed += check
                (
                    "source against Tufillaro et al. Eq. 125-127,"
                    " |S - f|/(1 + |f|)",
                    maxPaperError,
                    paperTol
                );
            }

            // hasSource() and source() must agree
            if (sol.hasSource())
            {
                Info<< "    hasSource true, max |source| = " << maxSource
                    << endl;

                if (maxSource > zeroSourceTol)
                {
                    Info<< "    PASS: source is non-zero (hasSource true): "
                        << maxSource << " > " << zeroSourceTol << endl;
                }
                else
                {
                    Info<< "    FAIL: source is zero but hasSource is true: "
                        << maxSource << " <= " << zeroSourceTol << endl;
                    nFailed++;
                }
            }
            else
            {
                nFailed += check
                (
                    "source is zero (hasSource false)",
                    maxSource,
                    zeroSourceTol
                );
            }

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
    }

    Info<< nl << nCases << " cases tested" << nl;

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
