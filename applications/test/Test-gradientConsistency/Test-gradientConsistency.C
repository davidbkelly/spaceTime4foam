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
    Test-gradientConsistency

Description
    Gradient-consistency diagnostic (CLAUDE.md section 6, S5 finding):
    the cell gradient of an exactly linear field with the OpenFOAM
    gradient schemes used by the cellCentred benchmark.

    Sets the cell field u = x + 2 t (x = mesh x, t = mesh y), with the
    exact face values u(x_f) on every non-empty boundary face, so the
    exact gradient is (1, 2, 0) everywhere. Evaluates the gradient with
        GaussLinear   Gauss linear
        leastSquares  leastSquares
    (selected by name with fv::gradScheme::New, as fvc::grad does from
    gradSchemes), and prints, for each scheme and for the interior cells
    (no face on a non-empty boundary patch) and all cells, the maximum
    and the RMS (unweighted mean over the cells) of |grad(u) - (1, 2, 0)|.

    Output: one line per scheme and subset,
        gradientError <scheme> <subset> <nCells> <max> <rms>
    Nothing is written to the case.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "emptyPolyPatch.H"
#include "gradScheme.H"
#include "IStringStream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Print the error of the Gauss linear and leastSquares cell gradients"
        " of the linear field u = x + 2 t"
    );

    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const vector exactGrad(1, 2, 0);

    volScalarField u
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    );

    const volVectorField& C = mesh.C();
    u.primitiveFieldRef() = C.component(vector::X)().primitiveField()
        + 2*C.component(vector::Y)().primitiveField();

    // Exact values on the boundary faces; interior cells are those with no
    // face on a non-empty boundary patch
    boolList interior(mesh.nCells(), true);

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if (isA<emptyPolyPatch>(patch.patch()))
        {
            continue;
        }

        const vectorField& Cf = patch.Cf();
        u.boundaryFieldRef()[patchi] ==
            Cf.component(vector::X)() + 2*Cf.component(vector::Y)();

        const labelUList& faceCells = patch.faceCells();
        forAll(faceCells, i)
        {
            interior[faceCells[i]] = false;
        }
    }

    const wordList schemeTags({"GaussLinear", "leastSquares"});
    const wordList schemeNames({"Gauss linear", "leastSquares"});

    Info<< "u = x + 2 t, exact gradient " << exactGrad << nl
        << "Columns: gradientError scheme subset nCells max rms" << nl;

    forAll(schemeTags, schemei)
    {
        IStringStream schemeData(schemeNames[schemei]);
        tmp<fv::gradScheme<scalar>> scheme
        (
            fv::gradScheme<scalar>::New(mesh, schemeData)
        );

        const volVectorField gradU
        (
            scheme().grad(u, "grad(u)" + schemeTags[schemei])
        );

        const scalarField err(mag(gradU.primitiveField() - exactGrad));

        for (const word subset : {"interior", "all"})
        {
            label n = 0;
            scalar maxErr = 0;
            scalar sumSqr = 0;

            forAll(err, celli)
            {
                if (subset == "all" || interior[celli])
                {
                    n++;
                    maxErr = max(maxErr, err[celli]);
                    sumSqr += sqr(err[celli]);
                }
            }

            const scalar rms = (n > 0 ? Foam::sqrt(sumSqr/n) : 0);

            Info<< "gradientError " << schemeTags[schemei] << " " << subset
                << " " << n << " " << setprecision(10) << maxErr << " "
                << rms << setprecision(6) << nl;
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
