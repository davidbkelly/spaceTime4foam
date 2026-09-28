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
    Test-spaceTimeLinearExtrapolation

Description
    Test of the spaceTimeLinearExtrapolation boundary condition on a solved
    case (latest time by default, or -time). On every face of every patch
    of u with this condition, with P the owner cell and the gradient
    recomputed here with fvc::grad and the patch's gradSchemeName:

        expected u_f = u_P + grad(u)_P & (C_f - C_P)

    Checks:
    1. Formula (round-off): from the field as read, compute the expected
       value, then update and evaluate the boundary condition once. The new
       patch values must equal the expected values to formulaTol.
    2. Written state: the value entry written by the solver must equal the
       expected value computed from the written field to convergedTol. It
       is not round-off because the solver evaluates the face value after
       the last solve with the gradient of the previous iteration (see the
       boundary condition header); the difference is bounded by the change
       of u in one iteration at convergence.
    3. Read-back: the patch values of the field as read (reconstructed
       from the gradient entry) must equal the written value entry to
       formulaTol.
    4. Not trivial: the written value must differ from the zeroGradient
       value u_P by more than nonTrivialMin on at least one face, so that
       the test cannot pass with zeroGradient behaviour.

    Prints the maximum difference of each check and exits with status 1 if
    any check fails.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "spaceTimeLinearExtrapolationFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// check <description> <maximum difference> <tolerance>
// Print the result and return 1 on failure, 0 on success
label check(const string& description, const scalar maxDiff, const scalar tol)
{
    if (maxDiff <= tol)
    {
        Info<< "    PASS: " << description.c_str() << ": max difference "
            << maxDiff << " <= " << tol << endl;
        return 0;
    }

    Info<< "    FAIL: " << description.c_str() << ": max difference "
        << maxDiff << " > " << tol << endl;
    return 1;
}


//- Expected face values u_P + grad(u)_P & (C_f - C_P) on a patch
tmp<scalarField> expectedValues
(
    const volScalarField& u,
    const volVectorField& gradU,
    const label patchi
)
{
    const fvPatch& p = u.mesh().boundary()[patchi];
    const labelUList& faceCells = p.faceCells();
    const vectorField& Cf = p.Cf();
    const vectorField& C = u.mesh().C().primitiveField();

    tmp<scalarField> tvalues(new scalarField(p.size()));
    scalarField& values = tvalues.ref();

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];

        values[facei] =
            u[celli] + (gradU[celli] & (Cf[facei] - C[celli]));
    }

    return tvalues;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test of the spaceTimeLinearExtrapolation boundary condition on a"
        " solved case (latest time by default)"
    );

    timeSelector::addOptions_singleTime();

    // Only serial runs have been tested
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    // Latest time by default
    if (!timeSelector::setTimeIfPresent(runTime, args))
    {
        const instantList times = runTime.times();
        runTime.setTime(times.last(), times.size() - 1);
    }

#   include "createMesh.H"

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Tolerances. The face values are O(1).
    // - formulaTol: round-off of the same expression evaluated twice
    // - convergedTol: lag of the gradient by one iteration at convergence
    //   (convergenceTolerance 1e-10 on the normalised initial residual)
    // - nonTrivialMin: the extrapolation correction is O(h) on these meshes
    const scalar formulaTol = 1e-12;
    const scalar convergedTol = 1e-8;
    const scalar nonTrivialMin = 1e-3;

    // The file as written, read as a plain dictionary for the value entries
    const dictionary uDict(IFstream(runTime.timePath()/"u")());

    volScalarField u
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Patches with the condition, and their written and read values
    labelList patchIDs;
    List<scalarField> writtenValues;
    List<scalarField> readValues;
    List<word> gradSchemeNames;

    forAll(u.boundaryField(), patchi)
    {
        if
        (
            isA<spaceTimeLinearExtrapolationFvPatchScalarField>
            (
                u.boundaryField()[patchi]
            )
        )
        {
            const spaceTimeLinearExtrapolationFvPatchScalarField& pf =
                refCast<const spaceTimeLinearExtrapolationFvPatchScalarField>
                (
                    u.boundaryField()[patchi]
                );

            const word& patchName = mesh.boundary()[patchi].name();
            const dictionary& patchDict =
                uDict.subDict("boundaryField").subDict(patchName);

            patchIDs.append(patchi);
            writtenValues.append
            (
                scalarField("value", patchDict, pf.size())
            );
            readValues.append(scalarField(pf));
            gradSchemeNames.append(pf.gradSchemeName());
        }
    }

    label nFailed = 0;

    if (patchIDs.empty())
    {
        Info<< "    FAIL: no patch of u has type "
            << spaceTimeLinearExtrapolationFvPatchScalarField::typeName
            << endl;
        nFailed++;
    }

    // Expected values from the field as read (the written state)
    List<scalarField> expected(patchIDs.size());

    forAll(patchIDs, i)
    {
        const tmp<volVectorField> tgradU =
            fvc::grad(u, gradSchemeNames[i]);
        expected[i] = expectedValues(u, tgradU(), patchIDs[i]);
    }

    // Update all the patches first, then evaluate, so that every patch
    // computes its gradient from the same field (as in the solver, where
    // the matrix construction updates all the patches)
    u.boundaryFieldRef().updateCoeffs();

    forAll(patchIDs, i)
    {
        u.boundaryFieldRef()[patchIDs[i]].evaluate();
    }

    forAll(patchIDs, i)
    {
        const label patchi = patchIDs[i];
        const fvPatch& p = mesh.boundary()[patchi];
        const scalarField uP(u.primitiveField(), p.faceCells());
        const scalarField& uUpdated = u.boundaryField()[patchi];

        Info<< nl << "Patch " << p.name() << " (" << p.size() << " faces,"
            << " gradSchemeName " << gradSchemeNames[i] << ", "
            << mesh.gradScheme(gradSchemeNames[i]).toString() << ")" << endl;

        const scalar formulaDiff = gMax(mag(uUpdated - expected[i]));
        const scalar convergedDiff =
            gMax(mag(writtenValues[i] - expected[i]));
        const scalar readBackDiff =
            gMax(mag(readValues[i] - writtenValues[i]));

        const scalarField correction(mag(writtenValues[i] - uP));
        const scalar maxCorrection = gMax(correction);
        label nNonTrivial = 0;

        forAll(correction, facei)
        {
            if (correction[facei] > nonTrivialMin)
            {
                nNonTrivial++;
            }
        }

        nFailed += check
        (
            "1. updated value against u_P + grad(u)_P & (C_f - C_P)",
            formulaDiff,
            formulaTol
        );
        nFailed += check
        (
            "2. written value against u_P + grad(u)_P & (C_f - C_P)",
            convergedDiff,
            convergedTol
        );
        nFailed += check
        (
            "3. read-back value against written value",
            readBackDiff,
            formulaTol
        );

        if (nNonTrivial > 0)
        {
            Info<< "    PASS: 4. ";
        }
        else
        {
            Info<< "    FAIL: 4. ";
            nFailed++;
        }
        Info<< nNonTrivial << " of " << p.size() << " faces have"
            << " |u_f - u_P| > " << nonTrivialMin << " (max |u_f - u_P| = "
            << maxCorrection << ")" << endl;
    }

    Info<< nl;

    if (nFailed > 0)
    {
        Info<< "Test-spaceTimeLinearExtrapolation: " << nFailed
            << " check(s) FAILED" << nl << endl;
        return 1;
    }

    Info<< "Test-spaceTimeLinearExtrapolation: all checks passed" << nl
        << endl;

    return 0;
}


// ************************************************************************* //
