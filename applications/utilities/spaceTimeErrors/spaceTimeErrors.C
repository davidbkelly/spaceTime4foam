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
    spaceTimeErrors

Description
    Error norms of the space-time solution u against the analytical
    solution selected in constant/spaceTimeProperties. The latest time
    (iteration) is used, unless -time or -latestTime is given.

    The branch is chosen by the spaceTimeModel keyword. Only cellCentred is
    implemented so far; vertexCentred will be added with that model.

    cellCentred:
    - Error at each cell centroid: e = u_P - u_exact(x_P).
    - Weight: the cell area A_P = V_P / Lz, where Lz is the z extent of the
      mesh (Lz = 1 on the tutorial meshes, so A_P = V_P).
    - Norms over a subset S of cells:
          L1   = sum_S |e| A / sum_S A
          L2   = sqrt(sum_S e^2 A / sum_S A)
          Linf = max_S |e|
      for S = all cells, interior cells and boundary cells (cells with at
      least one face on a non-empty, non-coupled patch).
    - Final-time error on the tEnd patch (t = T), weighted by the face
      length |S_f| / Lz, with the same norm definitions, for two face
      values:
          (a) extrapolated: u_P + grad(u)_P & (x_f - x_P), with grad(u)
              from fvc::grad using the grad(u) scheme of the case;
          (b) cell value u_P, which is the face value that zeroGradient
              gives.
    - Mesh size h = sqrt(A_domain/(nTriangles/2)), with A_domain the area of
      the front patch (the frontAndBack faces at the minimum z). This is
      1/N on the N x N triangle meshes. If the front faces are not all
      triangles, h = sqrt(A_domain/nFrontFaces) is used instead (1/N on an
      N x N quad mesh) and a note is printed.

    Output: postProcessing/spaceTimeErrors/errors.dat (one data line; the
    commented header names the columns) and a summary on screen.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "emptyPolyPatch.H"
#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- L1, L2 and Linf norms of an error with weights
class errorNorms
{
public:

    scalar L1;
    scalar L2;
    scalar Linf;
    scalar weightSum;
    label n;

    errorNorms()
    :
        L1(0),
        L2(0),
        Linf(0),
        weightSum(0),
        n(0)
    {}

    //- Add one error value with its weight
    void add(const scalar e, const scalar w)
    {
        L1 += mag(e)*w;
        L2 += sqr(e)*w;
        Linf = max(Linf, mag(e));
        weightSum += w;
        n++;
    }

    //- Normalise by the total weight (call once, after all add calls)
    void finalise()
    {
        reduce(L1, sumOp<scalar>());
        reduce(L2, sumOp<scalar>());
        reduce(Linf, maxOp<scalar>());
        reduce(weightSum, sumOp<scalar>());
        reduce(n, sumOp<label>());

        if (weightSum > 0)
        {
            L1 = L1/weightSum;
            L2 = Foam::sqrt(L2/weightSum);
        }
    }
};


//- Write the three norms to a stream, separated by spaces
Ostream& writeNorms(Ostream& os, const errorNorms& norms)
{
    os  << norms.L1 << " " << norms.L2 << " " << norms.Linf;
    return os;
}


//- Print one line of the screen summary
void printNorms(const word& name, const errorNorms& norms)
{
    Info<< "    " << name << " (n = " << norms.n << "):"
        << " L1 = " << norms.L1
        << "  L2 = " << norms.L2
        << "  Linf = " << norms.Linf << endl;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Error norms of u against the analytical solution in"
        " constant/spaceTimeProperties (latest time by default)"
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

    // Registered, so that spaceTimeAnalyticalFixedValue uses it
    IOdictionary spaceTimeProperties
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

    const word method(spaceTimeProperties.get<word>("spaceTimeModel"));

    autoPtr<analyticalSolution> analytical = analyticalSolution::New
    (
        spaceTimeProperties.subDict("analyticalSolution"),
        analyticalSolution::spaceTimeVelocity(spaceTimeProperties)
    );


    // Mesh size h from the front patch (the frontAndBack faces at min z)

    const scalar zMin = mesh.bounds().min().z();
    const scalar Lz = mesh.bounds().max().z() - zMin;

    const label frontAndBackID =
        mesh.boundaryMesh().findPatchID("frontAndBack");

    if (frontAndBackID < 0)
    {
        FatalErrorInFunction
            << "Cannot find the frontAndBack patch" << exit(FatalError);
    }

    const polyPatch& frontAndBack = mesh.boundaryMesh()[frontAndBackID];
    const vectorField fbCf(frontAndBack.faceCentres());
    const vectorField fbSf(frontAndBack.faceAreas());

    scalar domainArea = 0;
    label nFrontFaces = 0;
    label nFrontTriangles = 0;

    forAll(frontAndBack, facei)
    {
        if (mag(fbCf[facei].z() - zMin) < 1e-6*Lz)
        {
            domainArea += mag(fbSf[facei]);
            nFrontFaces++;

            if (frontAndBack[facei].size() == 3)
            {
                nFrontTriangles++;
            }
        }
    }

    scalar h = 0;

    if (nFrontFaces > 0 && nFrontTriangles == nFrontFaces)
    {
        h = Foam::sqrt(domainArea/(0.5*nFrontTriangles));
    }
    else if (nFrontFaces > 0)
    {
        h = Foam::sqrt(domainArea/nFrontFaces);

        Info<< "Note: the front faces are not all triangles: using"
            << " h = sqrt(A_domain/nFrontFaces)" << nl << endl;
    }
    else
    {
        FatalErrorInFunction
            << "No frontAndBack faces found at z = " << zMin
            << exit(FatalError);
    }


    // Errors

    errorNorms allNorms;
    errorNorms interiorNorms;
    errorNorms boundaryNorms;
    errorNorms tEndExtrapolatedNorms;
    errorNorms tEndCellNorms;
    label nUnknowns = 0;
    const scalar timeValue = runTime.value();

    if (method == "cellCentred")
    {
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

        // Mark the boundary cells: cells with a face on a non-empty,
        // non-coupled patch
        boolList isBoundaryCell(mesh.nCells(), false);

        forAll(mesh.boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (!isA<emptyPolyPatch>(pp) && !pp.coupled())
            {
                const labelUList& faceCells = pp.faceCells();

                forAll(faceCells, facei)
                {
                    isBoundaryCell[faceCells[facei]] = true;
                }
            }
        }

        const vectorField& C = mesh.C().primitiveField();
        const scalarField& V = mesh.V().field();
        const scalarField uExact(analytical().value(C));

        forAll(C, celli)
        {
            const scalar e = u[celli] - uExact[celli];
            const scalar area = V[celli]/Lz;

            allNorms.add(e, area);

            if (isBoundaryCell[celli])
            {
                boundaryNorms.add(e, area);
            }
            else
            {
                interiorNorms.add(e, area);
            }
        }

        nUnknowns = mesh.nCells();

        // Final-time error on the tEnd patch
        const label tEndID = mesh.boundaryMesh().findPatchID("tEnd");

        if (tEndID < 0)
        {
            FatalErrorInFunction
                << "Cannot find the tEnd patch" << exit(FatalError);
        }

        // Uses the grad(u) scheme of system/fvSchemes
        const volVectorField gradU(fvc::grad(u));

        const fvPatch& tEndPatch = mesh.boundary()[tEndID];
        const labelUList& faceCells = tEndPatch.faceCells();
        const vectorField& Cf = tEndPatch.Cf();
        const scalarField magSf(tEndPatch.magSf());
        const scalarField uExactFace(analytical().value(Cf));

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            const scalar length = magSf[facei]/Lz;

            const scalar uExtrapolated =
                u[celli] + (gradU[celli] & (Cf[facei] - C[celli]));

            tEndExtrapolatedNorms.add
            (
                uExtrapolated - uExactFace[facei],
                length
            );

            tEndCellNorms.add(u[celli] - uExactFace[facei], length);
        }
    }
    else
    {
        FatalErrorInFunction
            << "spaceTimeErrors is not implemented for spaceTimeModel "
            << method << " yet" << exit(FatalError);
    }

    allNorms.finalise();
    interiorNorms.finalise();
    boundaryNorms.finalise();
    tEndExtrapolatedNorms.finalise();
    tEndCellNorms.finalise();


    // Summary

    Info<< "Method " << method << ", h = " << h
        << ", nUnknowns = " << nUnknowns
        << ", domain area = " << domainArea << nl
        << "Errors e = u - u_exact at the unknowns, area weighted:" << endl;
    printNorms("all", allNorms);
    printNorms("interior", interiorNorms);
    printNorms("boundary", boundaryNorms);
    Info<< "Errors on tEnd (t = T), length weighted:" << endl;
    printNorms("extrapolated", tEndExtrapolatedNorms);
    printNorms("cell value", tEndCellNorms);
    Info<< endl;


    // Write the data file

    const fileName outputDir
    (
        runTime.globalPath()/"postProcessing"/"spaceTimeErrors"
    );
    mkDir(outputDir);

    OFstream os(outputDir/"errors.dat");
    os.precision(10);

    os  << "# spaceTimeErrors: errors of u against the analytical solution"
        << nl
        << "# time (iteration) " << timeValue << nl
        << "# Columns:" << nl
        << "#  1 method" << nl
        << "#  2 h = sqrt(A_domain/(nTriangles/2))" << nl
        << "#  3 nUnknowns" << nl
        << "#  4 nInterior" << nl
        << "#  5 nBoundary" << nl
        << "#  6- 8 L1 L2 Linf over all unknowns" << nl
        << "#  9-11 L1 L2 Linf over interior unknowns" << nl
        << "# 12-14 L1 L2 Linf over boundary unknowns" << nl
        << "# 15-17 L1 L2 Linf on tEnd, extrapolated face value" << nl
        << "# 18-20 L1 L2 Linf on tEnd, cell value" << nl;

    os  << method << " " << h << " " << nUnknowns << " "
        << interiorNorms.n << " " << boundaryNorms.n << " ";
    writeNorms(os, allNorms) << " ";
    writeNorms(os, interiorNorms) << " ";
    writeNorms(os, boundaryNorms) << " ";
    writeNorms(os, tEndExtrapolatedNorms) << " ";
    writeNorms(os, tEndCellNorms) << nl;

    Info<< "Written " << os.name() << nl << nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
