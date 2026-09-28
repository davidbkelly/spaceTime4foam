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

\*---------------------------------------------------------------------------*/

#include "spaceTimeLinearExtrapolationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word
Foam::spaceTimeLinearExtrapolationFvPatchScalarField::defaultGradSchemeName()
const
{
    return "grad(" + internalField().name() + ')';
}


bool
Foam::spaceTimeLinearExtrapolationFvPatchScalarField::gaussTypeGradient()
const
{
    // Allow-list: the first word is Gauss or iterativeGauss, or the first
    // word is a limited scheme and the second word is Gauss or
    // iterativeGauss. Anything else (leastSquares, fourth, a scheme from a
    // user library, ...) is not Gauss-type.
    const ITstream& is = internalField().mesh().gradScheme(gradSchemeName_);

    if (is.size() < 1 || !is[0].isWord())
    {
        return false;
    }

    const word& first = is[0].wordToken();

    if (first == "Gauss" || first == "iterativeGauss")
    {
        return true;
    }

    const bool limited =
    (
        first == "cellLimited"
     || first.starts_with("cellLimited<")
     || first == "cellMDLimited"
     || first == "faceLimited"
     || first == "faceMDLimited"
    );

    if (limited && is.size() >= 2 && is[1].isWord())
    {
        const word& second = is[1].wordToken();

        return (second == "Gauss" || second == "iterativeGauss");
    }

    return false;
}


const Foam::boolList&
Foam::spaceTimeLinearExtrapolationFvPatchScalarField::multiFaces() const
{
    if (!multiFacePtr_)
    {
        const volScalarField& u =
            refCast<const volScalarField>(internalField());

        // Number of faces with this condition in each cell, over all patches
        labelList nFaces(u.mesh().nCells(), Zero);

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
                const labelUList& faceCells =
                    u.mesh().boundary()[patchi].faceCells();

                forAll(faceCells, facei)
                {
                    nFaces[faceCells[facei]]++;
                }
            }
        }

        multiFacePtr_.reset(new boolList(size(), false));
        boolList& multi = multiFacePtr_.ref();

        const labelUList& faceCells = patch().faceCells();
        label nMulti = 0;

        forAll(faceCells, facei)
        {
            if (nFaces[faceCells[facei]] >= 2)
            {
                multi[facei] = true;
                nMulti++;
            }
        }

        reduce(nMulti, sumOp<label>());

        if (cornerFallback_)
        {
            Info<< type() << ": patch " << patch().name() << ": " << nMulti
                << " of " << returnReduce(size(), sumOp<label>())
                << " faces use the corner fallback (u_f = u_P)" << endl;
        }
        else if (nMulti > 0)
        {
            if (!gaussTypeGradient())
            {
                FatalErrorInFunction
                    << "Patch " << patch().name() << " of field "
                    << internalField().name() << ": " << nMulti
                    << " faces are in cells with two or more "
                    << type() << " faces, cornerFallback is off and"
                    << " the gradient scheme " << gradSchemeName_ << " ("
                    << internalField().mesh().gradScheme(gradSchemeName_)
                      .toString()
                    << ") is not a Gauss-type gradient scheme." << nl
                    << "    In such a cell the extrapolated face values"
                    << " feed one gradient component back into itself."
                    << " Without the fallback only Gauss-type schemes"
                    << " (Gauss, iterativeGauss, or a limited scheme"
                    << " wrapping one of them) are allowed: they give a"
                    << " neutral mode (gain exactly 1 for Gauss linear)."
                    << " With leastSquares the gain"
                    << " is 1.5 on the left-diagonal triangle meshes, so"
                    << " the face values diverge, possibly without any sign"
                    << " in the residual; other schemes have not been"
                    << " analysed." << nl
                    << "    Set cornerFallback on (the default)."
                    << exit(FatalError);
            }
            else
            {
                WarningInFunction
                    << "Patch " << patch().name() << " of field "
                    << internalField().name() << ": " << nMulti
                    << " faces are in cells with two or more "
                    << type() << " faces and cornerFallback is off." << nl
                    << "    In such a cell one gradient component is fed"
                    << " back into itself. With the Gauss-type scheme "
                    << gradSchemeName_ << " ("
                    << internalField().mesh().gradScheme(gradSchemeName_)
                      .toString()
                    << ") this is a neutral mode (gain exactly 1 for"
                    << " Gauss linear; iterativeGauss and cellLimited Gauss"
                    << " were also seen to be neutral), so the solution is"
                    << " not unique and depends on the starting field." << nl
                    << "    Set cornerFallback on (the default) for a"
                    << " unique solution." << endl;
            }
        }
    }

    return multiFacePtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    gradSchemeName_(defaultGradSchemeName()),
    cornerFallback_(true),
    multiFacePtr_()
{}


Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
#ifdef OPENFOAM_COM
    // The gradient entry is optional: if it is missing, the face value is
    // set to the cell value (zero gradient) and the gradient to zero
    fixedGradientFvPatchScalarField(p, iF, dict, IOobjectOption::LAZY_READ),
#else
    fixedGradientFvPatchScalarField(p, iF),
#endif
    gradSchemeName_
    (
        dict.getOrDefault<word>("gradSchemeName", defaultGradSchemeName())
    ),
    cornerFallback_(dict.getOrDefault<Switch>("cornerFallback", true)),
    multiFacePtr_()
{
#ifndef OPENFOAM_COM
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        evaluate();
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
        gradient() = 0;
    }
#endif
}


Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const spaceTimeLinearExtrapolationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    gradSchemeName_(ptf.gradSchemeName_),
    cornerFallback_(ptf.cornerFallback_),
    multiFacePtr_()
{}


Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const spaceTimeLinearExtrapolationFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    gradSchemeName_(ptf.gradSchemeName_),
    cornerFallback_(ptf.cornerFallback_),
    multiFacePtr_()
{}


Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const spaceTimeLinearExtrapolationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    gradSchemeName_(ptf.gradSchemeName_),
    cornerFallback_(ptf.cornerFallback_),
    multiFacePtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::spaceTimeLinearExtrapolationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    multiFacePtr_.clear();
}


void Foam::spaceTimeLinearExtrapolationFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
    multiFacePtr_.clear();
}


void Foam::spaceTimeLinearExtrapolationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volScalarField& u = refCast<const volScalarField>(internalField());

    // Gradient of the current field with the gradSchemes entry
    // gradSchemeName_ (reused if cached in system/fvSolution)
    const tmp<volVectorField> tgradU = fvc::grad(u, gradSchemeName_);

    const vectorField gradUP
    (
        patch().patchInternalField(tgradU().primitiveField())
    );

    // Full offset vector from the owner cell centre to the face centre
    const vectorField delta(patch().Cf() - patch().Cn());

    // fixedGradient sets u_f = u_P + gradient/deltaCoeffs: multiplying by
    // deltaCoeffs gives u_f = u_P + (grad(u)_P & delta)
    gradient() = patch().deltaCoeffs()*(gradUP & delta);

    // Corner fallback: zeroGradient value in cells with two or more faces
    // with this condition
    if (cornerFallback_)
    {
        const boolList& multi = multiFaces();

        forAll(multi, facei)
        {
            if (multi[facei])
            {
                gradient()[facei] = 0;
            }
        }
    }
    else
    {
        // Checks the configuration and stops or warns
        multiFaces();
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::spaceTimeLinearExtrapolationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);

#ifdef OPENFOAM_COM
    os.writeEntryIfDifferent<word>
    (
        "gradSchemeName",
        defaultGradSchemeName(),
        gradSchemeName_
    );
    os.writeEntryIfDifferent<Switch>("cornerFallback", true, cornerFallback_);
    fvPatchScalarField::writeValueEntry(os);
#else
    if (gradSchemeName_ != defaultGradSchemeName())
    {
        os.writeEntry("gradSchemeName", gradSchemeName_);
    }
    if (!cornerFallback_)
    {
        os.writeEntry("cornerFallback", cornerFallback_);
    }
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        spaceTimeLinearExtrapolationFvPatchScalarField
    );
}

// ************************************************************************* //
