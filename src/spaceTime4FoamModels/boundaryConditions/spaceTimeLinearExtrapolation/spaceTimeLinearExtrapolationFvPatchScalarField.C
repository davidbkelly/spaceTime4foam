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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    gradSchemeName_(defaultGradSchemeName())
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
    )
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
    gradSchemeName_(ptf.gradSchemeName_)
{}


Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const spaceTimeLinearExtrapolationFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    gradSchemeName_(ptf.gradSchemeName_)
{}


Foam::spaceTimeLinearExtrapolationFvPatchScalarField::
spaceTimeLinearExtrapolationFvPatchScalarField
(
    const spaceTimeLinearExtrapolationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    gradSchemeName_(ptf.gradSchemeName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
    fvPatchScalarField::writeValueEntry(os);
#else
    if (gradSchemeName_ != defaultGradSchemeName())
    {
        os.writeEntry("gradSchemeName", gradSchemeName_);
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
