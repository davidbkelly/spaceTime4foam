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

#include "spaceTimeAnalyticalFixedValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::analyticalSolution&
Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::solution() const
{
    if (!analyticalSolutionPtr_)
    {
        const Time& runTime = db().time();
        const word dictName("spaceTimeProperties");

        if (runTime.foundObject<IOdictionary>(dictName))
        {
            // Registered by the spaceTime4Foam solver
            const dictionary& spaceTimeProperties =
                runTime.lookupObject<IOdictionary>(dictName);

            analyticalSolutionPtr_ = analyticalSolution::New
            (
                spaceTimeProperties.subDict("analyticalSolution"),
                analyticalSolution::spaceTimeVelocity(spaceTimeProperties)
            );
        }
        else
        {
            // Not registered, e.g. in a utility: read the file
            const IOdictionary spaceTimeProperties
            (
                IOobject
                (
                    dictName,
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false  // Do not register
                )
            );

            analyticalSolutionPtr_ = analyticalSolution::New
            (
                spaceTimeProperties.subDict("analyticalSolution"),
                analyticalSolution::spaceTimeVelocity(spaceTimeProperties)
            );
        }
    }

    return analyticalSolutionPtr_();
}


Foam::tmp<Foam::scalarField>
Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::analyticalValues() const
{
    return solution().value(patch().Cf());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::
spaceTimeAnalyticalFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    analyticalSolutionPtr_()
{}


Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::
spaceTimeAnalyticalFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
#ifdef OPENFOAM_COM
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
#else
    fixedValueFvPatchScalarField(p, iF, dict, false),
#endif
    analyticalSolutionPtr_()
{
    fvPatchScalarField::operator=(analyticalValues());
}


Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::
spaceTimeAnalyticalFixedValueFvPatchScalarField
(
    const spaceTimeAnalyticalFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    analyticalSolutionPtr_()
{}


Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::
spaceTimeAnalyticalFixedValueFvPatchScalarField
(
    const spaceTimeAnalyticalFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    analyticalSolutionPtr_()
{}


Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::
spaceTimeAnalyticalFixedValueFvPatchScalarField
(
    const spaceTimeAnalyticalFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    analyticalSolutionPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::spaceTimeAnalyticalFixedValueFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchScalarField::operator==(analyticalValues());

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        spaceTimeAnalyticalFixedValueFvPatchScalarField
    );
}

// ************************************************************************* //
