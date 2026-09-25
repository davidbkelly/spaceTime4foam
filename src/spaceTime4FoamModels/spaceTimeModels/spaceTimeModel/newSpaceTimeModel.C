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

InClass
    spaceTimeModel

\*---------------------------------------------------------------------------*/

#include "spaceTimeModel.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::spaceTimeModel> Foam::spaceTimeModel::New
(
    Time& runTime
)
{
    // NB: dictionary must be unregistered to avoid adding to the database
    word modelType;
    {
        const IOdictionary props
        (
            IOobject
            (
                "spaceTimeProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false  // Do not register
            )
        );

        modelType = props.get<word>("spaceTimeModel");
    }

    Info<< "Selecting spaceTimeModel " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "spaceTimeModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }
#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("spaceTimeModel::New(Time&)")
            << "Unknown spaceTimeModel type " << modelType
            << endl << endl
            << "Valid spaceTimeModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<spaceTimeModel>(ctorPtr(runTime));
}


// ************************************************************************* //
