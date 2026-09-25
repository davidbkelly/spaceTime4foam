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
    analyticalSolution

\*---------------------------------------------------------------------------*/

#include "analyticalSolution.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::analyticalSolution> Foam::analyticalSolution::New
(
    const dictionary& dict,
    const vector& A
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting analyticalSolution " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "analyticalSolution",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn("analyticalSolution::New(...)", dict)
            << "Unknown analyticalSolution type " << modelType
            << endl << endl
            << "Valid analyticalSolution types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<analyticalSolution>(ctorPtr(dict, A));
}


// ************************************************************************* //
