/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::foamChemistryReader

Description
    Chemistry reader for OpenFOAM format

SourceFiles
    foamChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef foamChemistryReader_H
#define foamChemistryReader_H

#include "chemistryReader.H"
#include "fileName.H"
#include "typeInfo.H"
#include "HashPtrTable.H"
#include "labelList.H"
#include "speciesTable.H"
#include "atomicWeights.H"
#include "cantera/IdealGasMix.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class foamChemistry Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class foamChemistryReader
:
    public chemistryReader<ThermoType>
{
        Cantera::IdealGasMix gas_;
		
		//- Chemistry/reactions dictionary
        dictionary chemDict_;

        //- Thermo properties dictionary
        dictionary thermoDict_;

        //- List of elements
        DynamicList<word> elementNames_;

        //- Element indices
        HashTable<label> elementIndices_;

        //- Table of species
        speciesTable& speciesTable_;

        //- Table of species composition
        speciesCompositionTable speciesComposition_;

        //- Table of the thermodynamic data given in the foamChemistry file
        HashPtrTable<ThermoType> speciesThermo_;

        //- List of the reactions
        ReactionList<ThermoType> reactions_;


    // Private Member Functions

        //- Set the species list
        speciesTable& setSpecies(const dictionary& dict, speciesTable& species);

        //- Read the species composition
        void readSpeciesComposition();

        //- Disallow default bitwise copy construct
        foamChemistryReader(const foamChemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const foamChemistryReader&);


public:

    //- Runtime type information
    TypeName("foamChemistryReader");


    // Constructors

        //- Construct from foamChemistry and thermodynamics file names
        foamChemistryReader
        (
            const fileName& reactionsFileName,
            speciesTable& species,
            const fileName& thermoFileName
        );

        //- Construct by getting the foamChemistry and thermodynamics file names
        //  from dictionary
        foamChemistryReader
        (
            const dictionary& thermoDict,
            speciesTable& species
        );


    //- Destructor
    virtual ~foamChemistryReader()
    {}


    // Member functions

		Cantera::IdealGasMix& canteraGas()
		{
			return gas_;
		}
		
        //- List of elements
        const wordList& elementNames() const
        {
            return elementNames_;
        }

        //- Element indices
        const HashTable<label>& elementIndices() const
        {
            return elementIndices_;
        }

        //- Table of species
        const speciesTable& species() const
        {
            return speciesTable_;
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }

        //- Table of the thermodynamic data given in the foamChemistry file
        const HashPtrTable<ThermoType>& speciesThermo() const
        {
            return speciesThermo_;
        }

        //- List of the reactions
        const ReactionList<ThermoType>& reactions() const
        {
            return reactions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamChemistryReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
