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
    Foam::hePsiThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    hePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef canteraPsiThermo_H
#define canteraPsiThermo_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class canteraPsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class canteraPsiThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        canteraPsiThermo(const canteraPsiThermo<BasicPsiThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("canteraPsiThermo");


    // Constructors

        //- Construct from mesh and phase name
        canteraPsiThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~canteraPsiThermo();


    // Member functions

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "canteraPsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
