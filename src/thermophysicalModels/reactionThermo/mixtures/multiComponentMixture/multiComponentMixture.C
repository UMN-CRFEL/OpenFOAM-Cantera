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

\*---------------------------------------------------------------------------*/

#include "multiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
    }

    return speciesData_[0];
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& thermoData,
    const fvMesh& mesh,
    const word& phaseName,
	Cantera::IdealGasMix& gas 
)
:
    basicSpecieMixture(thermoDict, specieNames, mesh, phaseName),
	gas_(gas),
    speciesData_(species_.size()),
    mixture_("mixture", *thermoData[specieNames[0]]),
    mixtureVol_("volMixture", *thermoData[specieNames[0]]),
    LewisNumber_(species_.size(), 1.0)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*thermoData[species_[i]])
        );
    }
	
	if (thermoDict.found("transportModel"))
	{ 
		const word transportModel(thermoDict.lookup("transportModel"));
        int log = 0;
	    if (transportModel == "mixtureAveraged") 
	    {
	    	MmixtureAveraged_ = true;
			tran_ = Cantera::newTransportMgr("Mix", &gas_, log=0);
	    }
	    if (transportModel == "multiComponent") 
	    {
	    	MmultiComponent_ = true;
			tran_ = Cantera::newTransportMgr("Multi", &gas_, log=0);
	    }
		if (transportModel == "constantLewis")
		{
			MconstantLewis_ = true;
			// the tran_ pointer is still needed for alpha and Cp from Cantera
			tran_ = Cantera::newTransportMgr("Mix", &gas_, log=0);
			// read constant Lewis number of each species 

			dictionary LewisNumberDict_ = thermoDict.subDict("LewisNumber");
			forAll (species_, i)
			{	
			    LewisNumber_[i] = readScalar(LewisNumberDict_.lookup(specieNames[i]));
                if (!LewisNumberDict_.found(specieNames[i])) 
				{
					FatalErrorInFunction
                    << "Lewis number of  " << specieNames[i] 
					<< " not found in input" << abort(FatalError);
				}				
			}
			
		}
		if (transportModel == "unityLewis")
		{
			MunityLewis_ = true;
			tran_ = Cantera::newTransportMgr("Mix", &gas_, log=0);
			forAll (species_, i)
			{				
				LewisNumber_[i] = 1;
			}
		}
	    
	    tranMix_ = dynamic_cast<Cantera::MixTransport*>(tran_);
        tranMulti_ = dynamic_cast<Cantera::MultiTransport*>(tran_);
	}

	Info << "Dezhi in multiComponentMixture: " << gas.nSpecies() << endl;	
//    Info << "Use cantera transport ? " << useCantera << endl;
//    Info << "hehe" << speciesData_[0].W() << endl;	
    correctMassFractions();
}


template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict)),
    mixtureVol_("volMixture", speciesData_[0])
{
	correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Y_[0].boundaryField()[patchi][facei]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv += Y_[i][celli]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/speciesData_[0].rho(p, T)/rhoInv*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/speciesData_[n].rho(p, T)/rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/speciesData_[0].rho(p, T)/rhoInv
      * speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/speciesData_[n].rho(p,T)
          / rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}

template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::muCellMixture
(
    const scalar p,
	const scalar T,
	const label celli
)
{
	scalar muMixture(0.0);
    doublereal Y[species_.size()];
	forAll(species_, i)
	{
		Y[i] = Y_[i][celli];
	}	
	gas_.setState_TPY(T, p, Y);

    if (MmixtureAveraged_)
	{
		muMixture = tranMix_->viscosity();
	}
    if (MmultiComponent_)
	{
		muMixture = tranMulti_->viscosity();
	}		
	
    return muMixture;
}

template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::alphahCellMixture
(
    const scalar p,
	const scalar T,
	const label celli
)
{
	scalar alphahMixture(0.0);
    doublereal Y[species_.size()];
	forAll(species_, i)
	{
		Y[i] = Y_[i][celli];
	}	
	gas_.setState_TPY(T, p, Y);

    if (MmixtureAveraged_)
	{
		alphahMixture = tranMix_->thermalConductivity()/gas_.cp_mass();
	}
    if (MmultiComponent_)
	{
		alphahMixture = tranMulti_->thermalConductivity()/gas_.cp_mass();
	}		
	
    return alphahMixture;
}

template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::muPatchFaceMixture
(
    const scalar p,
	const scalar T,    
	const label patchi,
    const label facei
)
{
	scalar muMixture(0.0);
    doublereal Y[species_.size()];
	forAll(species_, i)
	{
		Y[i] = Y_[i].boundaryField()[patchi][facei];
	}	
	gas_.setState_TPY(T, p, Y);

    if (MmixtureAveraged_)
	{
		muMixture = tranMix_->viscosity();
	}
    if (MmultiComponent_)
	{
		muMixture = tranMulti_->viscosity();
	}		
	
    return muMixture;
}

template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::alphahPatchFaceMixture
(
    const scalar p,
	const scalar T,
	const label patchi,
    const label facei
)
{
	scalar alphahMixture(0.0);
    doublereal Y[species_.size()];
	forAll(species_, i)
	{
		Y[i] = Y_[i].boundaryField()[patchi][facei];
	}	
	gas_.setState_TPY(T, p, Y);

    if (MmixtureAveraged_)
	{
		alphahMixture = tranMix_->thermalConductivity()/gas_.cp_mass();
	}
    if (MmultiComponent_)
	{
		alphahMixture = tranMulti_->thermalConductivity()/gas_.cp_mass();
	}		
	
    return alphahMixture;
}


// ************************************************************************* //
