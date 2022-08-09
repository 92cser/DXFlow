/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 Dalian Ocean University 
						    National Marine Environment Monitoring Center 
							Dalian University of Technology
    Copyright (C) 2018-2022 Ocean University of China
	Copyright (C) 2022-     Ningbo University 
-------------------------------------------------------------------------------
License
    This file is part of DXFlow, a toolbox developed based on OpenFOAM.

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
#include "StokesII.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(StokesII, 0);
    addToRunTimeSelectionTable(waterWaves, StokesII, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

// Constructor
Foam::waterWaveModels::StokesII::StokesII
(	
    const IOdictionary& dict
)
:
    monochromaticWave(dict),
    etaStokesII(0.0)
{}

void Foam::waterWaveModels::StokesII::Cal_SOItem(scalar t)
{
    scalar CH = cosh(k_*waterDepth_);
    scalar SH = sinh(k_*waterDepth_);
    scalar setDown = -0.5*pow(0.5*waveHeight_,2.0)*k_/SH;
    scalar lockphaseWave = 0.25*pow(0.5*waveHeight_,2.0)*k_*CH*(2*CH*CH + 1)/pow(SH,3)*cos(2*Omega_*t);
    etaStokesII = setDown + lockphaseWave;
}

Foam::word Foam::waterWaveModels::StokesII::name()
{
    return "StokesII wave";
}

Foam::scalar Foam::waterWaveModels::StokesII::eta(scalar t)
{
    Cal_SOItem(t);
    // slowly start factor
    scalar SSF_ = (wavePeriod_ > 0.5*t) ? pow(0.5*t/wavePeriod_,2.0) : 1.0;
    return (monochromaticWave::eta(t) + etaStokesII)*SSF_;
}

Foam::scalar Foam::waterWaveModels::StokesII::Cal_C()
{
    return monochromaticWave::Cal_C();
}

Foam::scalar Foam::waterWaveModels::StokesII::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
    return monochromaticWave::displacement(t, deltaT, WaterLevel);
}

void Foam::waterWaveModels::StokesII::PrintWaveProperties()
{
    Info << "---------------------------------" << "\n"
         << "Wave properties of StokesII wave" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Wave Height (m) = " << waveHeight_ << "\n"
         << " Wave Period (s) = " << wavePeriod_ << "\n"
         << " Wave Length (m) = " << waveLength_ << "\n"
         << "---------------------------------" << endl;
}
