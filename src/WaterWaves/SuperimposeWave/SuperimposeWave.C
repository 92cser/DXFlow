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
#include "SuperimposeWave.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(SuperimposeWave, 0);
    addToRunTimeSelectionTable(waterWaves, SuperimposeWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

// Constructor
Foam::waterWaveModels::SuperimposeWave::SuperimposeWave
(	
    const IOdictionary& dict
)
:
    waterWaves(dict),
    NComponents_(readScalar(WWParameters_.lookup("NComponents"))),
    activeAbsorption_(WWParameters_.lookupOrDefault<bool>("ActiveAbsorption", false)),
    X_corr_(0.0),
    F_array_(NComponents_, 0.0),
    T_array_(F_array_),
    L_array_(T_array_),
    A_array_(T_array_),
    A_EF_array_(T_array_),
    C_array_(T_array_),
    k_array_(T_array_),
    Ef_array_(T_array_),
    Tf_array_(T_array_),
    WS_(T_array_),
    position_(0.0),
    SecondOrderCorrection_(WWParameters_.lookupOrDefault<bool>("SecondOrderCorrection", false))
{}

// Jonswap spectrum
Foam::scalar Foam::waterWaveModels::SuperimposeWave::beta(scalar gamma) 
{
    return 0.06238*(1.094-0.01915*log(gamma))/(0.23+0.0336*gamma-0.185/(1.9+gamma));
}

Foam::scalar Foam::waterWaveModels::SuperimposeWave::sigma(scalar f_, scalar f_peak)
{
    return (f_ > f_peak) ? 0.09: 0.07;
}

Foam::scalar Foam::waterWaveModels::SuperimposeWave::para_Jon(scalar f_, scalar f_peak)
{
    return exp(-0.5*pow((f_/f_peak-1)/sigma(f_,f_peak),2.0));
}

Foam::scalar Foam::waterWaveModels::SuperimposeWave::Jonswap_(scalar Hs, scalar Tp, scalar gamma, scalar f_)
{
    return beta(gamma)*pow(Hs,2.0)*pow(Tp,-4.0)*pow(f_,-5.0)*exp(-1.25*pow(Tp*f_,-4.0))*pow(gamma,para_Jon(f_,1/Tp));
}

// P-M spectrum, use f rather than omega
Foam::scalar Foam::waterWaveModels::SuperimposeWave::PM_(scalar f_peak, scalar f_)
{
    return 0.0005*pow(f_, -5)*exp(-1.25*pow(f_peak/f_, 4));
}

