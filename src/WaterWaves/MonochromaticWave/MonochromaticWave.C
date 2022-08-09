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
#include "MonochromaticWave.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(monochromaticWave, 0);
    addToRunTimeSelectionTable(waterWaves, monochromaticWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

// Constructor
Foam::waterWaveModels::monochromaticWave::monochromaticWave
(	
    const IOdictionary& dict
)
:
    waterWaves(dict),
    waveHeight_(readScalar(WWParameters_.lookup("WaveHeight"))),
    wavePeriod_(readScalar(WWParameters_.lookup("WavePeriod"))),
    Omega_(2*pi_/wavePeriod_),
    waveLength_(Cal_L()),
    k_(2*pi_/waveLength_),
    stroke_(Stroke()),
    activeAbsorption_(WWParameters_.lookupOrDefault<bool>("ActiveAbsorption",false)),
    SecondOrderCorrection_(WWParameters_.lookupOrDefault<bool>("SecondOrderCorrection",false)),
    X_corr_(0.0)
{}

Foam::word Foam::waterWaveModels::monochromaticWave::name()
{
    return "Monochromatic wave";
}

Foam::scalar Foam::waterWaveModels::monochromaticWave::eta(scalar t)
{
    // slowly start factor
    scalar SSF_ = (wavePeriod_ > 0.5*t) ? pow(0.5*t/wavePeriod_,2.0) : 1.0;
    return 0.5*waveHeight_*cos(Omega_*t)*SSF_;
}

Foam::scalar Foam::waterWaveModels::monochromaticWave::Cal_L()
{
    scalar mu0_ = pow(2*pi_,2)*waterDepth_/g_*pow(wavePeriod_,-2);
    scalar index0_ = 1.835+1.225*pow(mu0_,1.35);
    scalar mu_ = mu0_*(1+mu0_*exp(-index0_))*pow(tanh(mu0_),-0.5);
    // another iteration to improve precise
    scalar muFinal_ = (pow(mu_,2)+mu0_*pow(cosh(mu_),2))/(mu_+0.5*sinh(2*mu_));
    return 2*pi_*waterDepth_/muFinal_;
}
	
Foam::scalar Foam::waterWaveModels::monochromaticWave::Cal_C()
{
    return waveLength_/wavePeriod_;
}

Foam::scalar Foam::waterWaveModels::monochromaticWave::Stroke()
{
    return 0.25*waveHeight_*(2*k_*waterDepth_+sinh(2*k_*waterDepth_))/sinh(k_*waterDepth_)/sinh(k_*waterDepth_);
}

Foam::scalar Foam::waterWaveModels::monochromaticWave::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
    scalar secondOrderCorrection = 0.0;
    if(SecondOrderCorrection_)
    {
        waveComponent wave1(waveHeight_,waterDepth_,wavePeriod_);
        waveComponent wave2(waveHeight_,waterDepth_,wavePeriod_);
        SecondOrderWave CombinedWave(wave1, wave2);
        secondOrderCorrection = CombinedWave.X2_super(t) + CombinedWave.X2_sub(t);
    }
    return 0.5*stroke_*sin(Omega_*t) + secondOrderCorrection - Cal_X_corr(t, deltaT, WaterLevel);
}
	
Foam::scalar Foam::waterWaveModels::monochromaticWave::Cal_X_corr(scalar t, scalar deltaT, scalar WaterLevel)
{
    if(activeAbsorption_)
    {
        if(t > 3.0*wavePeriod_)
        {
            scalar correctWaterLevel = WaterLevel- waterDepth_ - eta(t);
            scalar U_corr_ = Omega_*0.25*(2*k_*waterDepth_+sinh(2*k_*waterDepth_))/sinh(k_*waterDepth_)/sinh(k_*waterDepth_)* correctWaterLevel;
            X_corr_ += U_corr_*deltaT;
        }
    }
    return X_corr_;
}

void Foam::waterWaveModels::monochromaticWave::PrintWaveProperties()
{
    Info << "-----------------------------------------------" << "\n"
         << "Wave properties of monochromatic wave" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Wave Height (m) = " << waveHeight_ << "\n"
         << " Wave Period (s) = " << wavePeriod_ << "\n"
         << " Wave Length (m) = " << waveLength_ << "\n"
         << "-----------------------------------------------" << endl;
    if(activeAbsorption_)
    {
            Info << "Active absorption will be adopted after the third wave period" << endl;
    }
}
