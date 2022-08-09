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
#include "BichromaticWave.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(BichromaticWave, 0);
    addToRunTimeSelectionTable(waterWaves, BichromaticWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

// Constructor
Foam::waterWaveModels::BichromaticWave::BichromaticWave
(
    const IOdictionary& dict
)
:
    waterWaves(dict),
    waveHeightI_(readScalar(WWParameters_.lookup("WaveHeightI"))),
    wavePeriodI_(readScalar(WWParameters_.lookup("WavePeriodI"))),
    waveHeightII_(readScalar(WWParameters_.lookup("WaveHeightII"))),
    wavePeriodII_(readScalar(WWParameters_.lookup("WavePeriodII"))),
    wave1_(waveHeightI_,waterDepth_,wavePeriodI_),
    wave2_(waveHeightII_,waterDepth_,wavePeriodII_),
    SecondOrderCorrection_(WWParameters_.lookupOrDefault<bool>("SecondOrderCorrection", false)),
    enhanceI_(1.0),
    enhanceII_(1.0),
    TDamp_(max(wavePeriodI_, wavePeriodII_))
{
    Cal_Enhancefactor();
}

void Foam::waterWaveModels::BichromaticWave::Cal_Enhancefactor()
{
    { // calculate enhanceI_
        scalar relativeWaterDepth = wave1_.k_[0].Re()*waterDepth_/2.0/3.1415926;
        scalar nonlinearity = waveHeightI_/waterDepth_;

        // theta value
        scalar theta = (1 - (1 - exp(4*nonlinearity) * nonlinearity)*relativeWaterDepth)*relativeWaterDepth;

        scalar thetaLimit = 0.5/(1 - exp(4*nonlinearity) * nonlinearity);

        // enhanced factor
        enhanceI_ = 183.88*pow3(theta) - 63.26*pow(theta,2) + 7.58*theta + 0.72;

        // the minimum value is 1.0, maximum value is 3.115
        enhanceI_ = min(max(1.0, enhanceI_), 3.115);

        if(relativeWaterDepth > thetaLimit)// deep water waves, theta decreases with the increase of D/L
        {
            enhanceI_ = 3.115;
        }
    }
    { // calculate enhanceII_
        scalar relativeWaterDepth = wave2_.k_[0].Re()*waterDepth_/2.0/3.1415926;
        scalar nonlinearity = waveHeightII_/waterDepth_;

        // theta value
        scalar theta = (1 - (1 - exp(4*nonlinearity) * nonlinearity)*relativeWaterDepth)*relativeWaterDepth;

        scalar thetaLimit = 0.5/(1 - exp(4*nonlinearity) * nonlinearity);

        // enhanced factor
        enhanceII_ = 183.88*pow3(theta) - 63.26*pow(theta,2) + 7.58*theta + 0.72;

        // the minimum value is 1.0, maximum value is 3.115
        enhanceII_ = min(max(1.0, enhanceII_), 3.115);

        if(relativeWaterDepth > thetaLimit)// deep water waves, theta decreases with the increase of D/L
        {
            enhanceII_ = 3.115;
        }
    }
    //Info << "Enhanced Factors are " << enhanceI_ << " and " << enhanceII_ << endl;
}

Foam::word Foam::waterWaveModels::BichromaticWave::name()
{
    return "BichromaticWave wave";
}

Foam::scalar Foam::waterWaveModels::BichromaticWave::eta(scalar t)
{
    scalar waveSurface1 = enhanceI_*wave1_.A(t).Re() * wave1_.C();
    scalar waveSurface2 = enhanceII_*wave2_.A(t).Re() * wave2_.C();

    // slowly start factor
    scalar SSF_ = (t < TDamp_) ? pow(t/TDamp_,2.0) : 1.0;
    return SSF_*(waveSurface1+waveSurface2);
}

Foam::scalar Foam::waterWaveModels::BichromaticWave::Cal_C()
{
    return 1.0;
}

Foam::scalar Foam::waterWaveModels::BichromaticWave::displacement(scalar t, scalar deltaT, scalar correctWaterLevel)
{
    scalar Dis_wave1 =  wave1_.A(t).Im()/wave1_.TF_[0].Re();
    scalar Dis_wave2 =  wave2_.A(t).Im()/wave2_.TF_[0].Re();

    scalar secondOrderCorrection = 0.0;
    if(SecondOrderCorrection_)
    {
        SecondOrderWave CombinedWave(wave1_, wave2_);
        secondOrderCorrection = CombinedWave.X2_super(t) + CombinedWave.X2_sub(t);
    }
    // slowly start factor
    scalar SSF_ = (t < TDamp_) ? pow(t/TDamp_,2.0) : 1.0;
    return SSF_*(Dis_wave1 + Dis_wave2 + secondOrderCorrection);
}

void Foam::waterWaveModels::BichromaticWave::PrintWaveProperties()
{
    Info << "-----------------------------------------------" << "\n"
         << "Wave properties of biochromatic wave" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Wave Height (m) = " << waveHeightI_ << ", " << waveHeightII_ << "\n"
         << " Wave Period (s) = " << wavePeriodI_ << ", " << wavePeriodII_ << "\n"
         << " Enhanced factors = " << enhanceI_ << ", " << enhanceII_ << "\n"
         << " Relative Water Depth (kd) = " << wave1_.k_[0].Re()*waterDepth_ << ", " << wave2_.k_[0].Re()*waterDepth_ << "\n"
         << "-----------------------------------------------" << endl;
}
