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

Class
    Foam::waterWaveModels::BichromaticWave

Description
    Class of biochromatic waves

Usage
    in inputDict
    WWParameters
    {
        WaveType            BichromaticWave;
        WaterDepth          0.7;
        WaveHeightI         0.06;
        WavePeriodI         3.0;
        WaveHeightII        0.04;
        WavePeriodII        2.0;
    }

SourceFiles
    BichromaticWave.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef BICHROMATICWAVE_H_INCLUDED
#define BICHROMATICWAVE_H_INCLUDED
#include "WaterWaves.H"
#include "SecondOrderPiston.H"
#include "WaveComponent.H"

namespace Foam
{
namespace waterWaveModels
{
class BichromaticWave
:
    public waterWaves
{

    scalar waveHeightI_;
    scalar wavePeriodI_;

    scalar waveHeightII_;
    scalar wavePeriodII_;

    waveComponent wave1_;
    waveComponent wave2_;

    bool SecondOrderCorrection_;

    scalar enhanceI_;
    scalar enhanceII_;

    scalar TDamp_;

    void Cal_Enhancefactor();

public:
    TypeName("BichromaticWave");

    // Constructor
    BichromaticWave
    (
        const IOdictionary& dict
    );


    // Return the name of target wave
    word name();

    // surface elevation
    scalar eta(scalar t);

    // phase velocity
    scalar Cal_C();

    // Calculate the displacement of piston-type wavemaker
    scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

    void PrintWaveProperties();

}; // end of class
} // end of namespace waterwaveModels
} // end of namespace foam
#endif
