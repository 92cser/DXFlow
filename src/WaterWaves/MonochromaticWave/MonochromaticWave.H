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
    Foam::waterWaveModels::monochromaticWave

Description
    Class of monochromati, linear waves

Usage
    in inputDict
    WWParameters
    {
        WaveType           monochromaticWave;
        WaterDepth         0.7;
        WaveHeight         0.14;
        WavePeriod         3.0;

        // optional
        // reference:
        // Hemming A. Schäffer, G. Klopman. (2000). Review of multidirectional
            active wave absorption methods. Journal of Waterway, Port, Coastal,
            and Ocean Engineering. 126(2), 88-97.
        ActiveAbsorption		true;

        // optional
        // reference:
        // Hemming A., Schäffer. (1996). Second-order wavemaker theory for
           irregular waves. Ocean Engineering. 23(1), 47-88.
        SecondOrderCorrection		true;
    }

SourceFiles
    monochromaticWave.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/
#ifndef MonochromaticWave_H
#define MonochromaticWave_H
#include "WaterWaves.H"
#include "SecondOrderPiston.H"
#include "WaveComponent.H"

namespace Foam
{
namespace waterWaveModels
{
class monochromaticWave
: 
    public waterWaves
{
protected:
    const scalar waveHeight_;
    const scalar wavePeriod_;
    const scalar waveLength_;
    const scalar k_;
    const scalar Omega_;
    const scalar stroke_;

    bool activeAbsorption_;
    // Displacement correction of active absorption
    scalar X_corr_;

    // second-order correction
    bool SecondOrderCorrection_;
	
public:

    TypeName("monochromaticWave");

    // Constructor
    monochromaticWave
    (
        const IOdictionary& dict
    );

// Destructor
    virtual ~monochromaticWave(){};

    // Member functions

    // wave length
    scalar Cal_L();

    // Stroke for piston-type wavemaker
    scalar Stroke();

    // Calculate the displacement correction of wave paddle
    scalar Cal_X_corr(scalar t, scalar deltaT, scalar WaterLevel);

    // Return the name of target wave
    word name();

    // surface elevation
    scalar eta(scalar t);

    // phase velocity
    scalar Cal_C();

    // Calculate the displacement of piston-type wavemaker
    scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

    void PrintWaveProperties();

};//end of class 

}//end of namespace waterWaveModels
}//end of namespace Foam
#endif
