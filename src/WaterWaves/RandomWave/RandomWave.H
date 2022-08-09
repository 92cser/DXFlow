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
    Foam::waterWaveModels::SuperimposeWave::RandomWave

Description
    Class of random waves.

Usage
    in inputDict
    WWParameters
    {
        WaveType                    RandomWave;
        SpectrumType                Jonswap; //PM
        NComponents                 100;
        WaterDepth                  0.4;
        SignificantWaveHeight       0.053;
        PeakPeriod                  1.48;
        AmplitudeDiscretization     AF; //AE
    }

SourceFiles
    RandomWave.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/
#ifndef RandomWave_H
#define RandomWave_H
#include "SecondOrderPiston.H"
#include "WaveComponent.H"
#include "SuperimposeWave.H"
#include "Random.H"

namespace Foam
{
namespace waterWaveModels
{
class RandomWave 
:
    virtual public SuperimposeWave
{
	// whether a JONSWAP or P-M spectrum is adopted
	word SpectrumType_;

    // significant wave height
    const scalar waveHeight_sig_;

    // peak period
    const scalar peakPeriod_;

    // peak frequency
    const scalar peakFrequency_;

    // elevation factor
    const scalar gamma_;

    // Start and end of the frequency range
    const scalar f_start_;
    const scalar f_end_;

    // approach to disperse amplitude
    // AE: averged energy
    // AF: averged freqency
    const word AmplitudeDiscretization_;

    // array for random phases
    List<scalar> RandomPhase_;

    // switch of active absorption
    bool activeAbsorption_;

    // displacement correction of active absorption
    scalar X_corr_;

    // member functions
      // calculate arrays
    void Cal_TF_array();
    void Cal_LkC_array();
    void Cal_A_array
    (
        scalar& deltaFreqency,
        scalar& deltaEnergy
    );
    void Cal_Tf_array();
    // Calculate the random phases
    void Cal_randomPhase();

    // calculate enhanced factor
    void Cal_EF_array();

    // Calculate all the arrays
    void updateCoeff();

    // arrays for transfer function
    List<List<complex>> TF_final_super_;
    List<List<complex>> TF_final_sub_;

    // store all the transfer functions
    void storeTFs();
	
public:
    TypeName("RandomWave");

    // Constructor
    RandomWave
    (
            const IOdictionary& dict
    );

    // Destructor
    virtual ~RandomWave(){};

    // Return the name of the target wave
    word name();

    // surface elevation
    scalar eta(scalar t);

    // phase velocity
    scalar Cal_C();

    // Calculate the displacement of piston-type wavemaker
    scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

    void PrintWaveProperties();

}; //end of class RandomWave 
} // end of namespace waterWaveModels
} // end of namespace Foam
#endif
