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
    Foam::waterWaveModels::SuperimposeWave::FocusWave

Description
    Class of focused waves.

Usage
    in inputDict
    WWParameters
    {
        WaveType                FocusWave;
        WaterDepth              0.7;
        FocusAmplitude          0.14;
        Fre_start               1.0;
        Fre_end                 2.0;
        FocusTime               30.0;
        FocusPosition           30.0;
        AmplitudeDiscretization CWA;//CWS, NewWaveJS, NewWavePM
        NComponents             30;
    }

SourceFiles
    FocusWave.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef FocusWave_H
#define FocusWave_H
#include "SecondOrderPiston.H"
#include "WaveComponent.H"
#include "SuperimposeWave.H"

namespace Foam
{
namespace waterWaveModels
{
class FocusWave
:
    virtual public SuperimposeWave
{
protected:

    // focus wave amplitude
    const scalar amplitude_;

    // Start and end of the frequency range
    const scalar f_start_;
    const scalar f_end_;

    // peak frequency
    const scalar f_peak_;

    // Focus time and position
    const scalar t_focus_;
    const scalar x_focus_;

    // approach to disperse amplitude
    // CWA: constant wave amplitude
    // CWS: constant wave steepness
    // NewWave: according to spectrum (JONSWAP)
    const word AmplitudeDiscretization_;

    // arrays for transfer function
    List<List<complex>> TF_final_super_;
    List<List<complex>> TF_final_sub_;

    // store all the transfer functions
    void storeTFs();

public:
    TypeName("FocusWave");
    FocusWave
    (
        const IOdictionary& dict
    );

    // Destructor
    virtual ~FocusWave(){};

    //- Member Functions
    void updateCoeff();

    // member functions
      // calculate arrays
    void Cal_TF_array();
    void Cal_LkC_array();
    void Cal_A_array();
    void Cal_Tf_array();

      // Return the name of the target wave
    word name();

      // surface elevation
    scalar eta(scalar t);

      // surface elevation, for active absorption
    scalar eta1st(scalar t);

      // phase velocity
    scalar Cal_C();

      // Calculate the displacement of piston-type wavemaker
    scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

	// active absorption
	scalar Cal_X_corr(scalar t, scalar deltaT, scalar WaterLevel);

    void PrintWaveProperties();
};
} // end of namespace waterWaveModels
} // end of namespace Foam
#endif
