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
    Foam::waterWaveModels::StreamFunctionWave

Description
    Class of monochromatic waves, calculated using Fourier approximation.

    Reference:
    Fenton, J. D. (1999) Numerical methods for nonlinear waves, in P. L.-F. Liu, ed.,
    Advances in Coastal and Ocean Engineering, Vol. 5, World Scientific, Singapore, pp. 241-324.

Usage
    in inputDict
    WWParameters
    {
        WaveType           StreamFunctionWave;
        WaterDepth         0.5;
        Ej
        (
                4.6796992e-02
                3.3881428e-03
                2.4702329e-04
                2.0289131e-05
                1.8311399e-06
                1.7581812e-07
        );
        Properties
        (
                6.6777 // nondimensional wavelength
                0.1    // nondimensional wave height
                7.5301 // nondimensional wave period
                0.8868 // nondimensional wave speed
        );
    }

    Input Ej and Properties are calculated using Fenton's latest program, see:
    http://www.johndfenton.com/Steady-waves/Fourier.html.

SourceFiles
    StreamFunction.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef StreamFunction_H
#define StreamFunction_H

#include "WaterWaves.H"

namespace Foam
{
namespace waterWaveModels
{
class StreamFunctionWave
:
        public waterWaves
{

    // coefficients, read from dictionary
    // calculated by Fenton's latest program:
    // http://www.johndfenton.com/Steady-waves/Fourier.html.

    // List<scalar> Bj_;  // Bj is not used if we don't calculate velocity
    List<scalar> Ej_;

    // first four elements are used (wave length, wave height, wave period, wave speed)
    List<scalar> Properties_;

    // wave properties
    // wave height
    scalar H_;
    // wave length
    scalar L_;
    // wave number
    scalar k_;
    // wave period
    scalar T_;
    // wave speed
    scalar C_;
    // these parameters are not used
    /*// Eulerian current
    scalar EC_;
    // Stokes current
    scalar SC_;
    // mean fluid speed in frame of wave
    scalar MS_;
    // discharge
    scalar Q_;
    // BournBernoulli constant
    scalar R_;*/
    // frequency to prevent slow drift
    scalar OmegaC_;
    // transfer function
    scalar ts_;

    // position of wavemaker(for shallow water condition)
    scalar position_;

    void updateWaveProperties();

    scalar Cal_L(scalar D, scalar T);

    // surface elevation for piston
    scalar etaForPiston(scalar x, scalar t);

public:

    TypeName("StreamFunctionWave");

    // Constructor
    StreamFunctionWave
    (
        const IOdictionary& dict
    );

// Destructor
    ~StreamFunctionWave(){};

    // Return the name of target wave
    word name();

    // surface elevation
    scalar eta(scalar t);

    // phase velocity
    scalar Cal_C();

    // Calculate the displacement of piston-type wavemaker
    scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

    void PrintWaveProperties();

};
}
}
#endif
