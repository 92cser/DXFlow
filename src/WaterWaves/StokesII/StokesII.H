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
    Foam::waterWaveModels::monochromaticWave::StokesII

Description
    Class of monochromatic, StokesII waves.

Usage
    in inputDict
    WWParameters
    {
        WaveType           StokesII;
        WaterDepth         0.7;
        WaveHeight         0.14;
        WavePeriod         3.0;
    }

SourceFiles
    StokesII.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef StokesII_H
#define StokesII_H

#include "MonochromaticWave.H"

namespace Foam
{
namespace waterWaveModels
{
class StokesII
:
	public monochromaticWave
{
	// second-order item based on the airy wave 
	scalar etaStokesII;
	
public:
	TypeName("StokesII");
	
	// Constructor
	StokesII
	(	
		const IOdictionary& dict
	);
	
	// Destructor
        virtual ~StokesII(){};

	// Member Functions
	
	// Calculate the second-order item 
	void Cal_SOItem(scalar t);

	// Return the name of the target wave
	word name();
	
	// surface elevation
	scalar eta(scalar t);
	
	// phase velocity
	scalar Cal_C();
	
	// Calculate the displacement of piston-type wavemaker
	scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

	void PrintWaveProperties();

}; //end of class  StokesII 
} // end of namespace waterWaveModels
} // end of namespace Foam
#endif

