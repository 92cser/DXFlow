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
    Foam::waterWaveModels::NoWave

Description
    No waves.

Usage
    in inputDict
    WWParameters
    {
        WaveType        	noWave;
        WaterDepth		0.0;
    }

SourceFiles
    NoWave.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef NoWave_H
#define NoWave_H
#include "WaterWaves.H"

namespace Foam
{
namespace waterWaveModels
{
class NoWave
:	
	public waterWaves
{

public:
	TypeName("noWave");

	// Constructor
	NoWave
	(	
		const IOdictionary& dict
	);

	// Destructor
	virtual ~NoWave(){};

	// Member functions

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
} // end of namespace waterWaveModels
} // end of namespace foam
#endif
