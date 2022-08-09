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
    Foam::waterWaveModels::solitary

Description
    Class of solitary waves.

    Reference:
    Fenton, J. D. (1972) A ninth-order solution for the solitary wave, Journal
    of Fluid Mechanics 53, 257-271.

Usage
    in inputDict
    WWParameters
    {
        WaveType           solitary;
        WaterDepth         0.5;
        SolitaryAmplitude  0.2;
    }

SourceFiles
    Solitary.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef Solitary_H
#define Solitary_H
#include "WaterWaves.H"

namespace Foam
{
namespace waterWaveModels
{
class Solitary
:
	public waterWaves
{
protected:
	// Solitary wave amplitude
	const scalar amplitude_;
	
	// h / d
	const scalar nonlinearity_;
	
	// phase velocity 
	const scalar phaseVelocity_;
	
	// wave number of solitary wave
	const scalar k_;
	
	// Inital Position Of the Wave Crest
	const scalar IPOWC_;
	
	// Position of piston-type wave paddle 
	scalar position_;
	
public:
	TypeName("solitary");

	// Constructor
	Solitary
	(	
		const IOdictionary& dict
	);

    // Destructor
    virtual ~Solitary(){};
	
	// Calculate the wave number
	scalar Cal_k(); 

	// Return the name of target wave 
	word name();
	
	// surface elevation
	scalar eta(scalar t);

	// calculate phase velocity
	scalar Cal_C_Solitary();

	// return the phase velocity
	scalar Cal_C();
	
	// Calculate the displacement of piston-type wavemaker
	scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

	void PrintWaveProperties();
};


}// end of namespace waterWaveModels
}// end of namespace Foam




#endif
