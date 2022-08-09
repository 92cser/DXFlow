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
    Foam::waterWaveModels::waterWaves

Description
    Basic class of water waves.

SourceFiles
    WaterWaves.C
    newWaterWaves.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/
#ifndef waterWaves_H
#define waterWaves_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
namespace waterWaveModels
{
class waterWaves
{
protected:
	// Sub-dictionary of water wave parameters
	const dictionary WWParameters_;

	// Water depth
	const scalar waterDepth_;

public:
	// Runtime type information
    TypeName("waterWaves");
	
	// Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        waterWaves,
        nameList,
        (const IOdictionary& dict),
        (dict)
    );
	
	// Constructor
	waterWaves
	(	
		const IOdictionary& dict
	);
	
	// Selector
    static autoPtr<waterWaves> New
    (
		const IOdictionary& dict
    );
	
    // Destructor
    virtual ~waterWaves(){};
	
	// Return the name of the target wave
	virtual word name() = 0;
	
	// Access for mass source wavemaker
	// Calculate surface elevations at each time step
	virtual scalar eta(scalar t) = 0;
	
	// Calculate phase velocity
	virtual scalar Cal_C() = 0;
	
	// Access for piston-type wavemaker
	// Calculate the displacement of piston-type wavemaker
	virtual scalar displacement(scalar t, scalar deltaT, scalar WaterLevel) = 0;

	// print wave properties
	virtual void PrintWaveProperties() = 0;

	// Constants
	const scalar g_ = 9.81;
	const scalar pi_ = 3.141592653589793;
	const scalar Vsmall_ = 0.001;
	
};
}// end of namespace waterWaveModels
}// end of namespace Foam
#endif
