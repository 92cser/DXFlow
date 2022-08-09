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
#include "NoWave.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
	defineTypeNameAndDebug(NoWave, 0);
    addToRunTimeSelectionTable(waterWaves, NoWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

// Constructor
Foam::waterWaveModels::NoWave::NoWave
(	
	const IOdictionary& dict
)
:
	waterWaves(dict)
{}

Foam::word Foam::waterWaveModels::NoWave::name()
{
	return "No wave";
}

Foam::scalar Foam::waterWaveModels::NoWave::eta(scalar t)
{
	return 0.0;
}

	
Foam::scalar Foam::waterWaveModels::NoWave::Cal_C()
{
	return 0.0;
}

Foam::scalar Foam::waterWaveModels::NoWave::displacement(scalar t, scalar deltaT, scalar correctWaterLevel)
{
	return 0.0;
}

void Foam::waterWaveModels::NoWave::PrintWaveProperties()
{
    Info << "-----------------------------------------------" << "\n"
         << "            No Wave model is used              " << "\n"
		 << "-----------------------------------------------" << endl;
}
