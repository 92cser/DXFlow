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
#include "WaterWaves.H"
namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(waterWaves, 0);
    defineRunTimeSelectionTable(waterWaves, nameList);
}
}

// Constructor
Foam::waterWaveModels::waterWaves::waterWaves
(
	const IOdictionary& dict
)
:
    WWParameters_(dict.subDict("WWParameters")),
    waterDepth_(readScalar(WWParameters_.lookup("WaterDepth")))
{}

