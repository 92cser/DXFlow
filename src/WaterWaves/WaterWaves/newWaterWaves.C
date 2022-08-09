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
Foam::autoPtr<waterWaves> waterWaves::New
(
	const IOdictionary& dict
)
{
        const word waveType_ = dict.subDict("WWParameters").get<word>("WaveType");
	nameListConstructorTable::iterator cstrIter
        = nameListConstructorTablePtr_->find(waveType_);
	
	if (cstrIter == nameListConstructorTablePtr_->end())
  	{
        FatalErrorInFunction
            << "Unknown " << waterWaves::typeName << " " << waveType_ << nl << nl
            << "Valid types are:" << nl
            << nameListConstructorTablePtr_->sortedToc()
            << exit(FatalError);
        }
	return cstrIter()(dict);
}
}
}
