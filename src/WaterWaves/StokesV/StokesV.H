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
    Foam::waterWaveModels::monochromaticWave::StokesV

Description
    Class of monochromatic, StokesV waves.

    Reference:
    Fenton, J.D. (1985) A fifth-order Stokes theory for steady waves,
    Journal of Waterway, Port, Coastal and Ocean Engineering. 111, 216-234.

Usage
    in inputDict
    WWParameters
    {
        WaveType           StokesV;
        WaterDepth         0.7;
        WaveHeight         0.14;
        WavePeriod         3.0;
    }

SourceFiles
    StokesV.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef StokesV_H
#define StokesV_H

#include "MonochromaticWave.H"

namespace Foam
{
namespace waterWaveModels
{
class StokesV
:
	public monochromaticWave
{
	// Wave length for StokesV waves
	scalar waveLength_V_;
	
	// Parameters
	scalar B_22;
	scalar B_24;
	scalar B_33;
	scalar B_35;
	scalar B_44;
	scalar B_55;
	scalar C1;
	scalar C2;
	
public:
	TypeName("StokesV");
	
	// Constructor
	// from dictionary, run time selector
	StokesV
	(	
		const IOdictionary& dict
	);

	// Destructor
    virtual ~StokesV(){};

	// Member Functions
	scalar K();    
	scalar CH();	//	cosh(kd)
	scalar SH();	//	sinh(kd)
	scalar K1();	//	K1 =  H/d
	scalar K2();    //  2*pi*d/g/T/T
	scalar Q();     //  Q = K2 / [(d/L)tanh(kd)]
	scalar Lambda(); 
	scalar d_L();   //  d/L
	void Cal_B_22();
	void Cal_B_24();
	void Cal_B_33();
	void Cal_B_35();
	void Cal_B_44();
	void Cal_B_55();
	void Cal_C1();
	void Cal_C2();
	void Cal_All();
	void Solve();

	// Return the name of the target wave
	word name();
	
	// surface elevation
	scalar eta(scalar t);
	
	// phase velocity
	scalar Cal_C();
	
	// Calculate the displacement of piston-type wavemaker
	scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

	void PrintWaveProperties();

}; //end of class  StokesV 

} // end of namespace waterWaveModels
} // end of namespace Foam
#endif

