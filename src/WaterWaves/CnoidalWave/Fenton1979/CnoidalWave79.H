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
    Foam::waterWaveModels::CnoidalWaveS

Description
    Class of Cnoidal waves (small amplitude only)

    Reference:
    Fenton, J. D. (1979) A high-order cnoidal wave theory, Journal of
    Fluid Mechanics 94, 129-161.

Usage
    in inputDict
    WWParameters
    {
        WaveType           CnoidalWaveS;
        WaterDepth         0.7;
        WaveHeight         0.05;
        UseT               true;
        WavePeriod         3.0;
    }
    or
    WWParameters
    {
        WaveType           CnoidalWaveS;
        WaterDepth         0.7;
        WaveHeight         0.05;
        UseT               false;
        WaveLength         5.0;
    }
SourceFiles
    CnoidalWave79.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef CnoidalWave79_H
#define CnoidalWave79_H
#include "WaterWaves.H"

namespace Foam
{
namespace waterWaveModels
{
class CnoidalWaveS
:
	public waterWaves
{
	// Wave properties
	scalar waveHeight_;
	scalar wavePeriod_;
	scalar waveLength_;
	scalar phaseVelocity_;
	// Epsl = waveHeight_/waterDepth_
	const scalar Epsl_;

	// Parameters, need to solve
	// actual module number
	scalar m_;
	// approximate module number 
	scalar mApproximate_;
	// straining factor
	scalar alpha_;
	// minimum water depth
	scalar troughWaterDepth_;
	// Epsl0 = waveHeight_/troughWaterDepth_
	scalar Epsl0_;

	// position of wavemaker
	scalar position_;

	// true: use H,D,T calculate L (default)
	// false: use H,D,L calculate T
	bool UseT;

	// the complete elliptic integral of first kind
	scalar K(double m);
	// the complete elliptic integral of second kind
	scalar E(double m);

	// Coefficients
	// wave length coefficients
	scalar LCoeffs_[5][5][5];
	// wave period coefficients
	scalar TCoeffs_[5][5][5];
	// trough water depth coefficients
	scalar TWDCoeffs_[6][6][6];
	// straining factor coefficients
	scalar AlphaCoeffs_[5][5];
	// phase velocity coefficients
	scalar PVCoeffs_[6][6][2];
	// surface elevation coefficients
	scalar EtaCoeffs_[6][6][6];

	// Initialization
	void InitializationAll();

	void InitializeWaveLengthCoeff_();
	void InitializeWavePeriodCoeff_();
	void InitializeTWDCoeff_();
	void InitializeAlphaCoeff_();
	void InitializePVCoeff_();
	void InitializeEtaCoeff_();

	// Calculation
	void CalculationAll();
	// calculate module number
	void Cal_m();
	// calculate trough water depth, Epsl0 calculation
	void Cal_TWD();
	// calculate straining factor
	void Cal_alpha();
	// calculate phase velocity
	void Cal_PV();

	// Arithmetic-Geometric Mean method, three initial values and time
	scalar AGM_Cn
	(
	    scalar a0, 
	    scalar b0, 
	    scalar c0, 
	    scalar t
	);

	// Access for mass source wavemaker
	// Calculate surface elevations at each time step
	scalar etaForPiston(scalar t);

public:
        TypeName("CnoidalWaveS");
	
	// Constructor
        CnoidalWaveS
	(	
		const IOdictionary& dict
	);

	// Destructor
        virtual ~CnoidalWaveS(){};

	// Return the name of the target wave
	word name();
	
	// Access for mass source wavemaker
	// Calculate surface elevations at each time step
	scalar eta(scalar t);
	
	// Return calculated phase velocity
	scalar Cal_C();
	
	// Access for piston-type wavemaker
	// Calculate the displacement of piston-type wavemaker
	scalar displacement(scalar t, scalar deltaT, scalar WaterLevel);

	// print wave properties
	void PrintWaveProperties();

}; // end of class CnoidalWave
} // end of namespace waterWaveModels
} // end of namespace foam
#endif
