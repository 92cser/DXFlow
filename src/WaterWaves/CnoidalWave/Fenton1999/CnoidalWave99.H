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
    Foam::waterWaveModels::CnoidalWave

Description
    Class of Cnoidal waves (a modification to Fenton 1979)

    Reference:
    Fenton, J. D. (1999) Numerical methods for nonlinear waves, in P. L.-F. Liu, ed.,
    Advances in Coastal and Ocean Engineering, Vol. 5, World Scientific, Singapore, pp. 241-324.

Usage
    in inputDict
    WWParameters
    {
        WaveType           CnoidalWave;
        WaterDepth         0.7;
        WaveHeight         0.05;
        WavePeriod         3.0;
    }

SourceFiles
    CnoidalWave99.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef CnoidalWave99_H
#define CnoidalWave99_H
#include "WaterWaves.H"

namespace Foam
{
namespace waterWaveModels
{
class CnoidalWave
:
	public waterWaves
{
    scalar waveHeight_;
    scalar waveLength_;
    scalar wavePeriod_;
    scalar phaseVelocity_;

    // trough depth
    scalar h_;
    // wave height /waterDepth_
    scalar Epsl_;
    // modulus
    scalar m_;
    // position of the wave paddle
    scalar position_;

    //functions
    // Complete elliptic integral of the first kind: K(m)
    scalar K();
    // Complementary complete elliptic integral of the first kind K'(m)
    scalar K0();
    // Complete elliptic integral of the second kind: K(m)
    scalar E();
    // Ratio of elliptic integrals, e = E(m)/K(m)
    scalar e();
    // Complementary nome
    scalar q1();
    // Trough depth /waterDepth_
    scalar hbyD();
    // Wave height / trough depth
    scalar Hbyh();
    // Coefficient
    scalar Alpha();
    // non-dimensional phase velocity, C*(gh)^(-0.5)
    scalar NDPhaseVelocity();
    // solve the equation for calculating the m_, with the use of biosection method
    void solve();
    // the tolerance of the equation solved
    scalar tolerance(double m0);
    // update coefficients after solving m_
    void updateCoeffs();
    // wavelength
    scalar wavelength();
    // jacobian elliptic function
    scalar cn(double z);

public:
        TypeName("CnoidalWave");
	
	// Constructor
        CnoidalWave
	(	
            const IOdictionary& dict
	);

	// Destructor
        virtual ~CnoidalWave(){};

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
