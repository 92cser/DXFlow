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
    Foam::waterWaveModels::SuperimposeWave

Description
    A midst class between waterWaveModels and FocusWave or RandomWave

SourceFiles
    SuperimposeWave.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef SuperimposeWave_H
#define SuperimposeWave_H

#include "WaterWaves.H"

namespace Foam
{
namespace waterWaveModels
{
class SuperimposeWave
:
    public waterWaves
{
protected:
    // Number of components
    const int NComponents_;
	// active absorprion or not 
	bool activeAbsorption_;
	// correction paddle displacement
	scalar X_corr_;
    // array for frequency
    List<scalar> F_array_;
    // array for period
    List<scalar> T_array_;
    // array for wave length
    List<scalar> L_array_;
    // array for amplitude (for piston-type wavemaker)
    List<scalar> A_array_;
    // array for amplitude (after enhanced, for mass source wavemaker)
    List<scalar> A_EF_array_;
    // array for phase velocity
    List<scalar> C_array_;
    // array for wave number
    List<scalar> k_array_;
    // array for enhance factors
    List<scalar> Ef_array_;
    // array for transfer function
    List<scalar> Tf_array_;
    // array for wave spectrum
    List<scalar> WS_;

    // wave paddle position
    scalar position_;
    // second-order correction
    bool SecondOrderCorrection_;

    // JONSWAP
    scalar beta(scalar gamma);
    scalar sigma(scalar f_, scalar f_peak);
    scalar para_Jon(scalar f_, scalar f_peak);
    scalar Jonswap_(scalar Hs, scalar Tp, scalar gamma, scalar f_);
	
	// PM
	scalar PM_(scalar f_peak, scalar f_);

public:
    TypeName("SuperimposeWave");

    // Constructor
    SuperimposeWave
    (
        const IOdictionary& dict
    );

    // Destructor
    virtual ~SuperimposeWave(){};

    // Member Functions

    // Return the name of the target wave
    // Because this class is the base class, the error function is used to stop running
    word name()
    {
            FatalErrorInFunction
            << "SuperimposeWave is a base class that cannot directly used."
            << exit(FatalError);
    };

    // surface elevation
    scalar eta(scalar t){return 0;};

    // phase velocity
    scalar Cal_C(){return 0;};

    // Calculate the displacement of piston-type wavemaker
    scalar displacement(scalar t, scalar deltaT, scalar WaterLevel){return 0;};

    void PrintWaveProperties(){};
}; //end of class SuperimposeWave 
} // end of namespace waterWaveModels
} // end of namespace Foam

#endif
