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
    Foam::massSourceWaveMaker

Description
    A class of mass source wavemaker. A line-shaped source region is recommended,
    and the velocity field inside source regions can be adjusted to avoid the
    undesired vortical phenomenon near the source region.

    Reference:
    Wang, D., Dong, S. (2021). Generating shallow-and intermediate-water waves
    using a line-shaped mass source wavemaker. Ocean Engineering, 220, 108493.

    Note: The mass source can also be uneven distributed. But it is not recommended
    because the effect of using an uneven mass source seems not obvious.

Usage
    MassSourceWaveMaker
    {
        AreaOfSource        		0.0255;
        EnhanceFactor       		1.1;
        Gamma				1e-6;
        Zt				-0.075;
        Zb				-0.5;
        AdjustU
        {
                AdjustUndesiredW    	true;
                AdjustHeight	    	1.0;
                IndexOfAdjustment	0.1;
                AdjustmentScheme	tri; // power, tri, exp
                EnsureSymmetry		false;
        }
    }

SourceFiles
    MassSource.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef MassSource_H
#define MassSource_H

#include "WaterWaves.H"
#include "volFields.H"

namespace Foam
{
class massSourceWaveMaker
{
    //- Mark field
    //  0: outside source region; 1: inside source region
    volScalarField& markField_;

    //- Access of mesh
    const fvMesh& mesh_;

    //- Wave source value at each time step
    volScalarField& waveSourceValue_;

    //- Velocity value (for adjustment)
    volVectorField& velocityValue_;

    //- Water wave
    autoPtr<waterWaveModels::waterWaves> wave_;

    // Used mass source wave maker or not according to the subDict
    // This is because an object will be stated no matter it is used or not
    bool useMassSource_;

    //- Basic components
        // Phase velocity
        scalar phaseVelocity_;

        // Surface elevation
        scalar waveSurface_;

        // Area of source region
        scalar areaOfSource_;

        // Enhanced factor
        scalar waveEnhanceFactor_;
	
    //- Distribution of mass source
        // Use uneven wavemaker or not (uniform wavemaker by default)
        bool uneven_;

        // Distribution field (for uneven mass source)
        volScalarField distributionField_;

        // Factor of the exponent function (similar to kd)
        scalar gamma_;

        // Adjust U inside the source region or not
        // This contains two parts:
        // 1. Undesired vertical velocity
        // 2. Symmetry of velocity field
        bool adjustU_;

    //- Adjustment on undesired vertical velocity
        // Adjust or not
        bool adjustUndesiredW_;

        // Top of mass source region (z-direction)
        scalar Zt_;

        // Bottom of mass source region (z-direction)
        scalar Zb_;

        // Constant item of distribution function
        scalar Zc_;

        // Start position of vertical velocity adjustment
        // Nondimensional, the value is in range (0,1), Zs = coeff*(Zt-Zb) + Zb
        scalar UDampingCoeff_;

        // Index of the adjustment function
        scalar UDampingIndex_;

        // Vertical velocity damping field (inside mass source)
        volScalarField uAdjustField_;

        // Adjust scheme
        word AdjustmentScheme_;

    //-Adjustment on the symmetry of velocity field
        // adjust or not
        //bool adjustSymmetry_;

        // two dynamic list saves the label of the symmetric cells
        //DynamicList<label> x_left_;
        //DynamicList<label> x_right_;

    //- Functions
        // Calculate phase velocity
        void CalculatePhaseVelocity();

        // Calculate wave surface
        void CalculateWaveSurface(scalar t);

        // Calculate distribution field for uneven mass source
        void CalDistributionField();

        // Calculate vertical velocity damping field
        void CalUAdjustField();

        // start position of vertical velocity adjustment (z-direction)
        scalar Zs();

        // Calculate the twow dynamiclist for U adjustment
        //void CalDynamicLists();

public:

    // Constuctor from the dictionary
    massSourceWaveMaker
    (
        const IOdictionary& dict,
        volScalarField& massSource,
        volScalarField& waveSourceField,
        volVectorField& velocityField
    );

    // Destructor
    ~massSourceWaveMaker(){};

    // Update wave source value at each time step
    void Update(scalar t);

    // Adjust velocity
    void AdjustVelocity();
};	

}//end of namespace
#endif
