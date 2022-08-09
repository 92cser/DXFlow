/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 Dalian Ocean University / National Marine Environment Monitoring Center / Dalian University of Technology
    Copyright (C) 2018-2022 Ocean University of China
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
    Foam::waveDamping

Description:
    Source terms for wave damping.

    Reference:
    Wang, D., Dong, S. (2022). A discussion of numerical wave absorption
    using sponge layer methods. Ocean Engineering, 247, 110732.

Usage:
    in inputDict
    WaveDamping
    {
            // optional: adding a wave damping source to the VOF equation or not
            AlphaSource				false;

            // optional: using uneven forcing strength or not
            UnevenForcingStrength
            {
                    // wave number
                    waveNumber                          0.1;

                    // z-position of the top of the uneven sponge layer
                    // usually set as A = H/2 above the SWL
                    Spt                                 0.02;

                    // z-position of the bottom of the uneven sponge layer
                    // usually set as the bottom of the computational domain
                    Spb                                 -0.4;

                    // z-position of still water level
                    zSWL                                0.0;

            }

            // optional: adding source terms at four different positions
            x+
            {
                    ForcingStrength     		10.0; //forcing strength
                    Index				2.0;   // optional: index
                    StartPosition	    		2.0;
                    EndPosition                         4.0;
            }
            x-
            {	// consistent with the direction of incident waves
                    // therefore the start point should be larger than the end point
                    ForcingStrength                     10.0;
                    Index				2.0;
                    StartPosition                       0.0;
                    endPosition                         -2.0;
            }
            y+
            {
                    ForcingStrength     		10.0;
                    Index				2.0;
                    StartPosition	    		2.0;
                    EndPosition                         4.0;
            }
            y-
            {	// consistent with the direction of incident waves
                    // therefore the start point should be larger than the end point
                    ForcingStrength                     10.0;
                    Index				2.0;
                    StartPosition                       0.0;
                    endPosition                         -2.0;
            }

    }

SourceFiles
    waveDamping.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef waveDamping_H
#define waveDamping_H

#include "volFields.H"

namespace Foam
{
class waveDamping
{
	class DampingLayer
	{
            // mesh
            const fvMesh& mesh_;

            const word location_;

            // use wave damping region or not
            const bool Valid_;

            // x-direction factors
            // start position
            scalar xs_;
            // start position
            scalar xe_;
            // x-directional width of the sponge layer
            scalar xWidth_;
            // blending function, calculated as per x coordinate.
            // Exponent form is considered, ranging from 0 to 1
            volScalarField blendingField_;
            // the index used for horizontal blending function
            scalar index_;
            // uniform forcing Strength, calculated using the analytical model of Peric (2018), Ocean engineering
            scalar gamma0_;


            // z-direction factors
            // wave number
            scalar k_;
            // top of the uneven sponge layer. The value is better to be greater than wave amplitude.
            scalar Spt_;
            // bottom of the uneven sponge layer. The value is set as the bottom of the wave flume.
            scalar Spb_;
            // z-directional width of the sponge layer
            scalar zWidth_;
            // the z-directional position of still water level
            scalar zSWL_;
            // uniform part of uneven forcing Strength, the final uneven sponge layer equals gamma1*distributionFunction
            scalar gamma1_;
            // water depth
            scalar d_;


            // final wave damping source
            volScalarField ForceAndBlendingField_;

	public:
		// constructor
		DampingLayer
		(
			const fvMesh& mesh,
			const dictionary& inputSubDict,
			const word& location,
                        bool  uneven
		);
		
		// access
		const volScalarField blendingField();

		// access
		const volScalarField ForceAndBlendingField();
	};


	// mesh
	const fvMesh& mesh_;	

        // adding source term to alpha equation or not
        bool alphaOrNot_;
	
        // adopt (vertical) uneven forcing strength or not
        bool unevenOrNot_;

        // four damping layers used at each end of the wave tank
	DampingLayer X_pos_;
	
	DampingLayer X_neg_;
	
	DampingLayer Y_pos_;
	
	DampingLayer Y_neg_;
	
	// total Blending field
	const volScalarField BlendingField_;

	// total daming field
	const volScalarField DampingField_;
	
	const dimensionedScalar DimensionBalance_;
	
public:
	// constructor
	waveDamping
	(
		const fvMesh& mesh,
		const dictionary& inputDict
	);
	
        //volScalarField alphaCurrRef(volScalarField& alphaRef0);

	// access
	const volScalarField DampingField();

	const volScalarField BlendingField();

        scalar AlphaSource();

}; // end of class  
} // end of namespace Foam
#endif
