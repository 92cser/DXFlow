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
    Foam::PistonWaveGeneration

Description
    Piston wavemaker class
    
    reference:
    Wang, D., Sun, J., Gui, J. et al. (2019). A numerical piston-type wave-maker 
    toolbox for the open-source library OpenFOAM. Journal of Hydrodynamics 31, 800–813. 

Usage
    in pointDisplacement
    inlet
    {
        type            PistonWaveGeneration;
        value           uniform (0 0 0);
    }

SourceFiles
    PistonWaveGeneration.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef PistonWaveGeneration_H
#define PistonWaveGeneration_H

#include "fixedValuePointPatchField.H"
#include "IOdictionary.H"
#include "WaterWaves.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class PistonWaveGeneration Declaration
\*---------------------------------------------------------------------------*/

class PistonWaveGeneration
:
    public fixedValuePointPatchField<vector>
{
    // Private data

	// input dictionary
	const IOdictionary dict_;

	// final displacement of the patch
	vector Displacement_;

	// Water wave
	autoPtr<waterWaveModels::waterWaves> wave;
	
	// water level
	scalar waterLevel_;

public:

    //- Runtime type information
    TypeName("PistonWaveGeneration");


    // Constructors

        //- Construct from patch and internal field
        PistonWaveGeneration
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        PistonWaveGeneration
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        PistonWaveGeneration
        (
            const PistonWaveGeneration&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new PistonWaveGeneration
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        PistonWaveGeneration
        (
            const PistonWaveGeneration&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new PistonWaveGeneration
                (
                    *this,
                    iF
                )
            );
        }
	
	// Destructor
	virtual ~PistonWaveGeneration(){};

	//Constant 
	const scalar pi_ = 3.141592653589793;
    const scalar g_ = 9.81;

    // Member functions

        // Evaluation functions
		// Calculate the water level at this patch
		void Cal_waterLevel
		(
			const scalarField& boundarySurface,
			const scalarField& boundaryAlpha,
			scalar zSpan
		);

		//- Update the coefficients associated with the patch field
		virtual void updateCoeffs();

		//- Write
		virtual void write(Ostream&) const;

	//#include "PistonWaveGenerationI.H"
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
