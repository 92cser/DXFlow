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
    Foam::wavePorosityMedia

Description
    This class is identical to the one in OLAFlow/IHFoam. See OLAFlow/IHFoam for more
    information.

    Note that currently this class works with DXFlow (due to MULES) ONLY .

    Reference:
    Higuera, P., Lara, J.L. and Losada, I.J. (2014). Three-dimensional interaction of
    waves and porous coastal structures using OpenFOAM. Part I: Formulation and validation.
    Coastal Engineering, 83, 243-258.

    Higuera, P., Lara, J.L. and Losada, I.J. (2014). Three-dimensional interaction of
    waves and porous coastal structures using OpenFOAM. Part II: Application.
    Coastal Engineering, 83, 259–270.

    Wang, D., Dong, S., Sun, J. (2019). Numerical modeling of the interactions between
    waves and a jarlan-type caisson breakwater using openfoam. Ocean Engineering,
    188, 106230.

    Wang, D., Sun, D., Dong, S. (2022). Numerical investigation into effect of the rubble
    mound inside perforated caisson breakwaters under random sea states. Proceedings of
    the Institution of Mechanical Engineers, Part M: Journal of Engineering for the Maritime
    Environment, 236(1), 48-61.

Usage
    See OLAFlow or IHFoam.

SourceFiles
    wavePorosityMedia.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef wavePorosityMedia_H
#define wavePorosityMedia_H

#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
class wavePorosityMedia
{
	// Mark the field of porosity media
	volScalarField& markField_;
	
	// porosity field
	volScalarField& porosity_Field_;
	
	// The number of porosity materials
	const int NPorosityMaterials_;
	
	// Switch of use the porosity or not
	bool usePorosity_;
	
	// Array of coefficients 
        scalarList a_poro_;
        scalarList b_poro_;
        scalarList c_poro_;
        scalarList KC_;
        scalarList D50_;
        scalarList porosity_;
	
	// coefficient field
	volScalarField a_poro_Field_;
	volScalarField b_poro_Field_;
	volScalarField c_poro_Field_;
	
	// KC number field
	volScalarField KC_Field_;
	
	// average dimension of the porosity media field
	volScalarField D50_Field_;
	
	// update the value of fields
	void updateField();
	
	dimensionedScalar balance_;

public:

	// Constuctor from the dictionary
	wavePorosityMedia
		(
			const IOdictionary& dict,
			volScalarField& markField,
			volScalarField& valueField
		);
	
	// Destructor
	~wavePorosityMedia(){};
	
	// Access
	const volScalarField A_item();
	const volScalarField B_item();
	const volScalarField C_item();
	const volScalarField Poro_inverse_();

};
} // end of the namespace

#endif
