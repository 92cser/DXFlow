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

\*---------------------------------------------------------------------------*/

#include "wavePorosityMedia.H"
#include "fvCFD.H"
#include "argList.H"

// constructor
Foam::wavePorosityMedia::wavePorosityMedia
(
	const IOdictionary& dict,
	volScalarField& markField,
	volScalarField& valueField
)
:
	markField_(markField),
	porosity_Field_(valueField),
	NPorosityMaterials_(gMax(markField_)),
	usePorosity_(dict.isDict("PorosityProperties")),
// Arrays:
	a_poro_(NPorosityMaterials_,0.0),
	b_poro_(NPorosityMaterials_,0.0),
	c_poro_(NPorosityMaterials_,0.0),
	KC_(NPorosityMaterials_,0.0),
	D50_(NPorosityMaterials_,0.0),
	porosity_(NPorosityMaterials_,0.0),
// Fields
	a_poro_Field_(valueField),
	b_poro_Field_(valueField),
	c_poro_Field_(valueField),
	KC_Field_(valueField),
	D50_Field_(valueField),
// dimension balance
	balance_("balance_",dimensionSet(0,1,0,0,0,0,0),1.0)
{
	a_poro_Field_*= 0.0;
	b_poro_Field_*= 0.0;
	c_poro_Field_*= 0.0;
	D50_Field_ *= balance_;
	if(usePorosity_)
	{
                a_poro_ = dict.subDict("PorosityProperties").lookup("a_coeff");
                b_poro_ = dict.subDict("PorosityProperties").lookup("b_coeff");
                c_poro_ = dict.subDict("PorosityProperties").lookup("c_coeff");
                KC_ = dict.subDict("PorosityProperties").lookup("KCNumber");
                D50_ = dict.subDict("PorosityProperties").lookup("meanDiameter");
                porosity_ = dict.subDict("PorosityProperties").lookup("porosity");
		updateField();
		Info << "The porosity model is used." << nl 
			 << "There are totally " << NPorosityMaterials_ << " porosity materials." << nl
			 << "The length of arrays in inputDict should also be " << NPorosityMaterials_ << "." << endl;
		Info << "Mean diameters are: " << D50_ << nl
                        << "Porosities are:" << porosity_ << endl;
	}
}

void Foam::wavePorosityMedia::updateField()
{
	forAll(markField_,i)
	{
		if(markField_[i] > 0)
		{
			a_poro_Field_[i] = a_poro_[markField_[i]-1];
			b_poro_Field_[i] = b_poro_[markField_[i]-1];
			c_poro_Field_[i] = c_poro_[markField_[i]-1];
			KC_Field_[i] = KC_[markField_[i]-1];
			D50_Field_[i] = D50_[markField_[i]-1];
			porosity_Field_[i] = porosity_[markField_[i]-1];
		}
	}
}

const Foam::volScalarField Foam::wavePorosityMedia::A_item()
{
	return a_poro_Field_*pow(1.0-porosity_Field_, 2)/pow3(porosity_Field_)/pow(D50_Field_,2);
}

const Foam::volScalarField Foam::wavePorosityMedia::B_item()
{
	return b_poro_Field_*(1.0-porosity_Field_) / pow3(porosity_Field_) / D50_Field_ * (1.0 + 7.5 / KC_Field_);
}

const Foam::volScalarField Foam::wavePorosityMedia::C_item()
{
	return (1+c_poro_Field_)/porosity_Field_;
}

const Foam::volScalarField Foam::wavePorosityMedia::Poro_inverse_()
{
	return 1/porosity_Field_;
}
