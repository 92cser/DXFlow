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

#include "Solitary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
	defineTypeNameAndDebug(Solitary, 0);
    addToRunTimeSelectionTable(waterWaves, Solitary, nameList);
}
}

// Constructor
Foam::waterWaveModels::Solitary::Solitary
(	
	const IOdictionary& dict
)
:
	waterWaves(dict),
	amplitude_(readScalar(WWParameters_.lookup("SolitaryAmplitude"))),
	nonlinearity_(amplitude_/waterDepth_),
	phaseVelocity_(Cal_C_Solitary()),
	k_(Cal_k()),
	IPOWC_(-2*pi_/k_),
	position_(0.0)
{
	Info << "The amplitude of the solitary wave is " << amplitude_ << endl;
}

Foam::word Foam::waterWaveModels::Solitary::name()
{
	return "solitary wave";
}

Foam::scalar Foam::waterWaveModels::Solitary::Cal_k()
{
	scalar para_1 = 0.75*nonlinearity_/waterDepth_/waterDepth_;
	scalar para_2 = 0.0;
	scalar Ki_[9] = {1.0, -0.625000, 0.554688, -0.561535, 0.567095, -0.602969, 0.624914, -0.670850, 0.700371};
	for(int i=0;i<9;i++)
	{
		para_2 += Ki_[i]*pow(nonlinearity_,i);
	}
	return pow(para_1,0.5)* para_2;
}

Foam::scalar Foam::waterWaveModels::Solitary::Cal_C_Solitary()
{
	scalar para_1 = g_*waterDepth_;
	scalar para_2 = 0.0;
	scalar Ci_[10] = {1.0, 1.0, -0.05, -0.042857, -0.034286, -0.031520, -0.029278, -0.026845, -0.030263, -0.021935};
	for(int i=0;i<10;i++)
	{
		para_2 += Ci_[i]*pow(nonlinearity_,i);
	}
	scalar C_2 = para_1* para_2;
	return pow(C_2, 0.5);
}

Foam::scalar Foam::waterWaveModels::Solitary::Cal_C()
{
	return phaseVelocity_;
}

Foam::scalar Foam::waterWaveModels::Solitary::eta(scalar t)
{
	scalar KX = k_*(position_-phaseVelocity_*t - IPOWC_);
	scalar S2 = 1.0/pow(cosh(KX),2);
	scalar para_1 = 0;
	scalar Eta_i_[10];
	Eta_i_[0] = 0.0; //water depth is not included
	Eta_i_[1] = S2;
	Eta_i_[2] = -0.75*S2 + 0.75*pow(S2,2);
	Eta_i_[3] = 0.625*S2 - 1.8875*pow(S2,2) + 1.2625*pow(S2,3);
	Eta_i_[4] = -1.36817*S2 + 3.88033*pow(S2,2) - 4.68304*pow(S2,3) + 2.17088*pow(S2,4);
	Eta_i_[5] = 1.86057*S2 - 7.45136*pow(S2,2) + 12.7637*pow(S2,3) - 11.4199*pow(S2,4) + 4.24687*pow(S2,5);
	Eta_i_[6] = -2.57413*S2 + 13.2856*pow(S2,2) -31.1191*pow(S2,3) + 40.1068*pow(S2,4) - 28.4272*pow(S2,5) + 8.728*pow(S2,6);
	Eta_i_[7] = 3.4572*S2 - 22.782*pow(S2,2) + 68.258*pow(S2,3) - 116.974*pow(S2,4) + 120.49*pow(S2,5) - 71.057*pow(S2,6) + 18.608*pow(S2,7);
	Eta_i_[8] = -4.6849*S2 + 37.67*pow(S2,2) - 139.28*pow(S2,3) + 301.442*pow(S2,4) - 411.416*pow(S2,5) + 355.069*pow(S2,6) - 180.212*pow(S2,7) + 41.412*pow(S2,8) ;
	Eta_i_[9] = 6.191*S2 - 60.57*pow(S2,2) + 269.84*pow(S2,3) -712.125*pow(S2,4) + 1217.98*pow(S2,5) - 1384.37*pow(S2,6) + 1023.07*pow(S2,7) - 450.29*pow(S2,8) + 90.279*pow(S2,9);
	for(int i=0;i<10;i++)
	{
		para_1 += Eta_i_[i]*pow(nonlinearity_,i);
	}
	scalar eta_ = waterDepth_*para_1;
	Info << "Current time and theoretical eta = " << t << ", " << eta_ << endl;
   	return eta_;
}


Foam::scalar Foam::waterWaveModels::Solitary::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
	scalar PistonVelocity_ = phaseVelocity_*eta(t) / (waterDepth_+eta(t));
	position_ += PistonVelocity_*deltaT;
	return position_;
}

void Foam::waterWaveModels::Solitary::PrintWaveProperties()
{
    Info << "---------------------------------" << "\n"
         << "Wave properties of Solitary Wave" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Ampitude (m) = " << amplitude_ << "\n"
         << "---------------------------------" << endl;
}
