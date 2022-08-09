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

#include "StreamFunction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(StreamFunctionWave, 0);
    addToRunTimeSelectionTable(waterWaves, StreamFunctionWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

Foam::waterWaveModels::StreamFunctionWave::StreamFunctionWave
(
    const IOdictionary& dict
)
:
    waterWaves(dict),
    //Bj_(WWParameters_.lookup("Bj")),
    Ej_(WWParameters_.lookup("Ej")),
    Properties_((WWParameters_.lookup("Properties"))),
    position_(0.0)
{
    updateWaveProperties();
}

Foam::word Foam::waterWaveModels::StreamFunctionWave::name()
{
    return "Stream function wave";
}

Foam::scalar Foam::waterWaveModels::StreamFunctionWave::Cal_L
(
    scalar D,
    scalar T
)
{
    scalar pi_ = 3.1415926;
    scalar g_ = 9.81;
    scalar mu0_ = pow(2*pi_,2)*D/g_*pow(T,-2);
    scalar index0_ = 1.835+1.225*pow(mu0_,1.35);
    scalar mu_ = mu0_*(1+mu0_*exp(-index0_))*pow(tanh(mu0_),-0.5);
    // another iteration to improve precise
    scalar muFinal_ = (pow(mu_,2)+mu0_*pow(cosh(mu_),2))/(mu_+0.5*sinh(2*mu_));
    return 2*pi_*D/muFinal_;
}

void Foam::waterWaveModels::StreamFunctionWave::updateWaveProperties()
{
    scalar pi_ = 3.1415926;
    scalar g_ = 9.81;
    scalar gd_ = g_*waterDepth_;
    L_ = Properties_[0]*waterDepth_;
    k_ = 2*pi_/L_;
    H_ = Properties_[1]*waterDepth_;
    T_ = Properties_[2]/sqrt(g_/waterDepth_);
    C_ = Properties_[3]*sqrt(gd_);
    /*EC_ = Properties_[4]*sqrt(gd_);
    SC_ = Properties_[5]*sqrt(gd_);
    MS_ = Properties_[6]*sqrt(gd_);
    Q_ = Properties_[7]*sqrt(g_*pow(waterDepth_,3));
    R_ = Properties_[8]*gd_;*/
    OmegaC_ = pi_/15;
    ts_ = k_*waterDepth_*(2*k_*waterDepth_ + sinh(2*k_*waterDepth_))/4.0/sinh(k_*waterDepth_)/sinh(k_*waterDepth_);
    if(2*pi_/T_ > 7.0/sqrt(waterDepth_/g_))
    {
        scalar T0 = 2*pi_/OmegaC_;
        scalar L0 = Cal_L(waterDepth_,T0);
        scalar k0 = 2*pi_/L0;
        ts_ = k0*waterDepth_*(2*k0*waterDepth_ + sinh(2*k0*waterDepth_))/4.0/sinh(k0*waterDepth_)/sinh(k0*waterDepth_);
    }
}

Foam::scalar Foam::waterWaveModels::StreamFunctionWave::eta(scalar t)
{
    scalar eta_ = 0.0;
    forAll(Ej_, i)
    {
        eta_ += Ej_[i]*cos((i+1)*k_*(-C_*t));
    }
    // slowly start factor
    scalar SSF_ = (t < 2*T_) ? pow(0.5*t/T_,2.0) : 1.0;
    return SSF_*eta_/k_;
}

Foam::scalar Foam::waterWaveModels::StreamFunctionWave::etaForPiston(scalar x, scalar t)
{
    scalar eta_ = 0.0;
    forAll(Ej_, i)
    {
        eta_ += Ej_[i]*cos((i+1)*k_*(x-C_*t));
    }
    return eta_/k_;
}

Foam::scalar Foam::waterWaveModels::StreamFunctionWave::Cal_C()
{
    return C_;
}

Foam::scalar Foam::waterWaveModels::StreamFunctionWave::displacement(scalar t, scalar deltaT, scalar correctWaterLevel)
{
    // store the position of last time step
    scalar position0_ = position_;
    // a guessed position
    scalar x_guess = position0_;
    // corresponding eta at the wave paddle
    scalar etaForPiston_ = etaForPiston(x_guess, t);
    // the calculated paddle position
    scalar x_new = (C_*etaForPiston_/(waterDepth_+ etaForPiston_) - OmegaC_*x_guess)*deltaT + position0_;
    // iteration
    while (mag((x_new-x_guess)/x_new) > 0.001)
    {
        x_guess = x_new;
        etaForPiston_ = etaForPiston(x_guess, t);
        x_new = (C_*etaForPiston_/(waterDepth_+ etaForPiston_) - OmegaC_*x_guess)*deltaT + position0_;
    }
    position_ = x_new;
    return position_*ts_;
}

void Foam::waterWaveModels::StreamFunctionWave::PrintWaveProperties()
{
    Info << "-----------------------------------------------" << "\n"
         << "Wave properties of stream function wave" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Wave Height (m) = " << H_ << "\n"
         << " Wave Period (s) = " << T_ << "\n"
         << " Wave Length (m) = " << L_ << "\n"
         << " Wave Speed (m/s) = " << C_ << "\n"
         << "-----------------------------------------------" << endl;
}
