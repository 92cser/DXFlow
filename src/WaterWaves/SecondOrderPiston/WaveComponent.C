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
#include "WaveComponent.H"
namespace Foam
{
waveComponent::waveComponent
(
    scalar waveHeight,
    scalar waterDepth,
    scalar wavePeriod,
    bool light
)
:
    H_(waveHeight),
    D_(waterDepth),
    T_(wavePeriod),
    Omega_(2*pi_/wavePeriod),
    L_(CalWaveLength(D_,T_)),
    k_(51, complex(0,0)),
    TF_(51, complex(0,0))
{
    if(light)
    {
        k_[0] = {2*pi_/L_,0.0};
    }
    else
    {
        Cal_k();
        Cal_TF();
    }
}

// calculate wave length
scalar waveComponent::CalWaveLength
(
    scalar D,
    scalar T
)
{
    scalar mu0_ = pow(2*pi_,2)*D/g_*pow(T,-2);
    scalar index0_ = 1.835+1.225*pow(mu0_,1.35);
    scalar mu_ = mu0_*(1+mu0_*exp(-index0_))*pow(tanh(mu0_),-0.5);
    // another iteration to improve precise
    scalar muFinal_ = (pow(mu_,2)+mu0_*pow(cosh(mu_),2))/(mu_+0.5*sinh(2*mu_));
    return 2*pi_*D/muFinal_;
}

void waveComponent::Cal_k()
{
    for(int i=0; i<51; i++)
    {
        if(i==0)//store real wave number
        {
            k_[i] = {2*pi_/L_,0.0};
        }
        else //solve and store imaginary wave number
        {
            k_[i] = {0.0, k_im_bisection((i-0.5)*pi_,(i+0.5)*pi_,pi_)};
        }
    }
}

scalar waveComponent::k_im_bisection // take A = kd as a variable
(
    scalar A_left,
    scalar A_right,
    scalar deltaA
)
{
    // prevent divergence
    A_left += deltaA/50;
    A_right -= deltaA/50;

    while(mag(A_left-A_right) > 1e-6)
    {
        scalar A_middle = 0.5*(A_left + A_right);
        scalar fx_left = pow(Omega_,2)*D_+g_*A_left*tan(A_left);
        scalar fx_right = pow(Omega_,2)*D_+g_*A_right*tan(A_right);
        scalar fx_mid = pow(Omega_,2)*D_+g_*A_middle*tan(A_middle);
        if(fx_left*fx_mid < 0)
        {
            A_right = A_middle;
        }
        else
        {
            A_left = A_middle;
        }
    }
    return 0.5*(A_left+A_right)/D_;
}

void waveComponent::Cal_TF()
{
    for(int i=0; i<51; i++)
    {
        complex kd = k_[i]*D_;
        TF_[i] = 4.0*sinhC(kd)*sinhC(kd)/(2.0*kd + sinhC(2.0*kd));
    }
}

complex waveComponent::A(scalar t)
{
    scalar phi = Omega_*t;
    scalar A0 = 0.5*H_;
    complex Amp(A0*cos(phi),A0*sin(phi));
    return Amp;
}
}// end of namespace
