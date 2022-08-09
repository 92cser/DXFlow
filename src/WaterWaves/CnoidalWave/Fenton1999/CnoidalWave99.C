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

#include "CnoidalWave99.H"
#include "addToRunTimeSelectionTable.H"
#include "DynamicList.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(CnoidalWave, 0);
    addToRunTimeSelectionTable(waterWaves, CnoidalWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

Foam::waterWaveModels::CnoidalWave::CnoidalWave
(	
    const IOdictionary& dict
)
:
    waterWaves(dict),
    waveHeight_(readScalar(WWParameters_.lookup("WaveHeight"))),
    waveLength_(0.0),
    wavePeriod_(readScalar(WWParameters_.lookup("WavePeriod"))),
    phaseVelocity_(0.0),
    m_(0.5),
    h_(0.0),
    Epsl_(waveHeight_/waterDepth_),
    position_(0.0)
{
    solve();
    updateCoeffs();
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::K()
{
    scalar para_1 = 1 + pow(m_, 0.25);
    scalar para_2 = 2.0*(1 + pow(m_, 0.25))/(1-pow(m_, 0.25));
    return 2.0/pow(para_1,2.0)*log(para_2);
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::K0()
{
    scalar para_1 = 1 + pow(m_, 0.25);
    return 2.0*pi_/pow(para_1,2.0);
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::e()
{
    scalar para_1 = pow(q1(),2.0);
    scalar K_ = K();
    scalar K0_ = K0();
    return (2.0-m_)/3 + 0.5*pi_/K_/K0_ + 2.0*pow(pi_/K0_,2)*(-1./24 + para_1/pow(1 - para_1,2.0));
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::q1()
{
    return exp(-pi_*K()/K0());
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::E()
{
    return K()*e();
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::hbyD()
{
    scalar e_ = e();
    scalar para_1 = 573./2000*e_ - 57./400*pow(e_,2.0) + 0.25*pow(e_,3.0);
    scalar para_2 = -302159./1470000*e_ + 1779./2000*pow(e_,2.0) - 123./400*pow(e_,3.0) + 0.25*pow(e_,4.0);
    return 1 + Epsl_*(-e_) + pow(Epsl_,2.0)*0.25*e_ + pow(Epsl_,3.0)*(-e_/25 + 0.25*pow(e_,2.0)) + pow(Epsl_,4.0)*para_1 + pow(Epsl_,5.0)*para_2;
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::Hbyh()
{
    return Epsl_/hbyD();
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::Alpha()
{
    scalar Hbyh_ = Hbyh();
    scalar para_1 = 0.75*Hbyh_;
    return pow(para_1,0.5)* (1 - 5./8*Hbyh_ + 71./128*pow(Hbyh_,2.0) - 100627./179200*pow(Hbyh_,3.0) + 16259737./28672000*pow(Hbyh_,4.0));
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::NDPhaseVelocity()
{
    scalar Hbyh_ = Hbyh();
    scalar e_ = e();
    return 1 + Hbyh_*(0.5-e_) + pow(Hbyh_,2.0)*(-3./20 + 5./12*e_) + pow(Hbyh_,3.0)*(3./56 - 19./600*e_) + pow(Hbyh_,4.0)*(-309./5600 + 3719./21000*e_) + pow(Hbyh_,5.0)*(12237./616000- 997699./8820000*e_);
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::tolerance(scalar m0)
{
    m_ = m0;
    scalar para_1 = wavePeriod_*pow(g_/waterDepth_,0.5);
    scalar hbyD_ = hbyD();
    return NDPhaseVelocity()*pow(hbyD_,0.5) - 1/para_1*hbyD_*2*K()/Alpha();
}

void Foam::waterWaveModels::CnoidalWave::solve()
{
    scalar m_left = 1e-10;
    scalar m_right = 1-m_left;
    while(mag(m_left-m_right)> 1e-6)
    {
        scalar tol_left = tolerance(m_left);
        scalar tol_right = tolerance(m_right);
        scalar m_middle = 0.5*(m_left + m_right);
        scalar tol_middle = tolerance(m_middle);
        if(tol_left*tol_middle > 0 && tol_right*tol_middle < 0)
        {
            m_left = m_middle;
        }
        else
        {
            m_right = m_middle;
        }
    }
    m_ = 0.5*(m_left + m_right);
    //cout << " m = " << m_ << endl;
}

void Foam::waterWaveModels::CnoidalWave::updateCoeffs()
{
    h_ = hbyD()*waterDepth_;
    phaseVelocity_ = NDPhaseVelocity()*pow(g_*h_,0.5);
    waveLength_ = wavelength()*waterDepth_;
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::cn(scalar z)
{
    scalar q1_ = q1();
    scalar w_ = pi_*z/2/K0();
    scalar para_1 = pow((1-m_)/m_/q1_, 0.25);
    if(1-m_ < 1e-8)
    {
        para_1 = 1;
    }
    return 0.5*para_1*(1 - 2.0*q1_*cosh(2.0*w_))/(cosh(w_) +pow(q1_,2.0)*cosh(3.0*w_));
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::wavelength()
{
    scalar para_1 = 4.0*K()/pow(3*Epsl_, 0.5);
    scalar e_ = e();
    scalar para_2 = -21./128+1./16*e_+3./8*pow(e_,2);
    scalar para_3 = 20127./179200 - 409./6400*e_ + 7./64*pow(e_,2) + 1./16*pow(e_,3);
    scalar para_4 = -1575087./28672000 + 1086367./1792000*e_ - 2679./25600*pow(e_,2) + 13./128*pow(e_,3) + 3./128*pow(e_,4);
    return para_1*(1 + Epsl_*(5./8-1.5*e_) + pow(Epsl_,2)*para_2 + pow(Epsl_,3)*para_3 + pow(Epsl_,4)*para_4);
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::eta(scalar t)
{
    // time shift
    int n = t/wavePeriod_;
    scalar t_ = -0.5*wavePeriod_ + t - n*wavePeriod_;

    scalar Hbyh_ = waveHeight_/h_;
    scalar phase_ = Alpha()*(position_-phaseVelocity_*t_)/h_;  //position = 0
    scalar cn_ = cn(phase_);
    scalar para_1 = -0.75*pow(cn_,2) + 0.75*pow(cn_,4);
    scalar para_2 = 5./8*pow(cn_,2) - 151./80*pow(cn_,4) + 101./80*pow(cn_,6);
    scalar para_3 = -8209./6000*pow(cn_,2) + 11641./3000*pow(cn_,4) - 112393./24000*pow(cn_,6) + 17367./8000*pow(cn_,8);
    scalar para_4 = 364671./196000*pow(cn_,2) - 2920931./392000*pow(cn_,4) + 2001361./156800*pow(cn_,6) - 17906339./1568000*pow(cn_,8) + 1331817./313600*pow(cn_,10);
    // slowly start factor
    scalar SSF_ = (wavePeriod_ > 0.5*t) ? pow(0.5*t/wavePeriod_,2.0) : 1.0;
    return SSF_*((1 + Hbyh_*pow(cn_,2)+ pow(Hbyh_,2)*para_1 + pow(Hbyh_,3)*para_2 + pow(Hbyh_,4)*para_3 + pow(Hbyh_,5)*para_4)*h_ - waterDepth_);
}


void Foam::waterWaveModels::CnoidalWave::PrintWaveProperties()
{
    Info << "-----------------------------------------------" << "\n"
         << "Wave properties of the Fifth-order Cnoidal Wave (Fenton 1999)" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Wave Height (m) = " << waveHeight_ << "\n"
         << " Wave Period (s) = " << wavePeriod_ << "\n"
         << " Wave Length (m) = " << waveLength_ << "\n"
         << " Minimum Water Depth (m) = " << h_ << "\n"
         << " Phase Velocity (m/s) = " << phaseVelocity_ << "\n"
         << " Actual Module Number = " << m_ << "\n"
         << " K = " << K() << "\n"
         << " Straining Factor = " << Alpha() << "\n"
         << "-----------------------------------------------" << endl;
}

Foam::word Foam::waterWaveModels::CnoidalWave::name()
{
    return "Cnoidal wave";
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::Cal_C()
{
    return phaseVelocity_;
}

Foam::scalar Foam::waterWaveModels::CnoidalWave::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
    // slowly start factor
    scalar SSF_ = (wavePeriod_ > 0.5*t) ? pow(0.5*t/wavePeriod_,2.0) : 1.0;
    scalar etaForPiston_ = eta(t)/SSF_;
    scalar PistonVelocity_ = phaseVelocity_*etaForPiston_ / (waterDepth_+etaForPiston_) -0.1*position_;
    position_ += PistonVelocity_*deltaT;
    Info << "TAndX = " << t << ", " << position_ << endl;
    return position_;
}
