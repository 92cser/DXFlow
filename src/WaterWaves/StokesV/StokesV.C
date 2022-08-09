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

#include "StokesV.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(StokesV, 0);
    addToRunTimeSelectionTable(waterWaves, StokesV, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

// Constructor
Foam::waterWaveModels::StokesV::StokesV
(	
	const IOdictionary& dict
)
:
    monochromaticWave(dict),
    waveLength_V_(waveLength_), // initialised using the linear wave length
    B_22(0.0),
    B_24(0.0),
    B_33(0.0),
    B_35(0.0),
    B_44(0.0),
    B_55(0.0),
    C1(0.0),
    C2(0.0)
{
    Solve();
}

Foam::scalar Foam::waterWaveModels::StokesV::K()
{
    return 2*pi_/waveLength_V_;
}

Foam::scalar Foam::waterWaveModels::StokesV::CH()
{
    return cosh(K()*waterDepth_);
}

Foam::scalar Foam::waterWaveModels::StokesV::SH()
{
    return sinh(K()*waterDepth_);
}

Foam::scalar Foam::waterWaveModels::StokesV::K1()
{
    return waveHeight_/waterDepth_;
}

Foam::scalar Foam::waterWaveModels::StokesV::K2()
{
    return 2*pi_*waterDepth_/g_/wavePeriod_/wavePeriod_;
}

Foam::scalar Foam::waterWaveModels::StokesV::Q()  
{
    return K2() / d_L() / tanh(K()*waterDepth_);
}

Foam::scalar Foam::waterWaveModels::StokesV::Lambda()
{
    scalar para_1 = C1*C1-4*C2*(1-Q());
    scalar para_2 = (-C1 + pow(para_1,0.5))/2.0/C2;
    if (para_2 < 0)
    {
         para_2 = 0;
    }
    return pow(para_2,0.5);
}

Foam::scalar Foam::waterWaveModels::StokesV::d_L()
{
    return waterDepth_ / waveLength_V_;
}

void Foam::waterWaveModels::StokesV::Cal_B_22()
{
    scalar ch = CH();
    scalar sh = SH();
    B_22 = 0.25*(2*ch*ch+1)*ch/pow(sh,3);
}

void Foam::waterWaveModels::StokesV::Cal_B_24()
{
    scalar ch = CH();
    scalar sh = SH();
    B_24 = (272*pow(ch,8) - 504*pow(ch,6) - 192*pow(ch,4) + 322*pow(ch,2) + 21)*ch/ 384.0 / pow(sh,9);
}

void Foam::waterWaveModels::StokesV::Cal_B_33()
{
    scalar ch = CH();
    scalar sh = SH();
    B_33 = 3*(8*pow(ch,6)+1)/64.0/pow(sh,6);
}

void Foam::waterWaveModels::StokesV::Cal_B_35()
{
    scalar ch = CH();
    scalar sh = SH();
    scalar para_1 = 88128*pow(ch,14)-208224*pow(ch,12)+70848*pow(ch,10)+54000*pow(ch,8)-21816*pow(ch,6)+6264*pow(ch,4)-54*pow(ch,2)-81;
    scalar para_2 = 12288*pow(sh,12)*(6*pow(ch,2)-1);
    B_35 = para_1 / para_2;
}

void Foam::waterWaveModels::StokesV::Cal_B_44()
{
    scalar ch = CH();
    scalar sh = SH();
    scalar para_1 = (768*pow(ch,10)-488*pow(ch,8)-48*pow(ch,6)+48*pow(ch,4)+106*pow(ch,2)-21)*ch;
    scalar para_2 = 384*pow(sh,9)*(6*pow(ch,2)-1);
    B_44 = para_1 / para_2;
}

void Foam::waterWaveModels::StokesV::Cal_B_55()
{
    scalar ch = CH();
    scalar sh = SH();
    scalar para_1 = 192000*pow(ch,16)-262720*pow(ch,14)+83680*pow(ch,12)+20160*pow(ch,10)-7280*pow(ch,8)+7160*pow(ch,6)-1800*pow(ch,4)-1050*pow(ch,2)+225;
    scalar para_2 = 12288*pow(sh,10)*(6*pow(ch,2)-1)*(8*pow(ch,4)-11*pow(ch,2)+3);
    B_55 = para_1 / para_2;
}

void Foam::waterWaveModels::StokesV::Cal_C1()
{
    scalar ch = CH();
    scalar sh = SH();
    C1 = (8*pow(ch,4)-8*pow(ch,2)+9)/8.0/pow(sh,4);
}

void Foam::waterWaveModels::StokesV::Cal_C2()
{
    scalar ch = CH();
    scalar sh = SH();
    scalar para_1 = 3840*pow(ch,12)-4096*pow(ch,10)+2592*pow(ch,8)-1008*pow(ch,6)+5944*pow(ch,4)-1830*pow(ch,2)+147;
    scalar para_2 = 512*pow(sh,10)*(6*pow(ch,2)-1);
    C2 = para_1 / para_2;
}

void Foam::waterWaveModels::StokesV::Cal_All()
{
    Cal_B_22();
    Cal_B_24();
    Cal_B_33();
    Cal_B_35();
    Cal_B_44();
    Cal_B_55();
    Cal_C1();
    Cal_C2();
}

void Foam::waterWaveModels::StokesV::Solve()
{
    Cal_All();
    scalar Lambda_ = Lambda();
    scalar fx_old;
    scalar fx = K1()*pi_*d_L()-( Lambda_ + pow(Lambda_,3) * B_33 + pow(Lambda_,5)*(B_35+B_55) );
    do
    {
        waveLength_V_ += Vsmall_;
        Cal_All();
        Lambda_ = Lambda();
        fx_old = fx;
        fx = K1()*pi_*d_L()-( Lambda_ + pow(Lambda_,3) * B_33 + pow(Lambda_,5)*(B_35+B_55) );
        /*cout << "    *_old    ="<< fx_old << endl;
        cout << "    *    ="<< fx << endl;
        cout << "    Lambda_    ="<< Lambda_ << endl;*/
    }while (fx*fx_old > 0.0);
    //cout << "the left value " << waveLength_V_ - Vsmall_ <<endl;
    //cout << "the right value " << waveLength_V_ << endl;
    waveLength_V_ -= 0.5*Vsmall_;
    Cal_All();
}


Foam::word Foam::waterWaveModels::StokesV::name()
{
	return "StokesV wave";
}

Foam::scalar Foam::waterWaveModels::StokesV::eta(scalar t)
{
    scalar Lambda_ = Lambda();
    scalar eta_ = 1/K()*(
                        Lambda_*cos(Omega_*t)
                        + (pow(Lambda_,2)*B_22 + pow(Lambda_,4)*B_24)*cos(2*Omega_*t)
                        + (pow(Lambda_,3)*B_33 + pow(Lambda_,5)*B_35)*cos(3*Omega_*t)
                        + pow(Lambda_,4)*B_44*cos(4*Omega_*t)
                        + pow(Lambda_,5)*B_55*cos(5*Omega_*t)
                       );
	// slowly start factor
    scalar SSF_ = (wavePeriod_ > 0.5*t) ? pow(0.5*t/wavePeriod_,2.0) : 1.0;
    return eta_*SSF_;
}

Foam::scalar Foam::waterWaveModels::StokesV::Cal_C()
{
    scalar Lambda_ = Lambda();
    scalar para_1 = g_*tanh(K()*waterDepth_)*(1 + pow(Lambda_,2)*C1 + pow(Lambda_,4)*C2);
    scalar para_2 = para_1 /K();
    return pow(para_2,0.5);
}

Foam::scalar Foam::waterWaveModels::StokesV::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
   return monochromaticWave::displacement(t, deltaT, WaterLevel);
}

void Foam::waterWaveModels::StokesV::PrintWaveProperties()
{
    Info << "-------------------------------" << "\n"
         << "Wave properties of StokesV wave" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Wave Height (m) = " << waveHeight_ << "\n"
         << " Wave Period (s) = " << wavePeriod_ << "\n"
         << " Wave Length (m) = " << waveLength_V_ << "\n"
         << "-------------------------------" << endl;
}
