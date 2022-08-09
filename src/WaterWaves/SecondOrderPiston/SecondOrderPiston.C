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
#include "SecondOrderPiston.H"
namespace Foam
{
// constructors
SecondOrderWave::SecondOrderWave
(
    waveComponent& wave1,
    waveComponent& wave2
)
:
    wave1_(&wave1),
    wave2_(&wave2),
    Omega_super_(wave1_->Omega() + wave2_->Omega()),
    Omega_sub_(wave1_->Omega() - wave2_->Omega()),
    Regular_(!Omega_sub_),
    delta_(Regular_ ? 0.5 : 1.0),
    K0_super_(K0_super()),
    K0_sub_(K0_sub()),
    E_super_(E_super()),
    E_sub_(E_sub()),
    F11h_super_(F11h_super()),
    F11h_sub_(F11h_sub()),
    F12h_super_(F12h_super()),
    F12h_sub_(F12h_sub()),
    F13h_super_(F13h_super()),
    F13h_sub_(F13h_sub()),
    iF22h_super_(iF22h_super()),
    iF22h_sub_(iF22h_sub()),
    iF23h_super_(iF23h_super()),
    iF23h_sub_(iF23h_sub()),
    iF24h_super_(iF24h_super()),
    iF24h_sub_(iF24h_sub()),
    TF_final_super_(F11h_super_+F12h_super_+F13h_super_-iF22h_super_-iF23h_super_-iF24h_super_),
    TF_final_sub_(F11h_sub_+F12h_sub_+F13h_sub_-iF22h_sub_-iF23h_sub_-iF24h_sub_)
{}

SecondOrderWave::SecondOrderWave
(
    waveComponent& wave1,
    waveComponent& wave2,
    word light
)
:
    wave1_(&wave1),
    wave2_(&wave2),
    Omega_super_(wave1_->Omega() + wave2_->Omega()),
    Omega_sub_(wave1_->Omega() - wave2_->Omega()),
    Regular_(!Omega_sub_)
{}

complex SecondOrderWave::K0_super()
{
    scalar T_super = 2*pi_/Omega_super_;
    scalar L_super = wave1_-> CalWaveLength(wave1_->D(),T_super);
    return complex(2*pi_/L_super, 0);
}

complex SecondOrderWave::K0_sub()
{
    if(Regular_)
    {
        return complex(0,0);
    }
    else
    {
        scalar T_sub = 2*pi_/Omega_sub_;
        scalar L_sub = wave1_-> CalWaveLength(wave1_->D(),T_sub);
        return complex(2*pi_/L_sub, 0);
    }
}

complex SecondOrderWave::E_super()
{
    return delta_*pow2C(K0_super_)*wave1_->D()/wave1_->TF_[0]/wave2_->TF_[0]/pow(Omega_super_,3);
}

complex SecondOrderWave::E_sub()
{
    if(Regular_)
    {
        return complex(0,0);
    }
    else
    {
        return delta_*pow2C(K0_sub_)*wave1_->D()/wave1_->TF_[0]/wave2_->TF_[0]/pow(Omega_sub_,3);
    }
}

complex SecondOrderWave::Hjl_super(bool normal,int j,int l)
{
    if(normal) //normal
    {
        complex k1j = wave1_->k_[j];
        complex k2l = wave2_->k_[l];
        complex OOmega = {wave1_->Omega()*wave2_->Omega(), 0};
        complex k2DivOmega_j = pow2C(k1j)/wave1_->Omega();
        complex k2DivOmega_l = pow2C(k2l)/wave2_->Omega();
        complex para_1 = {0.5*(pow(wave1_->Omega(),3) + pow(wave2_->Omega(),3)), 0.0};
        return Omega_super_*(OOmega - pow(g_,2)*k1j*k2l/OOmega) + para_1 - 0.5*pow(g_,2)*(k2DivOmega_j + k2DivOmega_l);
    }
    else  // reverse, exchange wave1 and wave2
    {
        complex k1j = wave2_->k_[j];
        complex k2l = wave1_->k_[l];
        complex OOmega = {wave1_->Omega()*wave2_->Omega(), 0};
        complex k2DivOmega_j = pow2C(k1j)/wave2_->Omega();
        complex k2DivOmega_l = pow2C(k2l)/wave1_->Omega();
        complex para_1 = {0.5*(pow(wave2_->Omega(),3) + pow(wave1_->Omega(),3)), 0.0};
        return Omega_super_*(OOmega - pow(g_,2)*k1j*k2l/OOmega) + para_1 - 0.5*pow(g_,2)*(k2DivOmega_j + k2DivOmega_l);
    }
}
complex SecondOrderWave::Hjl_sub(bool normal,int j,int l)
{
    if(Regular_)
    {
        return {0.0,0.0};
    }
    else
    {
        if(normal)
        {
            complex k1j = wave1_->k_[j];
            complex k2l = wave2_->k_[l];
            complex OOmega = {wave1_->Omega()*wave2_->Omega(), 0};
            complex k2DivOmega_j = pow2C(k1j)/wave1_->Omega();
            complex k2DivOmega_l = pow2C(k2l)/wave2_->Omega();
            complex para_1 = {0.5*(pow(wave1_->Omega(),3) - pow(wave2_->Omega(),3)), 0.0};
            return Omega_sub_*(-OOmega - pow(g_,2)*k1j*(k2l).conjugate()/OOmega) + para_1 - 0.5*pow(g_,2)*(k2DivOmega_j - k2DivOmega_l);
        }
        else
        {
            complex k1j = wave2_->k_[j];
            complex k2l = wave1_->k_[l];
            complex OOmega = {wave1_->Omega()*wave2_->Omega(), 0};
            complex k2DivOmega_j = pow2C(k1j)/wave2_->Omega();
            complex k2DivOmega_l = pow2C(k2l)/wave1_->Omega();
            complex para_1 = {0.5*(pow(wave2_->Omega(),3) - pow(wave1_->Omega(),3)), 0.0};
            return -Omega_sub_*(-OOmega - pow(g_,2)*k1j*(k2l).conjugate()/OOmega) + para_1 - 0.5*pow(g_,2)*(k2DivOmega_j - k2DivOmega_l);
        }
    }
}

complex SecondOrderWave::Djl_super(int i, int j)
{
    complex k_super = wave1_->k_[i] + wave2_->k_[j];
    return g_*(k_super)*tanhC(k_super*wave1_->D())-complex(pow(Omega_super_,2),0.0);
}

complex SecondOrderWave::Djl_sub(int i, int j)
{
    complex k_sub = wave1_->k_[i] - (wave2_->k_[j]).conjugate();
    return g_*(k_sub)*tanhC(k_sub*wave1_->D())-complex(pow(Omega_sub_,2),0.0);
}

complex SecondOrderWave::Ljl_super(int i, int j)
{
    complex kk = wave1_->k_[i]*wave2_->k_[j];
    complex OOmega = {wave1_->Omega()*wave2_->Omega(), 0};
    scalar Omega2_super = pow(wave1_->Omega(),2) + pow(wave2_->Omega(),2);
    return 0.5*(pow(g_,2)*kk/OOmega - OOmega - complex(Omega2_super,0.0));
}

complex SecondOrderWave::Ljl_sub(int i, int j)
{
    complex kk = wave1_->k_[i]*(wave2_->k_[j]).conjugate();
    complex OOmega = {wave1_->Omega()*wave2_->Omega(), 0};
    scalar Omega2_super = pow(wave1_->Omega(),2) + pow(wave2_->Omega(),2);
    return 0.5*(pow(g_,2)*kk/OOmega + OOmega - complex(Omega2_super,0.0));
}

complex SecondOrderWave::Gjl_super(int i, int j)
{
    return delta_/g_*(Omega_super_*Hjl_super(true,i,j)/Djl_super(i,j)-Ljl_super(i,j));
}

complex SecondOrderWave::Gjl_sub(int i, int j)
{
    return delta_/g_*(Omega_sub_*Hjl_sub(true,i,j)/Djl_sub(i,j)-Ljl_sub(i,j));
}

complex SecondOrderWave::F11h_super()
{
    complex k0_super = wave1_->k_[0] + wave2_->k_[0];
    return E_super_*wave1_->TF_[0]*wave2_->TF_[0]*k0_super*Hjl_super(true,0,0)/(pow2C(k0_super) - pow2C(K0_super_));
}

complex SecondOrderWave::F11h_sub()
{
    if(Regular_)
    {
        return {0.0,0.0};
    }
    else
    {
        complex k0_sub = wave1_->k_[0] - wave2_->k_[0];
        return E_sub_*wave1_->TF_[0]*wave2_->TF_[0]*k0_sub*Hjl_sub(true,0,0)/(pow2C(k0_sub) - pow2C(K0_sub_));
    }
}

complex SecondOrderWave::F12h_super()
{
    return -E_super_*(
                      0.5*g_/wave1_->Omega()*wave1_->TF_[0]*pow2C(wave1_->k_[0])/(pow2C(wave1_->k_[0])-pow2C(K0_super_))*(pow(wave1_->Omega(),2)-pow(Omega_super_,2))
                    + 0.5*g_/wave2_->Omega()*wave2_->TF_[0]*pow2C(wave2_->k_[0])/(pow2C(wave2_->k_[0])-pow2C(K0_super_))*(pow(wave2_->Omega(),2)-pow(Omega_super_,2))
                    );
}

complex SecondOrderWave::F12h_sub()
{
    if(Regular_)
    {
        return {0.0,0.0};
    }
    else
    {
        complex deltaSmall(1e-20, 0);
        return E_sub_*(
                      0.5*g_/wave1_->Omega()*wave1_->TF_[0]*pow2C(wave1_->k_[0])/(pow2C(wave1_->k_[0])-pow2C(K0_sub_)+deltaSmall)*(pow(wave1_->Omega(),2)-pow(Omega_sub_,2))
                    + 0.5*g_/wave2_->Omega()*wave1_->TF_[0]*pow2C(wave2_->k_[0])/(pow2C(wave2_->k_[0])-pow2C(K0_sub_)+deltaSmall)*(pow(wave2_->Omega(),2)-pow(Omega_sub_,2))
                    );
    }
}

complex SecondOrderWave::F13h_super()
{
    complex UnitImag(0,1);
    complex obj = {0,0};
    for(int i=1;i<51;i++)
    {
        complex obj_1 = wave1_->TF_[i]*wave2_->TF_[0];
        complex obj_1_reverse = wave2_->TF_[i]*wave1_->TF_[0];
        complex obj_2 = pow2C(pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_super_))-4.0*pow2C(wave1_->k_[i])*pow2C(wave2_->k_[0]);
        complex obj_2_reverse = pow2C(pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_super_))-4.0*pow2C(wave2_->k_[i])*pow2C(wave1_->k_[0]);
        complex obj_3 = wave1_->k_[i]*(pow2C(wave1_->k_[i])-pow2C(wave2_->k_[0])-pow2C(K0_super_))*Hjl_super(true,i,0).Re();
        complex obj_3_reverse = wave2_->k_[i]*(pow2C(wave2_->k_[i])-pow2C(wave1_->k_[0])-pow2C(K0_super_))*Hjl_super(false,i,0).Re();
        complex obj_4 = wave2_->k_[0]*(-pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_super_))*Hjl_super(true,i,0).Im()*UnitImag;
        complex obj_4_reverse = wave1_->k_[0]*(-pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_super_))*Hjl_super(false,i,0).Im()*UnitImag;
        obj += E_super_*obj_1*(obj_3 + obj_4)/obj_2 + E_super_*obj_1_reverse*(obj_3_reverse + obj_4_reverse)/obj_2_reverse;
    }
    return obj;
}

complex SecondOrderWave::F13h_sub()
{
    if(Regular_)
    {
        return {0.0,0.0};
    }
    else
    {
        complex UnitImag(0,1);
        complex obj = {0,0};
        for(int i=1;i<51;i++)
        {
            complex obj_1 = wave1_->TF_[i]*wave2_->TF_[0];
            complex obj_1_reverse = wave2_->TF_[i]*wave1_->TF_[0];
            complex obj_2 = pow2C(pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_sub_))-4.0*pow2C(wave1_->k_[i])*pow2C(wave2_->k_[0]);
            complex obj_2_reverse = pow2C(pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_sub_))-4.0*pow2C(wave2_->k_[i])*pow2C(wave1_->k_[0]);
            complex obj_3 = wave1_->k_[i]*(pow2C(wave1_->k_[i])-pow2C(wave2_->k_[0])-pow2C(K0_sub_))*Hjl_sub(true,i,0).Re();
            complex obj_3_reverse = wave2_->k_[i]*(pow2C(wave2_->k_[i])-pow2C(wave1_->k_[0])-pow2C(K0_sub_))*Hjl_sub(false,i,0).Re();
            complex obj_4 = wave2_->k_[0]*(-pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_sub_))*Hjl_sub(true,i,0).Im()*UnitImag;
            complex obj_4_reverse = wave1_->k_[0]*(-pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_sub_))*Hjl_sub(false,i,0).Im()*UnitImag;
            obj += E_sub_*obj_1*(obj_3 - obj_4)/obj_2 - E_sub_*obj_1_reverse*(obj_3_reverse - obj_4_reverse)/obj_2_reverse;
        }
        return obj;
    }
}

complex SecondOrderWave::iF22h_super()
{
    complex obj = {0,0};
    for(int i=1;i<51;i++)
    {
        for(int j=1;j<51;j++)
        {
            complex obj_1 = wave1_->k_[i] + wave2_->k_[j];
            obj += wave1_->TF_[i]*wave2_->TF_[j]*obj_1/(pow2C(obj_1)-pow2C(K0_super_))*Hjl_super(true,i,j);
        }
    }
    return E_super_*obj;
}

complex SecondOrderWave::iF22h_sub()
{
    if(Regular_)
    {
        return {0.0,0.0};
    }
    else
    {
        complex obj = {0,0};
        for(int i=1;i<51;i++)
        {
            for(int j=1;j<51;j++)
            {
                complex obj_1 = wave1_->k_[i] - (wave2_->k_[j]).conjugate();
                obj += wave1_->TF_[i]*(wave2_->TF_[j]).conjugate()*obj_1/(pow2C(obj_1)-pow2C(K0_sub_))*Hjl_sub(true,i,j);
            }
        }
        return E_sub_*obj;
    }
}

complex SecondOrderWave::iF23h_super()
{
    // first half part
    complex obj_fir = {0,0};
    for(int i=1;i<51;i++)
    {
        obj_fir += wave1_->TF_[i]*pow2C(wave1_->k_[i])/(pow2C(wave1_->k_[i])-pow2C(K0_super_))*(pow(wave1_->Omega(),2)-pow(Omega_super_,2));
    }
    obj_fir *= 0.5*g_/wave1_->Omega();
    // second half part
    complex obj_sec = {0,0};
    for(int i=1;i<51;i++)
    {
        obj_sec += wave2_->TF_[i]*pow2C(wave2_->k_[i])/(pow2C(wave2_->k_[i])-pow2C(K0_super_))*(pow(wave2_->Omega(),2)-pow(Omega_super_,2));
    }
    obj_sec *= 0.5*g_/wave2_->Omega();
    return -E_super_*(obj_fir + obj_sec);
}

complex SecondOrderWave::iF23h_sub()
{
    if(Regular_)
    {
        return {0.0,0.0};
    }
    else
    {
        // first half part
        complex obj_fir = {0,0};
        for(int i=1;i<51;i++)
        {
            obj_fir += wave1_->TF_[i]*pow2C(wave1_->k_[i])/(pow2C(wave1_->k_[i])-pow2C(K0_sub_))*(pow(wave1_->Omega(),2)-pow(Omega_sub_,2));
        }
        obj_fir *= 0.5*g_/wave1_->Omega();
        // second half part
        complex obj_sec = {0,0};
        for(int i=1;i<51;i++)
        {
            obj_sec += wave2_->TF_[i]*pow2C(wave2_->k_[i])/(pow2C(wave2_->k_[i])-pow2C(K0_sub_))*(pow(wave2_->Omega(),2)-pow(Omega_sub_,2));
        }
        obj_sec *= 0.5*g_/wave2_->Omega();
        return E_sub_*(obj_fir - obj_sec);
    }
}

complex SecondOrderWave::iF24h_super()
{
    complex UnitImag(0,1);
    complex obj = {0,0};
    for(int i=1;i<51;i++)
    {
        complex obj_1 = wave1_->TF_[i]*wave2_->TF_[0];
        complex obj_1_reverse = wave2_->TF_[i]*wave1_->TF_[0];
        complex obj_2 = pow2C(pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_super_))-4.0*pow2C(wave1_->k_[i])*pow2C(wave2_->k_[0]);
        complex obj_2_reverse = pow2C(pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_super_))-4.0*pow2C(wave2_->k_[i])*pow2C(wave1_->k_[0]);
        complex obj_3 = wave1_->k_[i]*(pow2C(wave1_->k_[i])-pow2C(wave2_->k_[0])-pow2C(K0_super_))*Hjl_super(true,i,0).Im()*UnitImag;
        complex obj_3_reverse = wave2_->k_[i]*(pow2C(wave2_->k_[i])-pow2C(wave1_->k_[0])-pow2C(K0_super_))*Hjl_super(false,i,0).Im()*UnitImag;
        complex obj_4 = wave2_->k_[0]*(-pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_super_))*Hjl_super(true,i,0).Re();
        complex obj_4_reverse = wave1_->k_[0]*(-pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_super_))*Hjl_super(false,i,0).Re();
        obj += E_super_*obj_1*(obj_3 + obj_4)/obj_2 + E_super_*obj_1_reverse*(obj_3_reverse + obj_4_reverse)/obj_2_reverse;
    }
    return obj;
}

complex SecondOrderWave::iF24h_sub()
{
    if(Regular_)
    {
        return {0.0,0.0};
    }
    else
    {
        complex UnitImag(0,1);
        complex obj = {0,0};
        for(int i=1;i<51;i++)
        {
            complex obj_1 = wave1_->TF_[i]*wave2_->TF_[0];
            complex obj_1_reverse = wave2_->TF_[i]*wave1_->TF_[0];
            complex obj_2 = pow2C(pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_sub_))-4.0*pow2C(wave1_->k_[i])*pow2C(wave2_->k_[0]);
            complex obj_2_reverse = pow2C(pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_sub_))-4.0*pow2C(wave2_->k_[i])*pow2C(wave1_->k_[0]);
            complex obj_3 = wave1_->k_[i]*(pow2C(wave1_->k_[i])-pow2C(wave2_->k_[0])-pow2C(K0_sub_))*Hjl_sub(true,i,0).Im()*UnitImag;
            complex obj_3_reverse = wave2_->k_[i]*(pow2C(wave2_->k_[i])-pow2C(wave1_->k_[0])-pow2C(K0_sub_))*Hjl_sub(false,i,0).Im()*UnitImag;
            complex obj_4 = wave2_->k_[0]*(-pow2C(wave1_->k_[i])+pow2C(wave2_->k_[0])-pow2C(K0_sub_))*Hjl_sub(true,i,0).Re();
            complex obj_4_reverse = wave1_->k_[0]*(-pow2C(wave2_->k_[i])+pow2C(wave1_->k_[0])-pow2C(K0_sub_))*Hjl_sub(false,i,0).Re();
            obj += E_sub_*obj_1*(obj_3 - obj_4)/obj_2 - E_sub_*obj_1_reverse*(obj_3_reverse - obj_4_reverse)/obj_2_reverse;
        }
        return obj;
    }
}
// regular wave: phase is 0
scalar SecondOrderWave::eta2_super(scalar t)
{
    complex AA = wave1_->A(t)*wave2_->A(t);
    return (Gjl_super(0,0)*AA).Re();
}

scalar SecondOrderWave::eta2_sub(scalar t)
{
    if(Regular_)
    {
        return 0.0;
    }
    else
    {
        complex AA = wave1_->A(t)*(wave2_->A(t)).conjugate();
        return (Gjl_sub(0,0)*AA).Re();
    }
}

scalar SecondOrderWave::X2_super(scalar t)
{
    complex ImagUnit(0,1);
    return (-ImagUnit*TF_final_super_*wave1_->A(t)*wave2_->A(t)/wave1_->D()).Re();
}

scalar SecondOrderWave::X2_sub(scalar t)
{
    if(Regular_)
    {
        return 0.0;
    }
    else
    {
        complex ImagUnit(0,1);
        return (-ImagUnit*TF_final_sub_*wave1_->A(t)*(wave2_->A(t)).conjugate()/wave1_->D()).Re();
    }
}

// Focus wave: phase is -k(super/sub)*x + k(super/sub)*xf - omega(super/sub)*tf
// x is the position, 0 for paddle, any values for surface elevation
complex SecondOrderWave::eta2_super
(
    scalar phi,
    scalar xf,
    scalar tf,
    scalar t
)
{
    complex AA = wave1_->A(t)*wave2_->A(t);
    scalar k_super = (wave1_->k_[0] + wave2_->k_[0]).Re();
    scalar phi_super = k_super*xf - Omega_super_*tf - phi;
    complex obj(cos(phi_super),sin(phi_super));
    return (AA*obj);
}

complex SecondOrderWave::eta2_sub
(
    scalar phi,
    scalar xf,
    scalar tf,
    scalar t
)
{
    if(Regular_)
    {
        complex ImagZero(0,0);
        return ImagZero;
    }
    else
    {
        complex AA = wave1_->A(t)*(wave2_->A(t)).conjugate();
        scalar k_sub = (wave1_->k_[0] - wave2_->k_[0]).Re();
        scalar phi_sub = k_sub*xf - Omega_sub_*tf -phi;
        complex obj(cos(phi_sub),sin(phi_sub));
        return (AA*obj);
    }
}

complex SecondOrderWave::X2_super
(
    scalar phi,
    scalar xf,
    scalar tf,
    scalar t
)
{
    complex ImagUnit(0,1);
    scalar k_super = (wave1_->k_[0] + wave2_->k_[0]).Re();
    scalar phi_super = k_super*xf - Omega_super_*tf - phi;
    complex obj(cos(phi_super),sin(phi_super));
    return (-ImagUnit*wave1_->A(t)*wave2_->A(t)/wave1_->D()*obj);
}

complex SecondOrderWave::X2_sub
(
    scalar phi,
    scalar xf,
    scalar tf,
    scalar t
)
{
    if(Regular_)
    {
        complex ImagZero(0,0);
        return ImagZero;
    }
    else
    {
        complex ImagUnit(0,1);
        scalar k_sub = (wave1_->k_[0] - wave2_->k_[0]).Re();
        scalar phi_sub = k_sub*xf - Omega_sub_*tf -phi;
        complex obj(cos(phi_sub),sin(phi_sub));
        return (-ImagUnit*wave1_->A(t)*(wave2_->A(t)).conjugate()/wave1_->D()*obj);
    }
}
/*
void SecondOrderWave::transfer(const SecondOrderWave& waveInput)
{
    waveComponent* wave1_input = waveInput.waveComponent_1();
    waveComponent* wave2_input = waveInput.waveComponent_2();
    this->wave1_ = wave1_input;
    this->wave2_ = wave2_input;
    this->Omega_super_ = wave1_input->Omega() + wave2_input->Omega();
    this->Omega_sub_ = wave1_input->Omega() - wave2_input->Omega();
    this->Regular_ = !Omega_sub_;
}*/

} //end of namespace
