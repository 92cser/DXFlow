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
#ifndef SECONDORDERPISTON_H_INCLUDED
#define SECONDORDERPISTON_H_INCLUDED

#include "WaveComponent.H"
namespace Foam
{
class SecondOrderWave
{
    // gravity acceleration
    scalar const g_ = 9.81;
    // pi
    scalar const pi_ = 3.1415926535897932384626;

    // pointers of pair of wave components
    waveComponent* wave1_;
    waveComponent* wave2_;

    // sum of circular frequencies
    scalar Omega_super_;
    // subtraction of circular frequencies
    scalar Omega_sub_;

    // switch of subharmonics
    // when regular wave is the target wave, items of subharmonics vanish
    bool Regular_;

    // coefficient
    scalar delta_;

    // real solution of dispersion equation(superharmonics)
    complex K0_super_;
    // real solution of dispersion equation(subharmonics)
    complex K0_sub_;

    // calculate K0
    complex K0_super();
    complex K0_sub();

    // Second-order items
    // E
    complex E_super_;
    complex E_sub_;
    // calculate E
    complex E_super();
    complex E_sub();

    // H, with number j,l
    complex Hjl_super
        (
            bool normal, // normal or reverse
            int j,
            int l
         );
    complex Hjl_sub
        (
            bool normal,
            int j,
            int l
         );

    // D, with number j,l
    complex Djl_super(int i, int j);
    complex Djl_sub(int i, int j);

    // L, with number j,l
    complex Ljl_super(int i, int j);
    complex Ljl_sub(int i, int j);

    // G, with number j,l
    complex Gjl_super(int i, int j);
    complex Gjl_sub(int i, int j);

    // F11h
    complex F11h_super_;
    complex F11h_sub_;
    // calculate F11h
    complex F11h_super();
    complex F11h_sub();

    // F12h
    complex F12h_super_;
    complex F12h_sub_;
    // calculate F12h
    complex F12h_super();
    complex F12h_sub();

    // F13h
    complex F13h_super_;
    complex F13h_sub_;
    // calculate F13h
    complex F13h_super();
    complex F13h_sub();

    // iF22h
    complex iF22h_super_;
    complex iF22h_sub_;
    // calculate F13h
    complex iF22h_super();
    complex iF22h_sub();

    // iF23h
    complex iF23h_super_;
    complex iF23h_sub_;
    // calculate F13h
    complex iF23h_super();
    complex iF23h_sub();

    // iF24h
    complex iF24h_super_;
    complex iF24h_sub_;
    // calculate F13h
    complex iF24h_super();
    complex iF24h_sub();

    // final transfer function
    complex TF_final_super_;
    complex TF_final_sub_;

    inline complex pow2C(const complex& z)
    {
        return complex
        (
            pow(z.Re(),2)-pow(z.Im(),2),
            2*z.Re()*z.Im()
        );
    }

    // complex exponential function
    inline complex expC(const complex& z)
    {
        return exp(z.Re())*complex(cos(z.Im()), sin(z.Im()));
    }

    // complex tanh
    inline complex tanhC(const complex& z)
    {
        return (expC(z) - expC(-z))/(expC(z) + expC(-z));
    }

public:
    // default constructor
    //SecondOrderWave();

    // constructor
    SecondOrderWave(waveComponent&, waveComponent&);

    // constructor, used for calculating 2nd surface and paddle motion during run time
    SecondOrderWave(waveComponent&, waveComponent&, word);
    
    // destructor
    ~SecondOrderWave(){};

    // regular wave
    // second-order surface elevation
    scalar eta2_super(scalar t);
    scalar eta2_sub(scalar t);

    // second-order paddle movement
    scalar X2_super(scalar t);
    scalar X2_sub(scalar t);

    // focus wave or irregular wave
    // second-order surface elevation
    // using complex form for accelerating the loop
    complex eta2_super
    (
        scalar phi, // random phase (0 for focus wave)
        scalar xf,  // focus position (0 for irregular wave)
        scalar tf,  // focus time (0 for irregular wave)
        scalar t    // current time
    );
    complex eta2_sub
    (
        scalar phi,
        scalar xf,
        scalar tf,
        scalar t
    );
    // second-order paddle movement
    complex X2_super
    (
        scalar phi,
        scalar xf,
        scalar tf,
        scalar t
    );
    complex X2_sub
    (
        scalar phi,
        scalar xf,
        scalar tf,
        scalar t
    );

    // access
    waveComponent* waveComponent_1() const {return wave1_;};
    waveComponent* waveComponent_2() const {return wave2_;};
    complex TF_final_super(){return TF_final_super_;};
    complex TF_final_sub(){return TF_final_sub_;};
/*
    const SecondOrderWave operator=(const SecondOrderWave& input)
    {
        this->transfer(input);
        return *this;
    }

    // transfer
    void transfer(const SecondOrderWave& waveInput);*/
};
} //end of namespace
#endif // SECONDORDERPISTON_H_INCLUDED
