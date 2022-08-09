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
#ifndef WAVECOMPONENT_H_INCLUDED
#define WAVECOMPONENT_H_INCLUDED
#include "IOdictionary.H"
#include "complex.H"
namespace Foam
{
class waveComponent
{
    // gravity acceleration
    scalar const g_ = 9.81;
    // pi
    scalar const pi_ = 3.1415926535897932384626;
    // wave height
    scalar H_;
    // water depth
    scalar D_;
    // wave period
    scalar T_;
    // circular frequency
    scalar Omega_;
    // wave length
    scalar L_;

    // calculate wave number
    void Cal_k();
    // bisection for solving imaginary k
    scalar k_im_bisection(scalar, scalar, scalar);
    // calculate transfer function
    void Cal_TF();

    // complex exponential function
    inline complex expC(const complex& z)
    {
        return exp(z.Re())*complex(cos(z.Im()), sin(z.Im()));
    }
    // complex sinh
    inline complex sinhC(const complex& z)
    {
        return 0.5*(expC(z) - expC(-z));
    }

public:
    // constructor, H,D,T
    waveComponent
    (
        scalar waveHeight,
        scalar waterDepth,
        scalar wavePeriod,
        bool light = false
    );
    // destructor
    ~waveComponent(){};

    // array of wave number
    // the first component k[0] is real wave number and
    // the other 50 components k[1]-k[50] are imaginary wave numbers
    List<complex> k_;

    // array of transfer function, corresponding to k array
    List<complex> TF_;

    // complex amplitude
    complex A(scalar t);

    // calculate wave length
    scalar CalWaveLength(scalar, scalar);

    // access
    // water depth
    const scalar D(){return D_;};
    // circular frequency
    const scalar Omega(){return Omega_;};
    // wave period
    const scalar T(){return T_;};
	// phase velocity
    const scalar C(){return L_/T_;};
};
}// end of namespace
#endif // WAVECOMPONENT_H_INCLUDED
