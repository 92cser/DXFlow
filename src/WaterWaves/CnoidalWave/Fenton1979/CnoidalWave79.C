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

#include "CnoidalWave79.H"
#include "addToRunTimeSelectionTable.H"
#include "DynamicList.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(CnoidalWaveS, 0);
    addToRunTimeSelectionTable(waterWaves, CnoidalWaveS, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

Foam::waterWaveModels::CnoidalWaveS::CnoidalWaveS
(	
	const IOdictionary& dict
)
:
	waterWaves(dict),
	UseT(WWParameters_.lookupOrDefault<bool>("UseT", true)),
	waveHeight_(readScalar(WWParameters_.lookup("WaveHeight"))),	
	waveLength_(0.0),
	wavePeriod_(0.0),
	phaseVelocity_(0.0),
	m_(0.5),
	mApproximate_(0.0),
	alpha_(0.0),
	troughWaterDepth_(0.0),
	Epsl_(waveHeight_/waterDepth_),
	Epsl0_(0.0),
	position_(0.0)
{
	InitializationAll();
	CalculationAll();
}

void Foam::waterWaveModels::CnoidalWaveS::InitializationAll()
{
	if(UseT) // use T calculate L
	{
		wavePeriod_ = readScalar(WWParameters_.lookup("WavePeriod"));
	}
	else     // use L calculate T
	{
		waveLength_ = readScalar(WWParameters_.lookup("WaveLength"));
	}
	InitializeWaveLengthCoeff_();
	InitializeWavePeriodCoeff_();
	InitializeTWDCoeff_();
	InitializeAlphaCoeff_();
	InitializePVCoeff_();
	InitializeEtaCoeff_();
}

void Foam::waterWaveModels::CnoidalWaveS::CalculationAll()
{
	Cal_m();
	Cal_TWD();
	Cal_alpha();
	Cal_PV();
}

Foam::scalar Foam::waterWaveModels::CnoidalWaveS::K(scalar m)
{
    scalar K_ = 0.0;
    scalar deltaTheta_ = 0.5*pi_/9999;
    for(int i = 0; i < 10000; i++)
    {
        scalar theta_i = i*deltaTheta_;
        K_ += pow(1 - m*pow(sin(theta_i),2), -0.5)*deltaTheta_;
    }
    return K_;
}

Foam::scalar Foam::waterWaveModels::CnoidalWaveS::E(double m)
{
    scalar E_ = 0.0;
    scalar deltaTheta_ = 0.5*pi_/9999;
    for(int i = 0; i < 10000; i++)
    {
        scalar theta_i = i*deltaTheta_;
        E_ += pow(1 - m*pow(sin(theta_i),2), 0.5)*deltaTheta_;
    }
    return E_;
}

void Foam::waterWaveModels::CnoidalWaveS::InitializeWaveLengthCoeff_()
{
	scalar LCoeffs_ijk[5][5][5] =
        //(0,0,0) - (0,4,4): 25 components
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        //(1,0,0) - (1,4,4): 25 components
         1.25,-1.5,0,0,0,
         -0.625,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        //(2,0,0) - (2,4,4): 25 components
         -0.46875,0.125,0.375,0,0,
         0.46875,-0.0625,0,0,0,
         -0.16406,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,
        //(3,0,0) - (3,4,4): 25 components
         1.01556,-0.91938,0.21875,0.0625,0,
         -1.5233,0.91983,-0.10983,0,0,
         0.73241,-0.06391,0,0,0,
         -0.11232,0,0,0,0,
         0,0,0,0,0,
        //(4,0,0) - (4,4,4): 25 components
         -2.79984,3.66490,-1.48453,0.20313,0.02344,
         -5.59969,-5.49735,1.48453,-0.10156,0,
         4.07395,3.04491,-0.10465,0,0,
         1.27410,-0.60623,0,0,0,
         -0.05493,0,0,0,0
        };
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			for(int k = 0; k < 5; k++)
			{
				LCoeffs_[i][j][k] = LCoeffs_ijk[i][j][k];
				//Info << LCoeffs_[i][j][k] << "(" << i << j << k <<"), ";
			}
		}
	}
}

void Foam::waterWaveModels::CnoidalWaveS::InitializeWavePeriodCoeff_()
{
    scalar TCoeffs_ijk[5][5][5] =
        //(0,0,0) - (0,4,4): 25 components
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        //(1,0,0) - (1,4,4): 25 components
         0.25,0,0,0,0,
         -0.125,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        //(2,0,0) - (2,4,4): 25 components
         0.01458,-0.108333,1,0,0,
         -0.01458,0.54167,0,0,0,
         -0.07656,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,
        //(3,0,0) - (3,4,4): 25 components
         0.36121,2.50417,-4.5,2.0,0,
         -0.54182,-2.50417,2.25,0,0,
         0.41216,0.33229,0,0,0,
         -0.11578,0,0,0,0,
         0,0,0,0,0,
        //(4,0,0) - (4,4,4): 25 components
         -1.86885,-4.22859,15.19111,-13.66667,4,
         3.73770,6.34288,-15.19111,6.8333,0,
         -2.73031,-1.88433,2.69111,0,0,
         0.86147,-0.11498,0,0,0,
         -0.07582,0,0,0,0
        };
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            for(int k = 0; k < 5; k++)
            {
                TCoeffs_[i][j][k] = TCoeffs_ijk[i][j][k];
                //cout << TCoeffs_[i][j][k] << "(" << i << j << k <<"), ";
            }
        }
    }
}

void Foam::waterWaveModels::CnoidalWaveS::InitializeTWDCoeff_()
{
    scalar TWDCoeffs_ijk[6][6][6] =
        //(0,0,0) - (0,5,5): 36 components
         {0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
        //(1,0,0) - (1,5,5): 36 components
         1,-1,0,0,0,0,
         -1,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
        //(2,0,0) - (2,5,5): 36 components
         -0.5,0.5,0,0,0,0,
         0.5,-0.25,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
        //(3,0,0) - (3,5,5): 36 components
         0.665,-1.165,0.5,0,0,0,
         -0.9975,1.165,-0.25,0,0,0,
         0.3325,-0.04,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
        //(4,0,0) - (4,5,5): 36 components
         -1.62667,3.20667,-2.08,0.5,0,0,
         3.25333,-4.81,2.08,-0.25,0,0,
         -2.454,2.17633,-0.1425,0,0,0,
         0.82733,-0.2865,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
        //(5,0,0) - (5,5,5): 36 components
         4.86659,-10.74409,8.6225,-3.245,0.5,0,
         -12.16647,21.48818,-12.93375,3.245,-0.25,0,
         -11.79929,-16.00776,6.09025,-0.30750,0,0,
         -5.53247,5.26368,-0.88950,0,0,0,
         1.03306,-0.20555,0,0,0,0,
         0,0,0,0,0,0,
        };
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 6; j++)
        {
            for(int k = 0; k < 6; k++)
            {
                TWDCoeffs_[i][j][k] = TWDCoeffs_ijk[i][j][k];
                //cout << TWDCoeffs_[i][j][k] << "(" << i << j << k <<"), ";
            }
        }
    }
}

void Foam::waterWaveModels::CnoidalWaveS::InitializeAlphaCoeff_()
{
    scalar AlphaCoeffs_ij[5][5] =
    {
		// (0,0) - (0,4): 5 components
		0,0,0,0,0,
		// (1,0) - (1,4): 5 components
		0.25,-0.875,0,0,0,
		// (2,0) - (2,4): 5 components
		0.03125,-0.34375,0.86719,0,0,
		// (3,0) - (3,4): 5 components
		-0.37743,0.51146,0.13743,-0.833,0,
		// (4,0) - (4,4): 5 components
		0.20322,0.44278,-1.38945,0.54282,0.76773
    };
    for (int i = 0;i < 5;i++)
    {
        for (int j = 0;j < 5;j++)
        {
            AlphaCoeffs_[i][j] = AlphaCoeffs_ij[i][j];
            //cout << AlphaCoeffs_[i][j] << "("<< i<<j<<"), " <<endl;
        }
    }
}

void Foam::waterWaveModels::CnoidalWaveS::InitializePVCoeff_()
{
    scalar PVCoeffs_ijk[6][6][2] =
		//(0,0,0) - (0,5,1): 12 components
		 {0,0,
		  0,0,
		  0,0,
		  0,0,
		  0,0,
		  0,0,
		//(1,0,0) - (1,5,1): 12 components
		 0.5,-1,
		 0,0,
		 0,0,
		 0,0,
		 0,0,
		 0,0,
		//(2,0,0) - (2,5,1): 12 components
		 -0.10833,0.33333,
		 -0.01667,0.08333,
		 -0.025,0,
		 0,0,
		 0,0,
		 0,0,
		//(3,0,0) - (3,5,1): 12 components
		 -0.17190,0.09333,
		 0.33911,-0.34333,
		 -0.16006,0.21833,
		 0.04643,0,
		 0,0,
		 0,0,
		//(4,0,0) - (4,5,1): 12 components
		 0.02097,0.37690,
		 0.17293,-0.68202,
		 -0.56238,1.04889,
		 0.39861,-0.56668,
		 -0.08531,0,
		 0,0,
		//(5,0,0) - (5,5,1): 12 components
		 0.11046,-0.94038,
		 -0.31285,1.22117,
		 -0.11262,0.35314,
		 0.91605,-1.75325,
		 -0.73881,1.00619,
		 0.15763,0
		};
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 6; j++)
        {
            for(int k = 0; k < 2; k++)
            {
                PVCoeffs_[i][j][k] = PVCoeffs_ijk[i][j][k];
                //cout << PVCoeffs_[i][j][k] << "(" << i << j << k <<"), ";
            }
        }
    }
}

void Foam::waterWaveModels::CnoidalWaveS::InitializeEtaCoeff_()
{
    scalar EtaCoeffs_ijk[6][6][6] =
        //(0,0,0) - (0,5,5): 36 components
         {0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
        //(1,0,0) - (1,5,5): 36 components
         0,0,0,0,0,0,
         0,1,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
        //(2,0,0) - (2,5,5): 36 components
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,-0.75,0.75,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
        //(3,0,0) - (3,5,5): 36 components
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,-0.7625,0.7625,0,0,0,
         0,1.38750,-2.65,1.2625,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
        //(4,0,0) - (4,5,5): 36 components
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,-0.80533,0.80533,0,0,0,
         0,2.48904,-4.33146,1.84242,0,0,
         0,-3.05188,7.40646,-6.52546,2.17088,0,
         0,0,0,0,0,0,
        //(5,0,0) - (5,5,5): 36 components
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0.43643,-0.43643,0,0,0,
         0,1.92280,-4.66167,2.73888,0,0,
         0,-7.04588,17.45561,-15.31697,4.90723,0,
         0,6.54722,-19.80887,25.34187,16.32709,4.24687,
        };
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 6; j++)
        {
            for(int k = 0; k < 6; k++)
            {
                EtaCoeffs_[i][j][k] = EtaCoeffs_ijk[i][j][k];
                //cout << EtaCoeffs_[i][j][k] << "(" << i << j << k <<"), ";
            }
        }
    }
}

void Foam::waterWaveModels::CnoidalWaveS::Cal_m()
{
    scalar mLeft_ = 0;
    scalar mRight_ = 1;
    if(UseT) // wave period is known
    {
        scalar wavePeriod0_ = 0;
        while(mag(wavePeriod0_ - wavePeriod_) > 0.001*wavePeriod_)
        {
            // start from m = 0.5
            if(wavePeriod0_ < wavePeriod_)
            {
                mLeft_ = m_;
                m_ = 0.5*(m_ + mRight_);
            }
            else if(wavePeriod0_ > wavePeriod_)
            {
                mRight_ = m_;
                m_ = 0.5*(m_ + mLeft_);
            }
            //Info << "m = " << m_;
            scalar para_1 = 4*K(m_)*pow((3*Epsl_/m_),-0.5);
            // initialization
            wavePeriod0_ *= 0;
            for(int i = 1; i < 5; i++)
            {
                for(int j = 0; j < i+1; j++)
                {
                    for(int k = 0; k < i+1; k++)
                    {
                        wavePeriod0_ += pow((Epsl_/m_),i)*pow(m_,j)*pow(E(m_)/K(m_),k)*TCoeffs_[i][j][k];
                    }
                }
            }
            wavePeriod0_ = para_1*(1+wavePeriod0_)*pow(g_/waterDepth_, -0.5);
            //Info << ", Period = " <<  wavePeriod0_ << endl;
        }
		// approximation 
		mApproximate_ = (m_ > 0.96)?1:m_; 

        // calculate wave length after m is calculated
        scalar para_2 = 4*K(m_)*pow((3*Epsl_/mApproximate_),-0.5);
        for(int i = 1; i < 5; i++)
        {
            for(int j = 0; j < i+1; j++)
            {
                for(int k = 0; k < i+1; k++)
                {
                    waveLength_ += pow((Epsl_/mApproximate_),i)*pow(mApproximate_,j)*pow(E(m_)/K(m_),k)*LCoeffs_[i][j][k];
                }
            }
        }
        waveLength_ = para_2*(1 + waveLength_)*waterDepth_;
    }
    else// wave length is known
    {
        scalar waveLength0_ = 0;
        while(mag(waveLength0_ - waveLength_) > 0.001*waveLength_)
        {
            // start from m = 0.5
            if(waveLength0_ < waveLength_)
            {
                mLeft_ = m_;
                m_ = 0.5*(m_ + mRight_);
            }
            else if(waveLength0_ > waveLength_)
            {
                mRight_ = m_;
                m_ = 0.5*(m_ + mLeft_);
            }
            //cout << "m = " << m_;
            scalar para_1 = 4*K(m_)*pow((3*Epsl_/m_),-0.5);
            // initialization
            waveLength0_ *= 0;
            for(int i = 1; i < 5; i++)
            {
                for(int j = 0; j < i+1; j++)
                {
                    for(int k = 0; k < i+1; k++)
                    {
                        waveLength0_ += pow((Epsl_/m_),i)*pow(m_,j)*pow(E(m_)/K(m_),k)*LCoeffs_[i][j][k];
                    }
                }
            }
            waveLength0_ = para_1*(1+waveLength0_)*waterDepth_;
            //cout << ", waveLength_ = " <<  waveLength0_ << endl;
        }
		// approximation 
		mApproximate_ = (m_ > 0.96)?1:m_; 

        // calculate wave period after m is calculated
        scalar para_2 = 4*K(m_)*pow((3*Epsl_/mApproximate_),-0.5);
        for(int i = 1; i < 5; i++)
        {
            for(int j = 0; j < i+1; j++)
            {
                for(int k = 0; k < i+1; k++)
                {
                    wavePeriod_ += pow((Epsl_/mApproximate_),i)*pow(mApproximate_,j)*pow(E(m_)/K(m_),k)*TCoeffs_[i][j][k];
                }
            }
        }
        wavePeriod_ = para_2*(1+wavePeriod_)*pow(g_/waterDepth_, -0.5);
    }
}

void Foam::waterWaveModels::CnoidalWaveS::Cal_TWD()
{
    for(int i = 1; i < 6; i++)
    {
        for(int j = 0; j < i+1; j++)
        {
            for(int k = 0; k < i+1; k++)
            {
                troughWaterDepth_ += pow((Epsl_/mApproximate_),i)*pow(mApproximate_,j)*pow(E(m_)/K(m_),k)*TWDCoeffs_[i][j][k];
                //cout << "(" << i << j << k <<"), ";
            }
        }
    }
    troughWaterDepth_ = (troughWaterDepth_ + 1)*waterDepth_;
// calculate epsl0
    Epsl0_ = waveHeight_/troughWaterDepth_;
}

void Foam::waterWaveModels::CnoidalWaveS::Cal_alpha()
{
    scalar para_1 = Epsl0_/mApproximate_;
    for (int i = 1;i < 5;i++)
    {
        for (int j = 0;j < 5;j++)
        {
            alpha_ += pow(para_1,i)*pow(mApproximate_,j)*AlphaCoeffs_[i][j];
            //Info << AlphaCoeffs_[i][j] << "("<< i<<j<<"), ";
        }
    }
    alpha_ = (alpha_ + 1)* pow(0.75*para_1, 0.5);
    //Info << "Alpha = " << alpha_ << endl;
}

void Foam::waterWaveModels::CnoidalWaveS::Cal_PV()
{
    scalar para_1 = Epsl0_/mApproximate_;
    for(int i = 1; i < 6; i++)
    {
        for(int j = 0; j < i+1; j++)
        {
            for(int k = 0; k < 2; k++)
            {
                phaseVelocity_ += pow(para_1,i)*pow(mApproximate_,j)*pow(E(m_)/K(m_),k)*PVCoeffs_[i][j][k];
                //Info << "(" << i << j << k <<"), ";
            }
        }
    }
    phaseVelocity_ = (phaseVelocity_ + 1)*pow(g_*troughWaterDepth_, 0.5);
    //Info << "phaseVelocity_ = " << phaseVelocity_ << endl;
}

Foam::scalar Foam::waterWaveModels::CnoidalWaveS::AGM_Cn
(
    scalar a0, 
    scalar b0, 
    scalar c0, 
    scalar t
)
{
    scalar alphaX = alpha_*(position_ - phaseVelocity_*t)/troughWaterDepth_;
    //vector<double> aN; //from a0 to aN
    DynamicList<scalar> aN;
    //vector<double> cNByaN;
    DynamicList<scalar> cNByaN;
    // initialization
    scalar a_ = 0.0;
    scalar b_ = 0.0;
    scalar c_ = 1.0;
    scalar a_old = a0;
    scalar b_old = b0;
    scalar c_old = c0;
    aN.append(a_old);
    cNByaN.append(c_old/a_old);
    while (fabs(c_) > 1e-6)
    {
        a_ = 0.5*(a_old + b_old);
        b_ = pow(a_old*b_old, 0.5);
        c_ = 0.5*(a_old - b_old);
        a_old = a_;
        b_old = b_;
        c_old = c_;
        aN.append(a_old);
        cNByaN.append(c_old/a_old);
    }
    DynamicList<scalar>phiN(aN);
    for(int i = 0; i < phiN.size();i++)
    {
        int location = phiN.size() - 1 - i;
        //Info <<"phiN[location]pre = " << phiN[location] << endl;
        if(i==0)
        {
            phiN[location] = pow(2, location)*aN[location]*alphaX;
        }
        else
        {
            scalar para_1 = cNByaN[location + 1]*sin(phiN[location + 1]);
            phiN[location] = 0.5*(asin(para_1) + phiN[location + 1]);
        }
        //Info <<"phiN[location]after = " << phiN[location] << endl;
    }
    return cos(phiN[0]);
}

Foam::scalar Foam::waterWaveModels::CnoidalWaveS::eta(scalar t)
{
    scalar eta = 0;
    scalar para_1 = Epsl0_/mApproximate_;
    for(int i = 1; i < 6; i++)
    {
        for(int j = 0; j < i+1; j++)
        {
            for(int k = 1; k < j+1; k++)
            {
                eta += pow(para_1,i)*pow(mApproximate_,j)*pow(AGM_Cn(1.0,pow(1-m_,0.5),pow(m_,0.5),t),2*k)*EtaCoeffs_[i][j][k];
                //cout << "(" << i << j << k <<"), ";
            }
        }
    }
    eta = (eta + 1)*troughWaterDepth_ -waterDepth_;
    // slowly start factor
    scalar SSF_ = (wavePeriod_ > 0.5*t) ? pow(0.5*t/wavePeriod_,2.0) : 1.0;
    return eta*SSF_;
}


void Foam::waterWaveModels::CnoidalWaveS::PrintWaveProperties()
{
    Info << "-----------------------------------------------" << "\n"
         << "Wave properties of the Fifth-order Cnoidal Wave (only accurate for small amplitude waves)" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Wave Height (m) = " << waveHeight_ << "\n"
         << " Wave Period (s) = " << wavePeriod_ << "\n"
         << " Wave Length (m) = " << waveLength_ << "\n"
         << " Minimum Water Depth (m) = " << troughWaterDepth_ << "\n"
         << " Phase Velocity (m/s) = " << phaseVelocity_ << "\n"
         << " Actual Module Number = " << m_ << "\n"
         << " Approximate Module Number = " << mApproximate_ << "\n"
         << " K = " << K(m_) << "\n"
         << " Straining Factor = " << alpha_ << "\n"
         << "-----------------------------------------------" << endl;
}

Foam::word Foam::waterWaveModels::CnoidalWaveS::name()
{
	return "Cnoidal wave";
}

Foam::scalar Foam::waterWaveModels::CnoidalWaveS::Cal_C()
{
	return phaseVelocity_;
}

Foam::scalar Foam::waterWaveModels::CnoidalWaveS::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
    scalar etaForPiston_ = etaForPiston(t);
    scalar PistonVelocity_ = phaseVelocity_*etaForPiston_ / (waterDepth_+etaForPiston_) -0.1*position_;
    position_ += PistonVelocity_*deltaT;
    Info << "TAndX = " << t << ", " << position_ << endl;
    return position_;
}

Foam::scalar Foam::waterWaveModels::CnoidalWaveS::etaForPiston(scalar t)
{
    scalar eta = 0;
    scalar para_1 = Epsl0_/mApproximate_;
    for(int i = 1; i < 6; i++)
    {
        for(int j = 0; j < i+1; j++)
        {
            for(int k = 1; k < j+1; k++)
            {
                eta += pow(para_1,i)*pow(mApproximate_,j)*pow(AGM_Cn(1.0,pow(1-m_,0.5),pow(m_,0.5),t),2*k)*EtaCoeffs_[i][j][k];
                //cout << "(" << i << j << k <<"), ";
            }
        }
    }
    eta = (eta + 1)*troughWaterDepth_ -waterDepth_;
    Info << "TAndEta = " << t << ", " << eta << endl;
    return eta;
}
