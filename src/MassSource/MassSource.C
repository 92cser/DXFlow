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
#include "MassSource.H"

// Constructor
Foam::massSourceWaveMaker::massSourceWaveMaker
(
    const IOdictionary& dict,
    volScalarField& massSource,
    volScalarField& waveSourceField,
    volVectorField& velocityField
)
:
    markField_(massSource),
    mesh_(markField_.mesh()),
    waveSourceValue_(waveSourceField),
    velocityValue_(velocityField),
    wave_(waterWaveModels::waterWaves::New(dict)),
    useMassSource_(dict.isDict("MassSourceWaveMaker")),
    phaseVelocity_(0.0),
    waveSurface_(0.0),
    areaOfSource_(0.0),
    waveEnhanceFactor_(1.0),
    uneven_(false),
    distributionField_(markField_),
    gamma_(0.0),
    adjustU_(false),//(dict.subDict("MassSourceWaveMaker").lookupOrDefault<bool>("AdjustU", false)),
    adjustUndesiredW_(false),
    Zt_(0.0),
    Zb_(0.0),
    Zc_(0.0),
    UDampingCoeff_(0.0),
    UDampingIndex_(0.0),
    uAdjustField_(markField_),
    AdjustmentScheme_("Default")
    //adjustSymmetry_(false),
    //x_left_(),
    //x_right_()

{
    if(useMassSource_)
    {
        Info << "Mass source wavemaker is adopted for " << wave_->name() << endl;

        // print
        wave_->PrintWaveProperties();

        // update phase velocity
        CalculatePhaseVelocity();

        // read basic parameters
        areaOfSource_ = readScalar(dict.subDict("MassSourceWaveMaker").lookup("AreaOfSource"));
        waveEnhanceFactor_ = readScalar(dict.subDict("MassSourceWaveMaker").lookup("EnhanceFactor"));

        // get top and bottom of source region
        Zt_ = readScalar(dict.subDict("MassSourceWaveMaker").lookup("Zt"));
        Zb_ = readScalar(dict.subDict("MassSourceWaveMaker").lookup("Zb"));

        // use uneven mass source or not
        gamma_ = readScalar(dict.subDict("MassSourceWaveMaker").lookup("Gamma"));

        if(gamma_ == 1e-6)
        {
            Info << "Uniform wavemaker is utilized" << endl;
        }
        else
        {
            Info << "Uneven wavemaker is utilized with the redistribution factor of "<< gamma_ << endl;
            uneven_ = true;
        }
        CalDistributionField();

        // adjust undesired vertical velocity or not
        adjustU_ = dict.subDict("MassSourceWaveMaker").isDict("AdjustU");

        if(adjustU_)  // adjust velocity inside the source region
        {
            // adjust vertical velocity or not
            adjustUndesiredW_ = dict.subDict("MassSourceWaveMaker").subDict("AdjustU").lookupOrDefault<bool>("AdjustUndesiredW", false);
            if(adjustUndesiredW_)
            {
                // read
                UDampingCoeff_ = readScalar(dict.subDict("MassSourceWaveMaker").subDict("AdjustU").lookup("AdjustHeight"));
                UDampingIndex_ = readScalar(dict.subDict("MassSourceWaveMaker").subDict("AdjustU").lookup("IndexOfAdjustment"));
                AdjustmentScheme_ = dict.subDict("MassSourceWaveMaker").subDict("AdjustU").lookupOrDefault<word>("AdjustmentScheme","tri");
                // calculate
                CalUAdjustField();
                Info << "The undesired veritical velocity will be adjusted using " << AdjustmentScheme_ << " scheme" << endl;
            }
            // ensure symmetry or not
            /*adjustSymmetry_ = dict.subDict("MassSourceWaveMaker").subDict("AdjustU").lookupOrDefault<bool>("EnsureSymmetry", false);
            if(adjustSymmetry_)
            {
                    // calculate
                    CalDynamicLists();
                    Info << "The symmetry of velocity field will be ensured"<< endl;
            }*/
        }
    };
}

// Member functions
void Foam::massSourceWaveMaker::CalculateWaveSurface(scalar t)
{
    waveSurface_ = wave_->eta(t);
}

void Foam::massSourceWaveMaker::CalculatePhaseVelocity()
{
    phaseVelocity_ = wave_->Cal_C();
}

void Foam::massSourceWaveMaker::Update(scalar t)
{			
    CalculateWaveSurface(t);
    scalar typicalMassSource = 2.0 * waveEnhanceFactor_* phaseVelocity_ * waveSurface_ / areaOfSource_;

    forAll(markField_, i)
    {
        if(markField_[i])
        {
            if(uneven_) //uneven mass source
            {
                waveSourceValue_[i] = typicalMassSource* distributionField_[i];
            }
            else // uniform mass source
            {
                waveSourceValue_[i] = typicalMassSource;
            }
        }
    }
}

void Foam::massSourceWaveMaker::CalDistributionField()
{
    Zc_ = gamma_/(exp(gamma_)-1);
    forAll(distributionField_, i)
    {
        if(distributionField_[i]) //inside mass source
        {
            scalar weight = Zc_*exp(gamma_*(mesh_.C()[i].component(2) - Zb_)/(Zt_ - Zb_));
            distributionField_[i] = weight;
        }
    }
}

void Foam::massSourceWaveMaker::CalUAdjustField()
{
    scalar adjustLength = Zs() - Zb_;
    forAll(uAdjustField_,i)
    {
        if (uAdjustField_[i]) //inside mass source
        {
            if(mesh_.C()[i].component(2) < Zs())
            {
                scalar verticalDis = mesh_.C()[i].component(2) - Zb_;
                scalar z_relative = verticalDis/adjustLength;
                if(AdjustmentScheme_ == "power")
                {
                    uAdjustField_[i] = Foam::pow(z_relative, UDampingIndex_);
                }
                if(AdjustmentScheme_ == "tri")//trigonometic
                {
                    uAdjustField_[i] = Foam::pow(sin(0.5*3.1415926*z_relative),UDampingIndex_);
                }
                if(AdjustmentScheme_ == "exp")//exponential
                {
                    uAdjustField_[i] = (exp(Foam::pow(z_relative, UDampingIndex_)) - 1.0)/(exp(1.0)-1.0);
                }
            }
        }
    }
}

Foam::scalar Foam::massSourceWaveMaker::Zs()
{
    return UDampingCoeff_*(Zt_-Zb_) + Zb_;
}
/*
void Foam::massSourceWaveMaker::CalDynamicLists()
{
	// cell label (inside the source region)
	DynamicList<label> CellLabel(1);
	// x position of cell center (inside the source region)
	DynamicList<scalar> CellCenterX(1);
	// z position of cell center (inside the source region)
	DynamicList<scalar> CellCenterZ(1);

	forAll(markField_, i)
	{
		if(markField_[i])
		{
			CellLabel.append(i);
			CellCenterX.append(mesh_.C()[i].component(0));
			CellCenterZ.append(mesh_.C()[i].component(2));
		}
	}

	if(CellLabel.size())
	{
	// calculate the cell size inside source region
	// the cell should be uniform
		Pout << "Pstream::myProcNo() = " << Pstream::myProcNo() << ", Cell number of source region = " << CellLabel.size() << endl;

		// face list 
		const faceList & FL = mesh_.faces();
		// point list
		const pointField & PL = mesh_.points();

		const cell & CL = mesh_.cells()[CellLabel[0]];

		// pLabels : labels of 8 points of the current cell 
		labelList pLabels(CL.labels(FL));

		// empty pointField 
		pointField pLocal(pLabels.size(), vector::zero);

		forAll (pLabels, pointi)
		{
			  pLocal[pointi] = PL[pLabels[pointi]];
		}

		scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
		//scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
		scalar zDim = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
		Pout << "xDim = " << xDim << endl;
		Pout << "zDim = " << zDim << endl;

		scalar NCell_x_ = (max(CellCenterX) - min(CellCenterX))/xDim + 1;
		scalar NCell_z_ = (max(CellCenterZ) - min(CellCenterZ))/zDim + 1;
		Pout << "nx, nz = " << NCell_x_ << ", " <<NCell_z_ << endl;

		scalar X_center = 0.0;
		// calculate x center
		forAll(CellCenterX, i)
		{
			X_center += CellCenterX[i];
		}
		X_center /= (NCell_x_*NCell_z_);
		Pout << "X_center = " << X_center << endl;

		forAll(CellLabel, i)
		{
			scalar xi = mesh_.C()[CellLabel[i]].component(0);
			//scalar yi = mesh_.C()[CellLabel[i]].component(1);
			scalar zi = mesh_.C()[CellLabel[i]].component(2);

			if(xi < X_center)
			{
				// x of the symmetrical cell
				scalar xi_p = 2.0*X_center - xi;
				// store 
				x_left_.append(CellLabel[i]);
				// find 
				forAll(CellLabel, j)
				{
					scalar xj = mesh_.C()[CellLabel[j]].component(0);
					//scalar yj = mesh_.C()[CellLabel[j]].component(1);
					scalar zj = mesh_.C()[CellLabel[j]].component(2);
					if(xj > X_center)
					{
						if( mag(xj-xi_p) < 1e-10  && mag(zj-zi)<1e-10)//(yi == yj && zi == zj))
						{
							x_right_.append(CellLabel[j]); //store
						}
					}
				}
			}
		}
		Pout << x_left_.size() << ", " << x_right_.size() <<endl;
	}
}
*/
void Foam::massSourceWaveMaker::AdjustVelocity()
{
    /*if(adjustSymmetry_)
    {
            Info << "Ensure the symmetry of velocity field"<< endl;
            forAll(x_left_, i)
            {
                    scalar L = velocityValue_[x_left_[i]].component(0)/(mag(velocityValue_[x_left_[i]].component(0))+1e-10);
                    scalar R = velocityValue_[x_right_[i]].component(0)/(mag(velocityValue_[x_right_[i]].component(0))+1e-10);
                    scalar Ux_adjust = 0.5*(mag(velocityValue_[x_left_[i]].component(0)) + mag(velocityValue_[x_right_[i]].component(0)));
                    velocityValue_[x_left_[i]].component(0) = L*Ux_adjust;
                    velocityValue_[x_right_[i]].component(0) = R*Ux_adjust;

                    scalar Uz_adjust = 0.5*(velocityValue_[x_left_[i]].component(2) + velocityValue_[x_right_[i]].component(2));
                    velocityValue_[x_left_[i]].component(2) = Uz_adjust;
                    velocityValue_[x_right_[i]].component(2) = Uz_adjust;
            }
    }*/
    if(adjustUndesiredW_)
    {
        Info << "Adjust undesired vertical velocity inside mass source region" << endl;
        forAll(markField_,i)
        {
            if(markField_[i] && mesh_.C()[i].component(2) < Zs())
            {
                velocityValue_[i] *= uAdjustField_[i];
            }
        }
    }
		
}




