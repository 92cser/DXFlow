/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 Dalian Ocean University / National Marine Environment Monitoring Center / Dalian University of Technology
    Copyright (C) 2018-2022 Ocean University of China
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

#include "waveDamping.H"
Foam::waveDamping::DampingLayer::DampingLayer
(
	const fvMesh& mesh,
	const dictionary& inputSubDict,
        const word& location,
        bool  uneven
)
:
        mesh_(mesh),
	location_(location),
        Valid_(inputSubDict.isDict(location)),
        xs_(0.0),
        xe_(0.0),
        xWidth_(0.0),
        index_(2),
        gamma0_(0.0),
        blendingField_(mesh.C().component(0)*0.0),
        ForceAndBlendingField_(mesh.C().component(0)*0.0),
        k_(0),
        Spt_(0),
        Spb_(0),
        zWidth_(0),
        zSWL_(0),
        gamma1_(0),
        d_(0)
{
        if(Valid_)
	{
            // mesh center vector
            const volVectorField& Centre = mesh.C();
            // x y z
            const scalarList x_ = Centre.component(0);
            const scalarList y_ = Centre.component(1);
            const scalarList z_ = Centre.component(2);

            // adopt uneven forcing strength
            if(uneven)
            {
                const dictionary& unevenSPLayerDict = inputSubDict.subDict("UnevenForcingStrength");
                k_ = unevenSPLayerDict.get<scalar>("waveNumber");
                Spt_ = unevenSPLayerDict.get<scalar>("Spt");
                Spb_ = unevenSPLayerDict.get<scalar>("Spb");
                zSWL_ = unevenSPLayerDict.get<scalar>("zSWL");
                zWidth_ = Spt_ - Spb_;
                d_ = zSWL_ - Spb_;
            }

            // calculate wave damping source
            if(location_ == "x+")
            {
                gamma0_ = inputSubDict.subDict(location).get<scalar>("ForcingStrength");
                xs_ = inputSubDict.subDict(location).get<scalar>("StartPosition");
                xe_ = inputSubDict.subDict(location).get<scalar>("EndPosition");
                index_ = inputSubDict.subDict(location).getOrDefault<scalar>("Index", 2.0);
                xWidth_ = mag(xe_ - xs_);

                if(uneven) //calculate gamma1
                {
                    gamma1_ = k_*gamma0_*zWidth_*sinh(k_*d_)/sinh(k_*(Spt_-zSWL_+d_));
                }
                // update the source
                forAll(blendingField_ , i)
                {
                    if(x_[i] > xs_ && x_[i] < xe_) // xs_ < xe_ because waves travel along positive x-direction
                    {
                        scalar x_inside_ = (x_[i] -  xs_)/xWidth_;
                        blendingField_[i] = (exp(pow(x_inside_, index_)) - 1) / (exp(1.0) - 1);
                        if(uneven)
                        {
                            // distribution function
                            scalar phiZ = cosh(k_*(z_[i] - zSWL_ + d_))/sinh(k_*d_);
                            if(z_[i] > (Spt_- zSWL_)) // above the Spt_, adopt uniform distribution
                            {
                                phiZ = cosh(k_*(Spt_ - zSWL_ + d_))/sinh(k_*d_);
                            }
                            ForceAndBlendingField_[i] = gamma1_* phiZ* blendingField_[i];
                        }
                        else
                        {
                            ForceAndBlendingField_[i] = gamma0_* blendingField_[i];
                        }
                    }
                }
            }
            else if(location_ == "x-")
            {
                gamma0_ = inputSubDict.subDict(location).get<scalar>("ForcingStrength");
                xs_ = inputSubDict.subDict(location).get<scalar>("StartPosition");
                xe_ = inputSubDict.subDict(location).get<scalar>("EndPosition");
                index_ = inputSubDict.subDict(location).getOrDefault<scalar>("Index", 2.0);
                xWidth_ = mag(xe_ - xs_);

                if(uneven) //calculate gamma1
                {
                    gamma1_ = k_*gamma0_*zWidth_*sinh(k_*d_)/sinh(k_*(Spt_-zSWL_+d_));
                }
                // update the source
                forAll(blendingField_ , i)
                {
                    if(x_[i] < xs_ && x_[i] > xe_) // xs_ > xe_ because waves travel along negative x-direction
                    {
                        scalar x_inside_ = (x_[i] -  xs_)/xWidth_;
                        blendingField_[i] = (exp(pow(x_inside_, index_)) - 1) / (exp(1.0) - 1);
                        if(uneven)
                        {
                            // distribution function
                            scalar phiZ = cosh(k_*(z_[i] - zSWL_ + d_))/sinh(k_*d_);
                            if(z_[i] > (Spt_- zSWL_)) // above the Spt_, adopt uniform distribution
                            {
                                phiZ = cosh(k_*(Spt_ - zSWL_ + d_))/sinh(k_*d_);
                            }
                            ForceAndBlendingField_[i] = gamma1_* phiZ* blendingField_[i];
                        }
                        else
                        {
                            ForceAndBlendingField_[i] = gamma0_* blendingField_[i];
                        }
                    }
                }
            }
            else if(location_ == "y+") // x <==> y
            {
                gamma0_ = inputSubDict.subDict(location).get<scalar>("ForcingStrength");
                scalar ys = inputSubDict.subDict(location).get<scalar>("StartPosition");
                scalar ye = inputSubDict.subDict(location).get<scalar>("EndPosition");
                index_ = inputSubDict.subDict(location).getOrDefault<scalar>("Index", 2.0);
                scalar yWidth_ = mag(ye - ys);

                if(uneven) //calculate gamma1
                {
                    gamma1_ = k_*gamma0_*zWidth_*sinh(k_*d_)/sinh(k_*(Spt_-zSWL_+d_));
                }
                // update the source
                forAll(blendingField_ , i)
                {
                    if(y_[i] > ys && y_[i] < ye) // ys < ye because waves travel along positive y-direction
                    {
                        scalar y_inside_ = (y_[i] -  xs_)/yWidth_;
                        blendingField_[i] = (exp(pow(y_inside_, index_)) - 1) / (exp(1.0) - 1);
                        if(uneven)
                        {
                            // distribution function
                            scalar phiZ = cosh(k_*(z_[i] - zSWL_ + d_))/sinh(k_*d_);
                            if(z_[i] > (Spt_- zSWL_)) // above the Spt_, adopt uniform distribution
                            {
                                phiZ = cosh(k_*(Spt_ - zSWL_ + d_))/sinh(k_*d_);
                            }
                            ForceAndBlendingField_[i] = gamma1_* phiZ* blendingField_[i];
                        }
                        else
                        {
                            ForceAndBlendingField_[i] = gamma0_* blendingField_[i];
                        }
                    }
                }
            }
            else if(location_ == "y-")
            {
                gamma0_ = inputSubDict.subDict(location).get<scalar>("ForcingStrength");
                scalar ys = inputSubDict.subDict(location).get<scalar>("StartPosition");
                scalar ye = inputSubDict.subDict(location).get<scalar>("EndPosition");
                index_ = inputSubDict.subDict(location).getOrDefault<scalar>("Index", 2.0);
                scalar yWidth_ = mag(ye - ys);

                if(uneven) //calculate gamma1
                {
                    gamma1_ = k_*gamma0_*zWidth_*sinh(k_*d_)/sinh(k_*(Spt_-zSWL_+d_));
                }

                // update the source
                forAll(blendingField_ , i)
                {
                    if(y_[i] < ys && y_[i] > ye) // ys > ye because waves travel along negative y-direction
                    {
                        scalar y_inside_ = (y_[i] -  ys)/yWidth_;
                        blendingField_[i] = (exp(pow(y_inside_, index_)) - 1) / (exp(1.0) - 1);
                        if(uneven)
                        {
                            // distribution function
                            scalar phiZ = cosh(k_*(z_[i] - zSWL_ + d_))/sinh(k_*d_);
                            if(z_[i] > (Spt_- zSWL_)) // above the Spt_, adopt uniform distribution
                            {
                                phiZ = cosh(k_*(Spt_ - zSWL_ + d_))/sinh(k_*d_);
                            }
                            ForceAndBlendingField_[i] = gamma1_* phiZ* blendingField_[i];
                        }
                        else
                        {
                            ForceAndBlendingField_[i] = gamma0_* blendingField_[i];
                        }
                    }
                }
            }
	}
	else
	{
            Info << "No Sponge layer is utilized at the position of " << location_
                 << endl;
	}
}

const Foam::volScalarField Foam::waveDamping::DampingLayer::blendingField()
{
	return blendingField_;
}

const Foam::volScalarField Foam::waveDamping::DampingLayer::ForceAndBlendingField()
{
        return ForceAndBlendingField_;
}

Foam::waveDamping::waveDamping
(
	const fvMesh& mesh,
	const dictionary& inputDict
)
:
	mesh_(mesh),
        alphaOrNot_(inputDict.subDict("WaveDamping").getOrDefault<bool>("AlphaSource", false)),
        unevenOrNot_(inputDict.subDict("WaveDamping").isDict("UnevenForcingStrength")),
        X_pos_
        (
            mesh_,
            inputDict.subDict("WaveDamping"),
            "x+",
            unevenOrNot_
        ),
        X_neg_
        (
            mesh_,
            inputDict.subDict("WaveDamping"),
            "x-",
            unevenOrNot_
            ),
        Y_pos_
        (
            mesh_,
            inputDict.subDict("WaveDamping"),
            "y+",
            unevenOrNot_
            ),
        Y_neg_
        (
            mesh_,
            inputDict.subDict("WaveDamping"),
            "y-",
            unevenOrNot_
        ),
	DimensionBalance_("DimensionBalance_", dimensionSet(0,-1,-1,0,0,0,0), 1.0),
        BlendingField_(X_pos_.blendingField() + X_neg_.blendingField() + Y_pos_.blendingField() + Y_neg_.blendingField()),
        DampingField_(X_pos_.ForceAndBlendingField() + X_neg_.ForceAndBlendingField() + Y_pos_.ForceAndBlendingField() + Y_neg_.ForceAndBlendingField())
{
    if(alphaOrNot_)
    {
        Info << "A wave damping source term is added in the VOF equation." << endl;
    }
    if(unevenOrNot_)
    {
        Info << "Uneven forcing strength is adopted. " << endl;
    }
}

const Foam::volScalarField Foam::waveDamping::DampingField()
{
	return DampingField_*DimensionBalance_;
}

const Foam::volScalarField Foam::waveDamping::BlendingField()
{
        return BlendingField_;
}
Foam::scalar Foam::waveDamping::AlphaSource()
{
    if(alphaOrNot_)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}
