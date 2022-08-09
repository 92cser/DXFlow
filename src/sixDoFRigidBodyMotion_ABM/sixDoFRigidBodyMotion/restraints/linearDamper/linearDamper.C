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

#include "linearDamper.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(linearDamper, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        linearDamper,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearDamper::linearDamper
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
    coeff_(),
    Cou_(false),
    CouForce_(vector::zero),
    DampingForce0_(vector::zero),
    vLptoOptMax_(vector::zero),
    VThreshold_(vector::zero),
    weightPTO_(0)
{
    read(sDoFRBMRDict);
    // coulomb
    Cou_ = sDoFRBMRDict.getOrDefault<bool>("Coulomb", false);
    if(Cou_)
    {
        CouForce_ = sDoFRBMRDict.get<vector>("CouForce");
        vLptoOptMax_ = sDoFRBMRDict.get<vector>("MaxVeloctyLptoOpt");
        weightPTO_ = sDoFRBMRDict.getOrDefault<scalar>("PTORelaxationWeight", 0.5);
        VThreshold_ = 0.05*vLptoOptMax_;
        Info << " Adopt Coulomb type damping, MaxVeloctyLptoOpt is " << vLptoOptMax_ << ", PTORelaxationWeight is " << weightPTO_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearDamper::~linearDamper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::linearDamper::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintForce = -coeff_*motion.v();

     // Coulomb force
    if (Cou_)
    {
        // previous iteration velocity (the velocity has not been updated)
        vector v0 = motion.v();
        // pre-previous iteration velocity
        vector v00 = motion.v00();
        forAll(v0,i)
        {
            // calculate Coulumb damping force
            scalar couForce = -sign(v0[i])*mag(CouForce_[i]);
            // if velocity is too small, or if the linear damping is greater than Coulomb damping force
            // do not use Coulomb damping force, use linear damping force instead
            scalar vSquare = v0[i]*v00[i];
            if(vSquare > VThreshold_[i]*VThreshold_[i] || mag(restraintForce) > mag(couForce))
            {
                restraintForce[i] = couForce;
            }
            // avoid suddent change in damping force
            scalar error = (restraintForce[i] - DampingForce0_[i])/(min(mag(DampingForce0_[i]), mag(restraintForce[i]))+1e-20);
            if(mag(error) > weightPTO_)
            {
                // further underrelaxation
                if(mag(error) > 2*weightPTO_)
                {
                    const scalar weightPTOSquare= weightPTO_*weightPTO_;
                    restraintForce[i] = (1-weightPTOSquare)*DampingForce0_[i] + weightPTOSquare*restraintForce[i];
                }
                else // common underrelaxation
                {
                    restraintForce[i] = (1-weightPTO_)*DampingForce0_[i] + weightPTO_*restraintForce[i];
                }
            }
        }
        DampingForce0_ = restraintForce;
    }

    restraintMoment = Zero;

    if (motion.report())
    {
        Info<< " force " << restraintForce
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::linearDamper::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("coeff", coeff_);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::linearDamper::write
(
    Ostream& os
) const
{
    os.writeEntry("coeff", coeff_);
}


// ************************************************************************* //
