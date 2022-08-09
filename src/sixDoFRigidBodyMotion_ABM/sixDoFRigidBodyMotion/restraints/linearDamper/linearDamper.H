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

Class
    Foam::sixDoFRigidBodyMotionRestraint::linearDamper

Description
    Linear damper class, including the method for modelling Coulomb type PTO
    damping force.

Usage
    in dynamicMeshDict
    restraints
    {
        linearSpring
        {
            sixDoFRigidBodyMotionRestraint linearDamper;
            coeff            		1.06;// 0.75*70.7*0.02
            Coulomb	  		true;
            CouForce	  		(0 0 0.223); // 0.75*lambdaOpt*Vmax
            MaxVeloctyLptoOpt 		(0 0 0.21); // maximum velocity under the condition of linear optimized PTO
            PTORelaxationWeight	  	 0.3;
        }
    }

SourceFiles
    linearDamper.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef linearDamper_H
#define linearDamper_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class linearDamper Declaration
\*---------------------------------------------------------------------------*/

class linearDamper
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        //- switch of using linear or Coulomb damping
        //- ture: use Coulomb damping
        //- false: use linear damping (by default)
        bool Cou_;

        //- Coulomb damping force (N)
        vector CouForce_;

        //- Damping force of the last iteration
        mutable vector DampingForce0_;

        //- Damping coefficient (Ns/m), for linear damping
        scalar coeff_;

        // - Maximum velocity under the condition of linear optimized PTO damping
        vector vLptoOptMax_;

        //- (0.05*vLptoOptMax_)^2
        vector VThreshold_;

        //- to avoid suddent change in PTO force, the current force might be blended with the previous iteration
        // this is the blending weight 0.5 by default
        scalar weightPTO_;

public:

    //- Runtime type information
    TypeName("linearDamper");


    // Constructors

        //- Construct from components
        linearDamper
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new linearDamper(*this)
            );
        }


    //- Destructor
    virtual ~linearDamper();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
