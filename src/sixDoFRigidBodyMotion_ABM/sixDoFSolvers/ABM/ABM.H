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

Class
    Foam::sixDoFSolvers::ABM

Description
    The class for updating the motion of floating structures, which consists
    of one predictor and several correctors.

    Reference:
    Wang, D., Dong, S., Ning, M., Incecik, A. (2020). Extended variable-time-step
    adams–bashforth–moulton method for strongly coupled fluid–structure
    interaction simulation. Ocean Engineering, 219, 108335.

Usage
    in dynamicMeshDict
    solver
    {
        type 		ABM;
    }
    The number of PIMPLE outerloop should be no less than 3.

SourceFiles
    ABM.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef ABM_H
#define ABM_H

#include "sixDoFSolver.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{

/*---------------------------------------------------------------------------*\
                           Class Newmark Declaration
\*---------------------------------------------------------------------------*/

class ABM
:
    public sixDoFSolver
{
    // Private data

        //- The first two coofficients are needee by Newmark method which is utilized in the first two time steps
        //- Coefficient for velocity integration (default: 0.5)
        const scalar gamma_;

        //- Coefficient for position and orientation integration (default: 0.25)
        const scalar beta_;

        //- The value of Aitken dynamic relax factor for acceleration
        scalar ADRF_a;

        //- The value of Aitken dynamic relax factor for torque
        scalar ADRF_tau;

        //- Number of iterations
        label nIteration_;

        //- The vectorList that stores accelerations/torques for relaxation factor calculation
        vectorList VL_;

        //- Wether the PISO-PIMPLE alogrithm is used or not
        bool dynamicOuterLoop_;

        //- If the PISO-PIMPLE alogrithm is used, read the time for PISO
        scalar timeForPISO_;

public:

    //- Runtime type information
    TypeName("ABM");

    // Constructors

        //- Construct from a dictionary and the body
        ABM
        (
            const dictionary& dict,
            sixDoFRigidBodyMotion& body
        );
       //- Construct and return a clone
        virtual autoPtr<sixDoFSolver> clone() const
        {
            return autoPtr<sixDoFSolver>::NewFrom<ABM>(dict_, body_);
        }

    //- Destructor
    virtual ~ABM();


    // Member Functions

        //- Drag coefficient
        virtual void solve
        (
            bool firstIter,
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT,
            scalar deltaT0,
            scalar currT
        );

        //- Explicit ABM method
        vector update_explicit
        (
            scalar deltaT,
            scalar deltaTBydeltaT0,
            vector value0,
            vector value00
        );
        //- Implicit ABM method
        vector update_implicit
        (
            scalar deltaT,
            scalar deltaTBydeltaT0,
            vector value_new,
            vector value0,
            vector value00
        );
        //- Residuals control
        /*inline scalar CalResiduals
        (
            vector a_pre,
            vector a_new
        );*/

        //- Residuals control (2022.5.8)
        inline vector CalResiduals
        (
            vector a_pre,
            vector a_new
        );

        //- Evaluate the initial Aitken dynamic relax factor
        scalar ADRF_intial
        (
            vector& a_i,
            vector& a_j,
            vector& a_m,
            vector& a_n
        );

        //- Calculate the Aitken dynamic relax factor
        scalar ADRF
        (
            vector& a_i,
            vector& a_j,
            vector& a_m,
            vector& a_n
        );

        //- Aitken dynamic relax factor in error vector form
        scalar ADRF_update
        (
            vector& firstError,
            vector& secondError
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sixDoFSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
