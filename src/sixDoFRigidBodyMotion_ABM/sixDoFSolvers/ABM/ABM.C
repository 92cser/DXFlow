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

#include "ABM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{
    defineTypeNameAndDebug(ABM, 0);
    addToRunTimeSelectionTable(sixDoFSolver, ABM, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::ABM::ABM
(
    const dictionary& dict,
    sixDoFRigidBodyMotion& body
)
:
    sixDoFSolver(dict, body),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 0.5)),
    beta_
    (
        max
        (
            0.25*sqr(gamma_ + 0.5),
            dict.lookupOrDefault<scalar>("beta", 0.25)
        )
    ),
    ADRF_a(0.0),
    ADRF_tau(0.0),
    VL_(8,vector::zero),
    nIteration_(0.0),
    timeForPISO_(dict.getOrDefault<scalar>("timeForPISO", 0.0))
{
    // wether the PISO-PIMPLE alogrithm is used
    if(timeForPISO_ > 0)
    {
        dynamicOuterLoop_ = true;
        Info << "Using the PISO-PIMPLE alogrithm... "<< nl
             << "Note that the timeForPISO in dynamicMeshDict should be equal to the one in fvSolution." << endl;
    }
    else
    {
        dynamicOuterLoop_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::ABM::~ABM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFSolvers::ABM::solve
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0,
    scalar currT
)
{
    // use PIMPLE algorithm, strong coulping
    if(currT >= timeForPISO_)
    {
        dynamicOuterLoop_ = false;
    }

    //Info << "deltaT0 = " << deltaT0 << endl;
    //Info << "a00() = " << a00() << endl;

    // PIMPLE only: reached at the initial two time steps
    // PISO-PIMPLE: reached before currT > timeForPISO_
    if (a00() == vector(0,0,0) || dynamicOuterLoop_)
    {
        //if(firstIter)
        {
            Info << "Using Newmark method for weak-coupling solving" << endl;

            // Update the linear acceleration and torque
            updateAcceleration(fGlobal, tauGlobal);

            // Correct linear velocity
            v() =
                    tConstraints()
              & (v0() + aDamp()*deltaT*(gamma_*a() + (1 - gamma_)*a0()));

            // Correct angular momentum
            pi() =
                    rConstraints()
              & (pi0() + aDamp()*deltaT*(gamma_*tau() + (1 - gamma_)*tau0()));

            // Correct position
            centreOfRotation() =
                    centreOfRotation0()
              + (
                            tConstraints()
                      & (
                                    deltaT*v0()
                              + aDamp()*sqr(deltaT)*(beta_*a() + (0.5 - beta_)*a0())
                            )
                    );

            // Correct orientation
            vector piDeltaT =
                    rConstraints()
              & (
                            deltaT*pi0()
                      + aDamp()*sqr(deltaT)*(beta_*tau() + (0.5 - beta_)*tau0())
                    );
            Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
            Q() = Qpi.first();
        }

        //else
        //{
            //Info << "Newmark method is only looped once for each time step" << endl;
        //}

    }
    else // use ABM method
    {
        if (firstIter) // Predictor (without relaxtion)
        {
            Info << "Using ABM method solve the rigid body motion (the predictor)" << endl;
            nIteration_ = 0;
            // Correct linear velocity
            v() =
                    tConstraints()
              & (v0() + update_explicit(deltaT, deltaT/deltaT0, a0(), a00()));
            // Correct angular momentum
            pi() =
                    rConstraints()
              & (pi0() + update_explicit(deltaT, deltaT/deltaT0, tau0(), tau00()));
            // Correct position
            centreOfRotation() =
                    centreOfRotation0()
              + (
                            tConstraints()
                      & (
                                    update_implicit
                                    (
                                            deltaT,
                                            deltaT/deltaT0,
                                            v(),
                                            v0(),
                                            v00()
                                    )
                            )
                    );

            // Correct orientation
            vector piDeltaT =
                    rConstraints()
              & (
                                    update_implicit
                                    (
                                            deltaT,
                                            deltaT/deltaT0,
                                            pi(),
                                            pi0(),
                                            pi00()
                                    )
                    );
            Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
            Q() = Qpi.first();

            Info << "Relative tolerance of linear acceleration =  " << vector::zero << endl;
            Info << "Relative tolerance of torque =  " << vector::zero << endl;

            Info << "Restraint force (ABM first iteration): " << vector::zero << endl;
        }
        else // corrector
        {
            // Update the number of iterations (this number equals to the number of corrector)
            nIteration_ += 1;

            Info << "Using ABM method solve the rigid body motion (the " << nIteration_ << "-th corrector)" << endl;

            // Update the linear acceleration and torque based on the result of the last iteration
             updateAcceleration_ABM(fGlobal, tauGlobal);

            // Store values of the last iteration and calculate the Aitken's dynamic under-relaxtion factor
            // the first corrector, now we only have a00, a0 and aPredictor. cannot relaxation.
            if (nIteration_ == 1)
            {
                // original

                VL_[0] = a00();
                VL_[1] = a0();
                VL_[2] = a();    //value obtained by predictor

                VL_[4] = tau00();
                VL_[5] = tau0();
                VL_[6] = tau();  //value obtained by predictor

                // since no relaxation, the error is evalueated using a0 and aPredictor
                Info << "Relative tolerance of linear acceleration =  " << CalResiduals(VL_[1] ,VL_[2]) << endl;
                Info << "Relative tolerance of torque =  " << CalResiduals(VL_[5], VL_[6]) << endl;

                // no relaxation, set to 0 to facilitate postprocessing
                ADRF_a = 0;
                ADRF_tau = 0;
            }
            // the second or subsequent iteration
            // Under-relaxation starts at the second cycle after predictor (The third cycle of PIMPLE outer loop)
            else
            {
                VL_[3] = a();    //value obtained by last corrector(before relaxation)
                VL_[7] = tau();  //value obtained by last corrector(before relaxation)

                // calculate relative error using the new value and the value of last iteration(after relaxation)
                Info << "Relative tolerance of linear acceleration =  " << CalResiduals(VL_[2] ,VL_[3]) << endl;
                Info << "Relative tolerance of torque =  " << CalResiduals(VL_[6], VL_[7]) << endl;

                if (nIteration_ == 2)
                {
                    // calculate the initial value of under-relaxtion factor
                    // using a00, a0, aPridictor and a_1stCorrector(before relaxation)
                    ADRF_a = ADRF_intial(VL_[0], VL_[1], VL_[2], VL_[3]);
                    ADRF_tau = ADRF_intial(VL_[4], VL_[5], VL_[6], VL_[7]);
                }
                else
                {
                    // calculate the under-relaxtion factor
                    ADRF_a *= -ADRF(VL_[0], VL_[1], VL_[2], VL_[3]);
                    ADRF_tau *= -ADRF(VL_[4], VL_[5], VL_[6], VL_[7]);
                }

                 // 0.1 - 1
                ADRF_a = max(0.1, min(ADRF_a, 1));
                ADRF_tau = max(0.1, min(ADRF_tau, 1));

                // Under-relaxtion
                a() =  VL_[2] + ADRF_a*(a()- VL_[2]);
                tau() = VL_[6] + ADRF_tau*(tau()-VL_[6]);

                // store the value after relaxation
                VL_[0] = VL_[2];
                VL_[1] = VL_[3];
                VL_[2] = a();
                VL_[4] = VL_[6];
                VL_[5] = VL_[7];
                VL_[6] = tau();
            }

            Info << "Aitken's dynamic under-relaxation factor for linear acceleration = " << ADRF_a << endl;
            Info << "Aitken's dynamic under-relaxation factor for torque = " << ADRF_tau << endl;

            // Correct linear velocity
            v() =
                    tConstraints()
              & (v0() + update_implicit(deltaT, deltaT/deltaT0, a(), a0(), a00()));
            // Correct angular momentum
            pi() =
                    rConstraints()
              & (pi0() + update_implicit(deltaT, deltaT/deltaT0, tau(), tau0(), tau00()));
            // Correct position
            centreOfRotation() =
                    centreOfRotation0()
              + (
                            tConstraints()
                      & (
                                    update_implicit
                                    (
                                            deltaT,
                                            deltaT/deltaT0,
                                            v(),
                                            v0(),
                                            v00()
                                    )
                            )
                    );
            // Correct orientation
            vector piDeltaT =
                    rConstraints()
              & (
                                    update_implicit
                                    (
                                            deltaT,
                                            deltaT/deltaT0,
                                            pi(),
                                            pi0(),
                                            pi00()
                                    )
                    );
            Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
            Q() = Qpi.first();
        }
    }
}

Foam::vector Foam::sixDoFSolvers::ABM::update_explicit
(
    scalar deltaT,
    scalar deltaTBydeltaT0,
    vector value0,
    vector value00
)
{
    return 0.5*deltaT*
            (
                (2.0 + deltaTBydeltaT0) * value0 - deltaTBydeltaT0 * value00
            );
}

Foam::vector Foam::sixDoFSolvers::ABM::update_implicit
(
    scalar deltaT,
    scalar deltaTBydeltaT0,
    vector value_new,
    vector value0,
    vector value00
)
{
    return 0.0625*deltaT/deltaTBydeltaT0*
            (
                (1.0 + 8.0*deltaTBydeltaT0) * value_new + (7.0*deltaTBydeltaT0 - 1.0) * value0 + deltaTBydeltaT0 * value00
            );
}
/*
inline Foam::scalar Foam::sixDoFSolvers::ABM::CalResiduals
(
    vector a_pre,
    vector a_new
)
{
    return Foam::mag(a_new-a_pre) / (Foam::mag(a_pre) + 1e-20);
}*/

inline Foam::vector Foam::sixDoFSolvers::ABM::CalResiduals
(
    vector a_pre,
    vector a_new
)
{
    //return Foam::mag(a_new-a_pre) / (Foam::mag(a_pre) + 1e-20);

    // 2022.5.8
    vector theResidual = vector::zero;
    forAll(theResidual,i)
    {
        theResidual[i] = mag(a_new[i] - a_pre[i])/mag(a_pre[i]+1e-20);
    }
    return theResidual;
}

Foam::scalar Foam::sixDoFSolvers::ABM::ADRF_intial
(
    vector& a_i, // a00
    vector& a_j, // a0
    vector& a_m, // aPredictor
    vector& a_n  // a_1stCorrector
)
{
    vector firstError = a_i - a_m;
    vector secondError = a_n - a_j;
    return ADRF_update(firstError, secondError);
}

Foam::scalar Foam::sixDoFSolvers::ABM::ADRF
(
    vector& a_i,
    vector& a_j,
    vector& a_m,
    vector& a_n
)
{
    vector firstError = a_j - a_i;
    vector secondError = a_n - a_m;
    return ADRF_update(firstError, secondError);
}

Foam::scalar Foam::sixDoFSolvers::ABM::ADRF_update
(
    vector& firstError,
    vector& secondError
)
{
     return firstError
                    & ( secondError -  firstError)
                    / (Foam::pow(Foam::mag(secondError -  firstError), 2.0) + 1e-20);
}
// ************************************************************************* //
