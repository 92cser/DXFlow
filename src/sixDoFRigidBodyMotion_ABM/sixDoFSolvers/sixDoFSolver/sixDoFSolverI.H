/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

inline Foam::point& Foam::sixDoFSolver::centreOfRotation()
{
    return body_.motionState_.centreOfRotation();
}

inline Foam::tensor& Foam::sixDoFSolver::Q()
{
    return body_.motionState_.Q();
}

inline Foam::vector& Foam::sixDoFSolver::v()
{
    return body_.motionState_.v();
}

inline Foam::vector& Foam::sixDoFSolver::a()
{
    return body_.motionState_.a();
}

inline Foam::vector& Foam::sixDoFSolver::pi()
{
    return body_.motionState_.pi();
}

inline Foam::vector& Foam::sixDoFSolver::tau()
{
    return body_.motionState_.tau();
}


inline const Foam::point& Foam::sixDoFSolver::centreOfRotation0() const
{
    return body_.motionState0_.centreOfRotation();
}

inline const Foam::tensor& Foam::sixDoFSolver::Q0() const
{
    return body_.motionState0_.Q();
}


inline const Foam::vector& Foam::sixDoFSolver::v0() const
{
    return body_.motionState0_.v();
}


inline const Foam::vector& Foam::sixDoFSolver::a0() const
{
    return body_.motionState0_.a();
}


inline const Foam::vector& Foam::sixDoFSolver::pi0() const
{
    return body_.motionState0_.pi();
}


inline const Foam::vector& Foam::sixDoFSolver::tau0() const
{
    return body_.motionState0_.tau();
}

//added 2019.2.20, ABM method
inline const Foam::point& Foam::sixDoFSolver::centreOfRotation00() const
{
    return body_.motionState00_.centreOfRotation();
}

inline const Foam::tensor& Foam::sixDoFSolver::Q00() const
{
    return body_.motionState00_.Q();
}


inline const Foam::vector& Foam::sixDoFSolver::v00() const
{
    return body_.motionState00_.v();
}


inline const Foam::vector& Foam::sixDoFSolver::a00() const
{
    return body_.motionState00_.a();
}


inline const Foam::vector& Foam::sixDoFSolver::pi00() const
{
    return body_.motionState00_.pi();
}


inline const Foam::vector& Foam::sixDoFSolver::tau00() const
{
    return body_.motionState00_.tau();
}

inline Foam::scalar Foam::sixDoFSolver::aDamp() const
{
    return body_.aDamp_;
}

inline Foam::tensor Foam::sixDoFSolver::tConstraints() const
{
    return body_.tConstraints_;
}

inline Foam::tensor Foam::sixDoFSolver::rConstraints() const
{
    return body_.rConstraints_;
}

//- Apply rotation tensors to Q0 for the given torque (pi) and deltaT
//  and return the rotated Q and pi as a tuple
inline Foam::Tuple2<Foam::tensor, Foam::vector> Foam::sixDoFSolver::rotate
(
    const tensor& Q0,
    const vector& pi,
    const scalar deltaT
) const
{
    return body_.rotate(Q0, pi, deltaT);
}

//- Update and relax accelerations from the force and torque
inline void Foam::sixDoFSolver::updateAcceleration
(
    const vector& fGlobal,
    const vector& tauGlobal
)
{
    body_.updateAcceleration(fGlobal, tauGlobal);
}

//- Update accelerations from the force and torque
inline void Foam::sixDoFSolver::updateAcceleration_ABM
(
    const vector& fGlobal,
    const vector& tauGlobal
)
{
    body_.updateAcceleration_ABM(fGlobal, tauGlobal);
}


// ************************************************************************* //
