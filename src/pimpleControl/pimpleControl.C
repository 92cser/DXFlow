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

Description
    Add some codes to ensure the weak-strong coupling, that is, a non-constant
    outer loop number. This can save time in floating structure simulations when
    strong coupling is considered.

    However, it should be noted that the use of the weak-strong coupling might
    unstabilize the simulation.

Usage
    in fvSolution:
    PIMPLE
    {
        // before this time, the case is run in PISO mode
        timeForPISO		7.5;  //5T
    }
    and simultaneously in dynamicMeshDict:
    solver
    {
        type 		ABM; // Newmark;
        timeForPISO	7.5; // should be equal to the time in fvSolution
    }

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/


#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::pimpleControl::read()
{
    solutionControl::read(false);

    const dictionary pimpleDict(dict());

    solveFlow_ = pimpleDict.getOrDefault("solveFlow", true);
    nCorrPIMPLE_ = pimpleDict.getOrDefault<label>("nOuterCorrectors", 1);
    nCorrPISO_ = pimpleDict.getOrDefault<label>("nCorrectors", 1);
    SIMPLErho_ = pimpleDict.getOrDefault("SIMPLErho", false);
    turbOnFinalIterOnly_ =
        pimpleDict.getOrDefault("turbOnFinalIterOnly", true);
    finalOnLastPimpleIterOnly_ =
        pimpleDict.getOrDefault("finalOnLastPimpleIterOnly", false);
    ddtCorr_ = pimpleDict.getOrDefault("ddtCorr", true);

    // for PISO-PIMPLE combined algorithm
    const scalar timeForPISO = pimpleDict.getOrDefault<scalar>("timeForPISO", 0);
    if(timeForPISO > 0)
    {
        const scalar currT = mesh_.time().value();
        if(currT < timeForPISO) // PISO is adopted
        {
            nCorrPIMPLE_ = 1;
        }
    }

    return true;
}


bool Foam::pimpleControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if ((corr_ == 1) || residualControl_.empty() || finalIter())
    {
        return false;
    }


    const bool storeIni = this->storeInitialResiduals();

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.solverPerformanceDict();
    forAllConstIters(solverDict, iter)
    {
        const entry& solverPerfDictEntry = *iter;

        const word& fieldName = solverPerfDictEntry.keyword();
        const label fieldi = applyToField(fieldName);

        if (fieldi != -1)
        {
            Pair<scalar> residuals = maxResidual(solverPerfDictEntry);

            checked = true;

            scalar relative = 0.0;
            bool relCheck = false;

            const bool absCheck =
                (residuals.last() < residualControl_[fieldi].absTol);

            if (storeIni)
            {
                residualControl_[fieldi].initialResidual = residuals.first();
            }
            else
            {
                const scalar iniRes =
                    (residualControl_[fieldi].initialResidual + ROOTVSMALL);

                relative = residuals.last() / iniRes;
                relCheck = (relative < residualControl_[fieldi].relTol);
            }

            achieved = achieved && (absCheck || relCheck);

            if (debug)
            {
                Info<< algorithmName_ << " loop:" << endl;

                Info<< "    " << fieldName
                    << " PIMPLE iter " << corr_
                    << ": ini res = "
                    << residualControl_[fieldi].initialResidual
                    << ", abs tol = " << residuals.last()
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << ", rel tol = " << relative
                    << " (" << residualControl_[fieldi].relTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


void Foam::pimpleControl::setFirstIterFlag(const bool check, const bool force)
{
    DebugInfo
        << "corr:" << corr_
        << " corrPISO:" << corrPISO_
        << " corrNonOrtho:" << corrNonOrtho_
        << endl;

    solutionControl::setFirstIterFlag(check && corrPISO_ <= 1, force);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleControl::pimpleControl
(
    fvMesh& mesh,
    const word& dictName,
    const bool verbose
)
:
    solutionControl(mesh, dictName),
    solveFlow_(true),
    nCorrPIMPLE_(0),
    nCorrPISO_(0),
    corrPISO_(0),
    SIMPLErho_(false),
    turbOnFinalIterOnly_(true),
    finalOnLastPimpleIterOnly_(false),
    ddtCorr_(true),
    converged_(false)
{
    read();

    if (verbose)
    {
        Info<< nl << algorithmName_;

        if (nCorrPIMPLE_ > 1)
        {
            if (residualControl_.empty())
            {
                Info<< ": no residual control data found. "
                    << "Calculations will employ " << nCorrPIMPLE_
                    << " corrector loops" << nl;
            }
            else
            {
                Info<< ": max iterations = " << nCorrPIMPLE_ << nl;

                for (const fieldData& ctrl : residualControl_)
                {
                    Info<< "    field " << ctrl.name << token::TAB
                        << ": relTol " << ctrl.relTol
                        << ", tolerance " << ctrl.absTol
                        << nl;
                }
            }
        }
        else
        {
            Info<< ": Operating solver in PISO mode" << nl;
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pimpleControl::loop()
{
    read();

    ++corr_;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    setFirstIterFlag();

    if (corr_ == nCorrPIMPLE_ + 1)
    {
        if (!residualControl_.empty() && (nCorrPIMPLE_ != 1))
        {
            Info<< algorithmName_ << ": not converged within "
                << nCorrPIMPLE_ << " iterations" << endl;
        }

        corr_ = 0;
        mesh_.data::remove("finalIteration");
        return false;
    }

    bool completed = false;
    if (converged_ || criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                << " iterations" << endl;

            mesh_.data::remove("finalIteration");
            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
            storePrevIterFields();

            mesh_.data::add("finalIteration", true);
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            mesh_.data::add("finalIteration", true);
        }

        if (corr_ <= nCorrPIMPLE_)
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
            storePrevIterFields();
            completed = false;
        }
    }

    return !completed;
}


// ************************************************************************* //
