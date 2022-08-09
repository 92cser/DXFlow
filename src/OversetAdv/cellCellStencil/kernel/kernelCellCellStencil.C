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

#include "kernelCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(kernel, 0);
    addToRunTimeSelectionTable(cellCellStencil, kernel, mesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::cellCellStencils::kernel::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{
    // Number of donors
    label nD = donorCcs.size();

    weights.setSize(nD);

    bool shortC = false;
    forAll(donorCcs, i)
    {
        scalar d = mag(donorCcs[i] - sample);
        if(d < 1e-6)
        {
            shortC = true;
            weights.setSize(nD);
            weights = 0;
            weights[i] = 1.0;
            return;
        }
    }
    if(!shortC)
    {
        // distances
        scalarList d(weights);
        // coefficients
        scalarList R(weights);

        // calculate distance
        forAll(donorCcs, i)
        {
            d[i] = mag(sample - donorCcs[i]);
        }

        // smoothing length
        scalar h = max(d);

        // coefficient
        scalar a = 105/16/3.1415926/h/h/h;

        scalar allWeights = 0;

        forAll(R, i)
        {
            R[i] = d[i]/h;
            if(R[i] <= 1)
            {
                scalar R_cur = R[i];
                weights[i] = a*(1+3*R_cur)*pow3(1-R_cur);
            }
            else
            {
                weights[i] = 0;
            }
            allWeights += weights[i];
        }

        forAll(weights, i)
        {
            weights[i] = weights[i]/allWeights;
        }
        //Info << "weights = " << weights << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::kernel::kernel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    inverseDistance(mesh, dict, false)
{
    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::kernel::~kernel()
{}


// ************************************************************************* //
