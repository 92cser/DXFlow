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

#include "polynomialCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "Matrix.H"
#include "QRMatrix.H"
#include "MatrixBlock.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(polynomial, 0);
    addToRunTimeSelectionTable(cellCellStencil, polynomial, mesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::cellCellStencils::polynomial::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{
    // Number of donors
    label nD = donorCcs.size();
    // set the length of weights
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
        // calculate
        bool cubicPolynomial = true;
        label index = 3;

        // too many donors, use quartic polynomials to yield an over-determined system
        if(nD > 20)
        {
            cubicPolynomial = false;
            index = 4;
        }
        // number of row
        label nRow = cubicPolynomial? 20 : 35;

        // AX = B, X the weights
        scalarRectangularMatrix A(nRow, nD), X(nD, 1), B(nRow, 1);

        // update matrix and source
        updateMatrix(A, donorCcs, index);
        updateSource(B, sample, index);

        // solve
        X = pinv(A)*B;

        // transfer result to scalarList
        weights.setSize(nD);
        forAll(weights, i)
        {
            weights[i] = X[i][0];
        }

        // adjust weights to avoid the problem of undesired volume fraction values
        // for floating strcture simulation
        if(adjustWeights_)
        {
            const scalar tolForAdj = 0.1;
            scalarList weights_backup = weights;
            // block the weight of the main donor
            weights_backup[0] *= 0;
            // find and average
            label minWeight = findMin<scalarList>(weights_backup);
            label maxWeight = findMax<scalarList>(weights_backup);
            if(weights[minWeight] < 0 && weights[maxWeight] > 0)
            {
                scalar moderate = weights[minWeight] + weights[maxWeight];
                if(mag(moderate)< tolForAdj)
                {
                    //Info << "weights[minWeight] = " << weights[minWeight] << endl;
                    //Info << "weights[maxWeight] = " << weights[maxWeight] << endl;
                    weights[minWeight] =  0.5*moderate;
                    weights[maxWeight] =  0.5*moderate;
                }
            }
        }
        //Info << "weights: " << weights << endl;
        //Info << "sum: " << sum(weights) << endl;
    }
}

void Foam::cellCellStencils::polynomial::updateMatrix
(
    scalarRectangularMatrix& matrixInput,
    const pointList& donorCcs,
    const label& index
) const
{
    label nCol = matrixInput.n();
    for(label col = 0; col < nCol; col++) // loop each column
    {
        label currentPosition = 0;
        for(label i = 0; i <= index; i++)
        {
            for(label j = 0; j <= index; j++)
            {
                for(label k = 0; k <= index; k++)
                {
                    if(i+j+k <= index)
                    {
                        scalar x = donorCcs[col].x();
                        scalar y = donorCcs[col].y();
                        scalar z = donorCcs[col].z();
                        matrixInput[currentPosition++][col] = pow(x,i)*pow(y,j)*pow(z,k);
                    }
                }
            }
        }
    }
}

void Foam::cellCellStencils::polynomial::updateSource
(
    scalarRectangularMatrix& matrixInput,
    const point& sample,
    const label& index
) const
{
    label currentPosition = 0;
    for(label i = 0; i <= index; i++)
    {
        for(label j = 0; j <= index; j++)
        {
            for(label k = 0; k <= index; k++)
            {
                if(i+j+k<=index)
                {
                    scalar x = sample.x();
                    scalar y = sample.y();
                    scalar z = sample.z();
                    matrixInput[currentPosition++][0] = pow(x,i)*pow(y,j)*pow(z,k);
                }
            }
        }
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::polynomial::polynomial
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    inverseDistance(mesh, dict, false),
    adjustWeights_(dict.getOrDefault<bool>("adjustWeights", false))
{
    if (doUpdate)
    {
        update();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::polynomial::~polynomial()
{}


// ************************************************************************* //
