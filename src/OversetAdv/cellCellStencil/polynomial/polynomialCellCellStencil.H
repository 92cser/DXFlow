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
    Foam::cellCellStencils::polynomial

Description
    Polynomial intepolation method. The method is not stable in floating
    structure simulation. Therefore, an adjustment procedure is provided.
    The method is sufficiently stable for single-phase simulations or
    two-phase simulations of submerged structures.

    Reference:
    Chandar, D. D. (2019). On overset interpolation strategies and conservation
    on unstructured grids in openfoam. Computer Physics Communications, 239, 72-83.

Usage
    in fvSchemes
    oversetInterpolation
    {
        method           polynomial;
        adjustWeights    true; // default: false
    }

SourceFiles
    polynomialCellCellStencil.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef polynomialCellCellStencils_H
#define polynomialCellCellStencils_H

#include "inverseDistanceCellCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                        Class leastSquares Declaration
\*---------------------------------------------------------------------------*/

class polynomial
:
    public inverseDistance
{    
    // Private Member Functions

        //- No copy construct
        polynomial(const polynomial&) = delete;

        //- No copy assignment
        void operator=(const polynomial&) = delete;

        // update coefficient matrix
        void updateMatrix
        (
            scalarRectangularMatrix& matrixInput,
            const pointList& donorCcs,
            const label& index
        ) const;

        // update source
        void updateSource
        (
            scalarRectangularMatrix& matrixInput,
            const point& sample,
            const label& index
        ) const;

        // adjust weights or not
        // the adjustment is required for floating structure simulation
        bool adjustWeights_;

public:

    //- Runtime type information
    TypeName("polynomial");

    // Constructors

        //- Construct from fvMesh
        polynomial(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~polynomial();


    // Member Functions

        //- Calculate polynomial weights for single acceptor
        virtual void stencilWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCellStencils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
