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
    Foam::cellCellStencils::kernel

Description
    Kernel intepolation method.

Usage
    in fvSchemes
    oversetInterpolation
    {
        method           kernel;
    }

SourceFiles
    kernelCellCellStencil.C

Author
    Wang Dongxu, Ocean University of China. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef kernelCellCellStencils_H
#define kernelCellCellStencils_H

#include "inverseDistanceCellCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                        Class leastSquares Declaration
\*---------------------------------------------------------------------------*/

class kernel
:
    public inverseDistance
{    
    // Private Member Functions

        //- No copy construct
        kernel(const kernel&) = delete;

        //- No copy assignment
        void operator=(const kernel&) = delete;


public:

    //- Runtime type information
    TypeName("kernel");


    // Constructors

        //- Construct from fvMesh
        kernel(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~kernel();


    // Member Functions

        //- Calculate kernel weights for single acceptor
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
