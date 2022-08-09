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

#include "PistonWaveGeneration.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PistonWaveGeneration::
PistonWaveGeneration
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    dict_
    (
        IOobject
        (
            "inputDict",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    wave(waterWaveModels::waterWaves::New(dict_)),
    Displacement_(Zero),
    waterLevel_(0.0)
{}

Foam::PistonWaveGeneration::
PistonWaveGeneration
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    dict_
	(
        IOobject
        (
            "inputDict",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    wave(waterWaveModels::waterWaves::New(dict_)),
    Displacement_(Zero),
    waterLevel_(0.0)
{
    Info << "The piston-type wavemaker is adopted for " << wave->name() << endl;
    wave->PrintWaveProperties();

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::PistonWaveGeneration::
PistonWaveGeneration
(
    const PistonWaveGeneration& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    dict_
    (
    IOobject
        (
            "inputDict",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    wave(waterWaveModels::waterWaves::New(dict_)),
    Displacement_(ptf.Displacement_),
    waterLevel_(0.0)
{}


Foam::PistonWaveGeneration::
PistonWaveGeneration
(
    const PistonWaveGeneration& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    dict_
    (
        IOobject
        (
            "inputDict",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    wave(waterWaveModels::waterWaves::New(dict_)),
    Displacement_(ptf.Displacement_),
    waterLevel_(0.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PistonWaveGeneration::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volScalarField& alpha = db().lookupObject<volScalarField>("alpha.water");
    const fvMesh& mesh = alpha.mesh();
    const Time& t = mesh.time();
    scalar deltaT = db().time().deltaTValue();
	
    label patchId = mesh.boundaryMesh().findPatchID(this->patch().name());
	
    const scalarField& boundarySurface_ = mesh.magSf().boundaryField()[patchId];
    const scalarField& boundaryAlpha_ = alpha.boundaryField()[patchId]; 
    const scalar zMin = gMin(this->patch().localPoints().component(2)); 
    const scalar zMax = gMax(this->patch().localPoints().component(2)); 
    const scalar zSpan = zMax-zMin;  

    Cal_waterLevel(boundarySurface_, boundaryAlpha_, zSpan);

    Displacement_.component(0) = wave->displacement(t.value(), deltaT, waterLevel_);

    Info << "Time and final displacement is " << t.value() <<  "  "  <<  Displacement_.component(0) << endl;

    Field<vector>::operator=(Displacement_);

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::PistonWaveGeneration::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("Displacement_")
        << Displacement_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

void Foam::PistonWaveGeneration::Cal_waterLevel
(
   const scalarField& boundarySurface,
   const scalarField& boundaryAlpha,
   scalar zSpan
)
{
   scalar sumOfArea = 0;
   scalar sumOfAreaAlpha = 0;
   
   forAll(boundarySurface , i)
   {
       sumOfArea += boundarySurface[i];                                    
       sumOfAreaAlpha += boundaryAlpha[i]*boundarySurface[i];             
   }
   reduce(sumOfArea, sumOp<scalar>());
   reduce(sumOfAreaAlpha, sumOp<scalar>());
   scalar fraction = sumOfAreaAlpha/sumOfArea;
   waterLevel_ = fraction * zSpan;  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        PistonWaveGeneration
    );
}

// ************************************************************************* //
