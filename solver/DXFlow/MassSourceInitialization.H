// Mark the mass source region 
volScalarField MassSource                 
(
	IOobject
	(
		"MassSource",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	0.0,
	"zeroGradient"
);

// Wave source field
volScalarField WaveSourceField                 
(
    IOobject
    (
         "WaveSourceField",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("WaveSourceField",dimensionSet(0,0,-1,0,0,0,0),0.0),
	"zeroGradient"
);
	
massSourceWaveMaker MS
(
	inputDict,
	MassSource,
	WaveSourceField,
	U
);


