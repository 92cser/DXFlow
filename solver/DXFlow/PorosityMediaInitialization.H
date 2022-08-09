volScalarField PorosityMedia                  
(
	IOobject
	(
		"PorosityMedia",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	0.0,
	"zeroGradient"
);

volScalarField PorosityField         
(
	IOobject
	(
		"PorosityField",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	1.0,
	"zeroGradient"
);

// wave porosity media 
// modify the field value at the same time 
wavePorosityMedia WPM(inputDict,PorosityMedia, PorosityField);

surfaceScalarField porositySF = fvc::interpolate(PorosityField);



