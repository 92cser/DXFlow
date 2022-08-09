zeroField Su;
zeroField Sp;
zeroField divU;

// current alpha source term
volScalarField alphaCurr = alphaRef - alpha1;

// Source for wave damping and mass source wavemaker
volScalarField Su_Waves                 
(
	IOobject
	(
		 "Su_Waves",
		 runTime.timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::NO_WRITE
	),
	WaveSourceField*alpha1 + WD.DampingField()*alphaCurr*WD.AlphaSource()
);

