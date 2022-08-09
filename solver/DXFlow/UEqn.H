    MRF.correctBoundaryVelocity(U);

    surfaceScalarField muEff
    (
        "muEff",
        mixture.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        WPM.C_item()* fvm::ddt(rho, U)
	  + WPM.Poro_inverse_()*fvm::div(rhoPhi/porositySF, U)
      + WPM.C_item()* MRF.DDt(rho, U)
      - fvm::laplacian(muEff/porositySF, U)
      - WPM.Poro_inverse_()*(fvc::grad(U) & fvc::grad(muEff))
	  + WPM.A_item()* mixture.mu()*U
	  + WPM.B_item()* rho * mag(U) * U
      //+ WPM.Poro_inverse_()*U*rho*WD.DampingField()
     ==
        fvOptions(rho, U)
		- fvm::Sp(WPM.Poro_inverse_()*rho*WD.DampingField(), U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
