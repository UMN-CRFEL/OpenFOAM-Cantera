{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        Qdot
      + fvOptions(rho, he)
    );
	
	forAll(Y, k)
	{
		EEqn -= fvc::laplacian(turbulence->alphaEff()*hsi[k], Y[k]);
		EEqn -= fvc::div(J[k], hsi[k], "div(Ji,hsi)");
	}

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info << "min/max(T) = "
         << min(T).value() << ", " << max(T).value() << endl;
}

// write Lewis Number out @ maximum temperature cell
/*forAll(mesh.cells(), celli)
{
	if ( abs( T[celli]-max(T).value() ) < 1.0e-10 )
	{
		forAll (Y, k)
		{
            scalar Diff =  dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
                              (composition).D(p, T, Y)[k][celli];
			if (Diff > 1.0e-20) 
			{
				scalar Lewis = thermo.alpha()[celli]/Diff;
	            LewisNumber
                << k << tab << Lewis << endl;	
			}

		}
	}
}*/
