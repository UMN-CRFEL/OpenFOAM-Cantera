tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);


PtrList<volVectorField> Dm =
dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
(composition).Dm(p, T, Y, rho);

// - Diffusivity (m2/s)
PtrList<volScalarField> D =
dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
(composition).D(p, T, Y, rho);

PtrList<volScalarField> hsi(Y.size());
PtrList<surfaceScalarField> J(Y.size());
PtrList<volScalarField> D(Y.size());
forAll(Y, i)
{
    hsi.set
	(
        i,        
        new volScalarField
        (
            IOobject
            (
        	    "hsi",
                mesh.time().timeName(),
        	    mesh
            ),
            mesh,
            dimensionedScalar
            (
                "hsi", 
        	    dimEnergy/dimMass, 
                Zero
            )
        )
	);

    J.set
	(
        i, 		
        linearInterpolate(Dm[i]*rho) & mesh.Sf()
	);	
}

forAll (Y, i)
{
	volScalarField& tHsi = hsi[i];	
	forAll(tHsi, celli)
	{
		tHsi[celli] = composition.Hs(i, p[celli], T[celli]);
	}
	volScalarField::Boundary& Bf = tHsi.boundaryFieldRef();	
	forAll(Bf, patchi)
	{
		forAll(Bf[patchi], facei)	
		{
			Bf[patchi][facei] = 
			composition.Hs
			(
			    i,
				p.boundaryField()[patchi][facei],
				T.boundaryField()[patchi][facei]
			);
		}
	}
}


{
    volScalarField alphat = turbulence->alphat();
	const scalar Let(thermo.lookupOrDefault<scalar>("Let", 1.0));
	volScalarField Dt = alphat/(rho*Let);


forAll(Y, i)
{
    D[i] = D[i] + Dt;
}
	
	reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
			  - fvm::laplacian(rho*D[i], Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
