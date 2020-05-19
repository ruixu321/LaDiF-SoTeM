// Info << "Extinction Model" << endl;


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


{
    volScalarField dEff
    (
        "dEff",
        turbulence->alphaEff()
    );

    if (!constD)
    {
        dEff -= turbulence->alpha()*(1-1.0/lewisNo);
    }
    else
    {
        dEff = dEff - turbulence->alpha() + DM*rho;
    }

    // radiation->correct();

    Info << " Correct combustion " << endl;

    combustion->correct();

    #include "ZEqn.H" // solve for combustion source term

    Info<< "min/max(dQ) = "
        << min(dQ).value() << ", " << max(dQ).value() << endl;

    

    Info<< "min/max(omegaFuel) = "
        << min(omegaFuel).value() << ", " << max(omegaFuel).value() << endl;

    // dQ is only for diagnostic purposes.
    // TODO: set dQ = combustoin->dQ(); once OpenCFD updates
    //       greyAbsorptionEmission model in OpenFOAM source
    dQ_EDC = combustion->Sh(); // kvm

    Info<< "min/max(dQ_EDC) = "
        << min(dQ_EDC).value() << ", " << max(dQ_EDC).value() << endl;
        
    // dQ_reIgn = combustion->Sh_reIgn(); // kvm

    
    // Info << " dQ_reIgn =  " << dQ_reIgn << endl;

    volScalarField Yt(0.0*Y[0]);

    Info << " Calculating species mass fraction Yi " << endl;

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            const label specieI = thermo.composition().species()[Yi.member()];

            // Info << Yi.member() << "  specieStoichCoeffs = " << singleMixture.specieStoichCoeffs()[specieI] << endl;

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              // FIXME: species uses dEff and energy uses alphaEff
              - fvm::laplacian(dEff, Yi)
              ==
                parcels.SYi(i, Yi)
              + surfaceFilm.Srho(i)
              // + combustion->R(Yi)
              + omegaFuel * singleMixture.specieStoichCoeffs()[specieI] // e.g. for CH4, specieStoichCoeffs()[specieI] will be -1, for O2 will be -4
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

    // Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex] == scalar(1) - Yt; // kvm
    Y[inertIndex].max(0.0);

    forAll(Y, i){
        Info << setw(8) << Y[i].name() << "\tmin/ave/max "<< " =\t"
             << setw(12) << setprecision(5) << gMin(Y[i]) << "\t"
             << setw(12) << setprecision(5) << gAverage(Y[i]) << "\t"
             << setw(12) << setprecision(5) << gMax(Y[i]) << nl;
    }

    Info << " Correct radiation " << endl;
    radiation->correct();

    const volScalarField::Internal T4 = pow(thermo.T(),4).ref();
    const volScalarField::Internal rad = radiation->Ru() - radiation->Rp().ref()*T4;
    volScalarField::Internal& radField = dQrad.ref();
    radField = rad;

    const volScalarField::Internal radRu = radiation->Ru();
    volScalarField::Internal& radFieldRu = dQradabs.ref();
    radFieldRu = radRu;

    const volScalarField::Internal radRp = radiation->Rp().ref()*T4;
    volScalarField::Internal& radFieldRp = dQradem.ref();
    radFieldRp = radRp;



    volScalarField& he = thermo.he();

    Info << " Calculating thermal enthalpy " << endl;

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
        // combustion->Sh()
        dQ
      + radiation->Sh(thermo)
      + parcels.Sh(he)
      + surfaceFilm.Sh()
      + fvOptions(rho, he)
    );


    fvScalarMatrix EadEqn
    (
        fvm::ddt(rho, hth_ad) + mvConvection->fvmDiv(phi, hth_ad)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      -dpdt
      - fvm::laplacian(turbulence->alphaEff(), hth_ad)
     ==
        // combustion->Sh()
        dQ
      //+ radiation->Sh(thermo)
      //+ parcels.Sh(he)
      //+ surfaceFilm.Sh()
      //+ fvOptions(rho, he)
    );

    fvScalarMatrix EmEqn
    (
        fvm::ddt(rho, hth_m) + mvConvection->fvmDiv(phi, hth_m)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      - fvm::laplacian(turbulence->alphaEff(), hth_m)
     ==
      dpdt
    );


    EEqn.relax();
    EadEqn.relax();
    EmEqn.relax();

    fvOptions.constrain(EEqn);

    Info << " EEqn " << endl;
    EEqn.solve();
    Info << " EadEqn " << endl;
    EadEqn.solve();
    Info << " EmEqn" << endl;
    EmEqn.solve();

    fvOptions.correct(he);

    Info << " thermo correct " << endl;
    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;

    forAll(T, celli)
    {
      if (T[celli] > 4000)
      {
        Info << " T[celli] = " << T[celli] << endl;
        Info << " omegaFuel[celli] = " << omegaFuel[celli] << endl;
        Info << " qFuel = " << qFuel << endl;
        Info << " rho[celli] * fu[celli] / deltaT() = " << rho[celli] * fu[celli] / runTime.deltaT().value() << endl;
        Info << " rho[celli] * O2[celli] / 4.0 / deltaT() = " << rho[celli] * O2[celli] / 4.0 / runTime.deltaT().value() << endl;
        Info << " chi_st[celli] = " << chi_st[celli] << endl;
        Info << " pow(mesh.V()[celli], 1/3.) = " << pow(mesh.V()[celli], 1/3.) << endl;
        Info << " W_HRR[celli] = " << W_HRR[celli] << endl;
      }
    }
}


