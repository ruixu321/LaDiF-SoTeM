
	const volScalarField& Gcst = mesh.lookupObject<volScalarField>("G");
       //solve the thermocouple equation
       Info<< "Creating field TC\n" << endl;
       //ht=h0*(T/tInf)^0.8, but use ^0.75 for now.
       htc = htc0*(T/tInf)/sqrt(sqrt(T/tInf));

	//Solve the convective heat transfer coefficient for any kind of fire
       
       volScalarField mu1 = thermo.mu();
       volScalarField kappa1 = thermo.kappa();
       
       
       volScalarField Re_ = rho*mag(U)*d_TC/mu1;
       scalar Pr_ = 0.7;
       volScalarField Nu_ = 2.00+0.6*pow(Re_,1/2)*pow(Pr_,1/3);
       h_Nu = Nu_*kappa1/d_TC;

      //Calculate the thermocouple temperature
       Info << "Calculate the thermocouple temperature" << endl;
       TC = TCOld + (1.00*Gcst/4.0 + h_Nu*(T-TCOld)-1.00*physicoChemical::sigma*pow4(TCOld))*runTime.time().deltaT()/Ctc;
       //   dTdt = (TC - TCOld)/runTime.time().deltaT();

       Info << "TC calculated" << endl;

       TCOld = TC;
