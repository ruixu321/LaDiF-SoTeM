
	//define some constant for thermocouple measurement
	dimensionedScalar Ctc("Ctc", dimensionSet(1, 0, -2, -1, 0, 0, 0), 653.25);  //real value 985 //200 for wallfire
	dimensionedScalar htc0("htc0", dimensionSet(1, 0, -3, -1, 0, 0, 0), 100);
	dimensionedScalar d_TC("d_TC", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.001); // diameter of the thermocouple
	dimensionedScalar tInf("tInf", T.dimensions(), 294.75); //303 wall fire //294.75 compartment fire



            volScalarField htc
            (
                IOobject
                (
                    "htc",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
		dimensionedScalar("htc",htc0.dimensions(), 0)
            );

            volScalarField h_Nu
            (
                IOobject
                (
                    "h_Nu",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
		dimensionedScalar("h_Nu",htc0.dimensions(), 0)
            );

	IOobject TCHeader
        (
		"TC",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	);
	volScalarField TC
	(
                TCHeader,
		mesh,
		dimensionedScalar("TC",T.dimensions(), 0)
	);
	if(TCHeader.headerOk())
	{
/*
            volScalarField TC
            (
                TCHeader,
		mesh,
		dimensionedScalar("TC",T.dimensions(), 0)
            );

            volScalarField TC_Temp
            (
                IOobject
                (
                    "TC",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
		mesh
            );
*/
		Info<<"Initialize Thermocouple Temperature with TC..."<<endl;
		//TC = TC_Temp;
	}
	else
	{
/*
		volScalarField TC
		(
			TCHeader,
			T
		);
*/
		TC = T;
		Info<<"NO Initial Thermocouple Temperature..."<<endl;
		Info<<"Initialize Thermocouple Temperature with T..."<<endl;
	}

            volScalarField TCOld
            (
                IOobject
                (
                    "TCOld",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
		TC
            );

            volScalarField dTdt
            (
                IOobject
                (
                    "dTdt",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
		dimensionedScalar("dTdt",T.dimensions()/dimTime, 0)
	    );

	    const radiation::radiationModel& rModel
	    = mesh.lookupObject<radiation::radiationModel>("radiationProperties");
	    Info<< "Creating radiation Model   " <<rModel.name()<< endl;

	    //const fvDOM& rFvDOM = dynamic_cast<const radiation::fvDOM&>(rModel);
