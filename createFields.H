Info<< "Creating combustion model\n" << endl;

autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New
    (
        mesh
    )
);

Info<< "Reading thermophysical properties\n" << endl;

psiReactionThermo& thermo = combustion->thermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}


volScalarField hth_ad // kvm
(
    IOobject
    (
        "hth_ad",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField hth_m // kvm
(
    IOobject
    (
        "hth_m",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField h_thermal // kvm
(
    IOobject
    (
        "h_thermal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.he()
);


#include "readAdditionalThermo.H" // kvm

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

#include "createMRF.H"


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

thermo.correct(); // kvm
rho = thermo.rho(); // kvm

mesh.setFluxRequired(p_rgh.name());

#include "phrghEqn.H"


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

IOdictionary additionalControlsDict
(
    IOobject
    (
        "additionalControls",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Switch solvePrimaryRegion
(
    additionalControlsDict.lookup("solvePrimaryRegion")
);

// extinction model

//Switch extinctionModel
//(
//    additionalControlsDict.lookupOrDefault<bool>("extinctionModel", false)
//);




Switch solvePyrolysisRegion
(
    additionalControlsDict.lookupOrDefault<bool>("solvePyrolysisRegion", true)
);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime/dimVolume, 0.0) // kvm
);

volScalarField dQ_EDC
(
    IOobject
    (
        "dQ_EDC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ_EDC", dimEnergy/dimTime/dimVolume, 0.0) // kvm
);

volScalarField omegaFuel
(
    IOobject
    (
        "omegaFuel",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaFuel", dimMass/dimTime/dimVolume, 0.0) // kvm
);

//extinction model
// volScalarField dQ_reIgn
// (
//     IOobject
//     (
//         "dQ_reIgn",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("dQ_reIgn", dimEnergy/dimTime/dimVolume, 0.0) // kvm
// );



Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


volScalarField kappa // kvm
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()*thermo.alpha()
);

volScalarField cp // kvm
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()
);

singleStepReactingMixture<gasHThermoPhysics>& singleMixture // kvm
(
    dynamic_cast<singleStepReactingMixture<gasHThermoPhysics>&>
    (thermo)
);

// stoichiometric O2 to fuel ratio
scalar s(singleMixture.s().value()); // kvm

// stoichiometric air to fuel ratio
scalar stoicRatio(singleMixture.stoicRatio().value());   // kvm

// heat of combustion [J/kg]
scalar qFuel(singleMixture.qFuel().value());   // kvm
Info << "heat of combustion = " << qFuel << nl;

label fuelIndex(singleMixture.fuelIndex()); // kvm
const label inertIndex(composition.species()[inertSpecie]);
// label inertIndex(singleMixture.inertIndex()); // kvm

const volScalarField* O2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
if (thermo.composition().contains("O2"))  // kvm
{
    O2Ptr = &thermo.composition().Y("O2"); // kvm
}
const volScalarField& O2 = *O2Ptr; // kvm

const volScalarField* CO2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
if (thermo.composition().contains("CO2"))  // kvm
{
    CO2Ptr = &thermo.composition().Y("CO2"); // kvm
}
const volScalarField& CO2 = *CO2Ptr; // kvm

const volScalarField& fu = Y[fuelIndex];  // kvm

scalar YO2Inf = 0.23301; //hardcode for now

// Fuel molecular weight
scalar Wu_ = singleMixture.speciesData()[fuelIndex].W(); // kvm
// Fuel enthalpy of formation
scalar Hu_ = singleMixture.speciesData()[fuelIndex].hc(); // kvm

// compute stoichiometric mixture fraction
scalar ftSt = 1.0 / ( 1.0 + stoicRatio ); // kvm
Info << "stoichiometric mixture fraction is = " << ftSt << nl; // kvm


// create fileds for mixture fraction
volScalarField ft // kvm
(
    IOobject
    (
        "ft",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (fu*s-O2+YO2Inf)/(s+YO2Inf)
);


// soot fraction

volScalarField fv(
		IOobject
		(
		 "fv",
		 runTime.timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero",dimless,0.0)
		);

volScalarField soot(
		IOobject
		(
		 "soot",
		 runTime.timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero",dimless,0.0)
		);

#include "createClouds.H"
#include "createSurfaceFilmModel.H"
#include "createPyrolysisModel.H"
#include "createRadiationModel.H"


Info<< "Creating variables for radiation model\n" << endl;

volScalarField dQrad
(
 IOobject
 (
  "dQrad",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("dQrad", dimEnergy/dimTime/dimVolume,0.0)
 );

volScalarField dQradabs
(
 IOobject
 (
  "dQradabs",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("dQradabs", dimEnergy/dimTime/dimVolume,0.0)

 );

volScalarField dQradem
(
 IOobject
 (
  "dQradem",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("dQradem", dimEnergy/dimTime/dimVolume,0.0)
 );

Info<< "Creating Z, Zvar, chi_st, and W_HRR\n" << endl;

volScalarField Z
(
  IOobject
   (
      "Z",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
   ),
   mesh
   // dimensionedScalar("Z", dimless,0.0)
 );

volScalarField Zvar
(
   IOobject
   (
      "Zvar",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
   ),
   mesh
   // dimensionedScalar("Zvar", dimless,0.0)
 );

volScalarField chi_st
(
    IOobject
    (
        "chi_st",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    // dimensionedScalar("chi_st",dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

volScalarField chi_stInstan
(
    IOobject
    (
        "chi_stInstan",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("chi_stInstan",dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

// volScalarField IntermF
// (
//  IOobject
//  (
//   "IntermF",
//   runTime.timeName(),
//   mesh,
//   IOobject::NO_READ,
//   IOobject::AUTO_WRITE
//  ),
//  mesh,
//  dimensionedScalar("IntermF", dimless,0.0)
//  );

volScalarField Z_left
(
 IOobject
 (
  "Z_left",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("Z_left", dimless,0.0)
 );

volScalarField Z_right
(
 IOobject
 (
  "Z_right",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("Z_right", dimless,0.0)
 );

volScalarField W_HRR
(
 IOobject
 (
  "W_HRR",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("W_HRR", dimless,0.0)
 );

Info<< "Reading combustionProperties dictionary\n" << endl;   
IOdictionary combustionProperties
(
    IOobject
    (
        "combustionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar C_chem(readScalar(combustionProperties.lookup("C_chem")));

Info<< "C_chem =  " << C_chem << endl;



// volScalarField ZgradCells
// (
//  IOobject
//  (
//   "ZgradCells",
//   runTime.timeName(),
//   mesh,
//   IOobject::NO_READ,
//   IOobject::AUTO_WRITE
//  ),
//  mesh,
//  dimensionedScalar("ZgradCells", dimensionSet(0,-1,0,0,0,0,0), 0.0)
//  );

// volScalarField deltaCells
// (
//  IOobject
//  (
//   "deltaCells",
//   runTime.timeName(),
//   mesh,
//   IOobject::NO_READ,
//   IOobject::AUTO_WRITE
//  ),
//  mesh,
//  dimensionedScalar("deltaCells", dimensionSet(0,1,0,0,0,0,0), 0.0)
//  );



// #include "createRTI.H"
