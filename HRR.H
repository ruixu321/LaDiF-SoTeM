        volScalarField::Boundary& O2Bf = O2.boundaryFieldRef();

        scalar TotalsumFlux = 0.0;

        forAll(O2Bf, patchi){

            word patchlabelName = mesh.boundaryMesh().names()[patchi];


            const scalarField & phiw = phi.boundaryField()[patchi];
            const scalarField & O2w = O2.boundaryField()[patchi];

            // Initialise the sumFlux variable

            scalar sumFlux = 0.0;

            // O2 mass flow rate

            const scalarField wO2 = phiw*O2w;

            // Sum all of the face values and assign it to sumFlux
            sumFlux = Foam::sum(wO2);

            Info << "Mass flow rate on patch " << patchlabelName << " : " << sumFlux <<  " kg/s" << endl;

            TotalsumFlux += sumFlux;
        }

        Info << "Mass flow rate O2 : " << TotalsumFlux <<  " kg/s " << endl;

        scalar HRR = 12.5E6 * TotalsumFlux/1000.;


