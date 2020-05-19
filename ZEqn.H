#include <math.h>
#include <algorithm>

//-  Solving equation of mixture fraction
// {
//     volScalarField dEff
//     (
//         "dEff",
//         turbulence->alphaEff()
//     );

//     if (!constD)
//     {
//         dEff -= turbulence->alpha()*(1-1.0/lewisNo);
//     }
//     else
//     {
//         dEff = dEff - turbulence->alpha() + DM*rho;
//     }

//     fvScalarMatrix ZEqn
//     (
//         fvm::ddt(rho, Z)
//       + fvm::div(phi, Z)
//       - fvm::laplacian(dEff, Z)
//     );
//     ZEqn.relax();
//     ZEqn.solve(mesh.solver("Z"));
//     Z.correctBoundaryConditions();
// }

{
    //- Solving equation of variance of mixture fraction
    fvScalarMatrix ZvarEqn
    (
        (
        fvm::ddt(rho, Zvar)
          + fvm::div(phi, Zvar)
          - fvm::laplacian(turbulence->mut()/0.5 + thermo.alpha(), Zvar)
          - 2.0*(turbulence->mut()/0.5)*magSqr(fvc::grad(Z))
          + 2.0*rho*(turbulence->epsilon()/turbulence->k())*Zvar
        )
    );

    ZvarEqn.relax();

    ZvarEqn.solve(mesh.solver("Zvar"));
    
    Zvar.correctBoundaryConditions();
    
    bound(Zvar, 0.0);
}

/* if (varianceEquation == false)
{
    //- Solving equation of variance of mixture fraction
    fvScalarMatrix ZvarEqn
    (
        (
	    fvm::ddt(rho, Zvar)
          + fvm::div(phi, Zvar)
          - fvm::laplacian(turbulence->mut()/sigmat, Zvar)
          - Cg*turbulence->mut()*magSqr(fvc::grad(Z))
          + Cd*rho*turbulence->epsilon()/(turbulence->k()+k_small)*Zvar
        )
    );

    ZvarEqn.relax();

    ZvarEqn.solve(mesh.solver("Zvar"));
}
else
{
    Info << "Algebraic equation for Zvar" << endl;

    scalarField& ZvarCells = Zvar.ref();

    ZvarCells =
        Cg*turbulence->mut()*magSqr(fvc::grad(Z))*turbulence->k()
      / (Cd*rho*(turbulence->epsilon()+epsilon_small));
}
 */


//- Updating variance of mixture fraction
// scalarField& ZvarCells = Zvar.ref();
// ZvarCells = C_Z*magSqr(sgs.delta())*magSqr(fvc::grad(Z));


//- Updating instantaneous and effective stoichiometric dissipation rate

scalarField& chi_stCells = chi_stInstan.ref();    // instantaneous chi_st
chi_stCells=2.0*turbulence->muEff()/(rho)*magSqr(fvc::grad(Z));

forAll(chi_st, celli){
    if (mesh.C()[celli].component(vector::Z) < 0.003){
        chi_st[celli] = chi_stCells[celli];
    }
 }

{
    //- Solving transport equation of effective chi_st
    fvScalarMatrix chi_stEqn
    (
        (
        fvm::ddt(rho, chi_st)
          + fvm::div(phi, chi_st)
          - fvm::laplacian(turbulence->muEff(), chi_st)
          -rho*2.0/0.024781*chi_st*(2.0*turbulence->muEff()/(rho)*magSqr(fvc::grad(Z))-chi_st)
        )
    );
    chi_stEqn.relax();
    chi_stEqn.solve(mesh.solver("chi_st"));
    chi_st.correctBoundaryConditions();
}


Info<< "Calculating Z " << endl;

forAll(Z, celli)
{
    Z[celli] = (fu[celli] - O2[celli] / 4.0 + 0.2314 / 4.0) / (1.0 + 0.2314 / 4.0);
}


Info<< "Calculating Z_left and Z_right" << endl;

forAll(Z_left, celli) 
{
    Z[celli] = max(Z[celli], 0.0);

    Z[celli] = min(Z[celli], 1.0);

    Zvar[celli] = max(Zvar[celli], 0.0);

    Zvar[celli] = min(Zvar[celli], Z[celli] * (1 - Z[celli]));

    if (Zvar[celli] == 0.0)
    {
        Z_left[celli] = Z[celli];
        Z_right[celli] = Z[celli];
    }
    else if (Zvar[celli] > 0 && Zvar[celli] <= pow(Z[celli], 2.0) / 3.0 &&  Zvar[celli] <= pow(1-Z[celli], 2.0) / 3.0)
    {
        Z_left[celli] = Z[celli] - sqrt(3.0 * Zvar[celli]);
        Z_right[celli] = Z[celli] + sqrt(3. * Zvar[celli]);
    }
    else if (Zvar[celli] > pow(Z[celli], 2.) / 3. && Zvar[celli] <= 2./3. * Z[celli] - pow(Z[celli], 2.))
    {
        Z_left[celli] = 0.0;
        Z_right[celli] = 3. * (Zvar[celli] + pow(Z[celli], 2.)) / (2. * Z[celli]);
    }
    else if (Zvar[celli] > pow(1.-Z[celli], 2.) / 3. && Zvar[celli] <= 2./3.*(1.-Z[celli]) - pow(1.-Z[celli], 2.))
    {
        Z_left[celli] = 1. - 3.*(Zvar[celli]+pow(1.-Z[celli], 2.)) / (2.*(1.-Z[celli]));
        Z_right[celli] = 1.0;
    }
    else if (Zvar[celli] > 2./3.*Z[celli] - pow(Z[celli], 2.) && Zvar[celli] > 2./3.*(1.-Z[celli]) - pow(1.-Z[celli], 2.) 
                && Zvar[celli] <= Z[celli] * (1. - Z[celli]))
    {
        Z_left[celli] = 0.0;
        Z_right[celli] = 1.0;
    }
    else 
    {
        Info<< "Error in calculating Z_left and Z_right" << endl;
    }
}


Info<< "Calculating W_HRR" << endl;

forAll(W_HRR, celli) 
{
    W_HRR[celli] = 0.0;
    if (Z[celli] < 1E-3 || 1 - Z[celli] < 1E-3){
        W_HRR[celli] = 0.0;
    }
    else if (0.5 * ftSt <= Z_right[celli] && Z_right[celli] <= 2.*ftSt){
        if (Z_left[celli] < 0.5 * ftSt){
            W_HRR[celli] = (Z_right[celli] - 0.5 * ftSt) / (2. * ftSt - 0.5 * ftSt);
        }
        else{
            W_HRR[celli] = (Z_right[celli] - Z_left[celli]) / (2. * ftSt - 0.5 * ftSt);
        }
    }
    else if (2.*ftSt < Z_right[celli] && Z_left[celli] <= 2.*ftSt){
        if (Z_left[celli] < 0.5 * ftSt){
            W_HRR[celli] = 1.;
        }
        else {
            W_HRR[celli] = (2.*ftSt - Z_left[celli]) / (2.*ftSt - 0.5 * ftSt);
        }
    }
}

Info<< "Calculating dQ and omegaFuel " << endl;

forAll(dQ, celli) 
{
    if (O2[celli] > 0 && fu[celli] > 0)
    {

        chi_st[celli] = max(chi_st[celli], 0.0);

        // dQ[celli] = 0.004 * qFuel * sqrt(chi_st[celli]) / pow(mesh.V()[celli], 1/3.) * W_HRR[celli] * C_chem;  // positive
        omegaFuel[celli] = 0.004 * sqrt(chi_st[celli]) / pow(mesh.V()[celli], 1/3.) * W_HRR[celli];  // positive
        omegaFuel[celli] = min(omegaFuel[celli], rho[celli] * fu[celli] / runTime.deltaT().value());
        omegaFuel[celli] = min(omegaFuel[celli], rho[celli] * O2[celli] / 4.0 / runTime.deltaT().value());
        dQ[celli] = omegaFuel[celli] * qFuel;  // positive

    }
    else 
    {
        dQ[celli] = 0.0;
        omegaFuel[celli] = 0.0;
    }
    
}
