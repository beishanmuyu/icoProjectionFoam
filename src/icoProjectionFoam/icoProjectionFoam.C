/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    pisoControl piso(mesh);

#include "createFields.H"
#include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n"
         << endl;

    dimensionedScalar dt = dimensionedScalar("dt", dimless * dimTime, 1.0);

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

#include "CourantNo.H"

        // Momentum predictor
        fvVectorMatrix UEqn(fvm::laplacian(nu, U));
        solve(UEqn);

        // Pressure corrector
        while (piso.correctNonOrthogonal())
        {
            phi = Foam::fvc::interpolate(U) & mesh.Sf();
            fvScalarMatrix pEqn(
                -fvm::laplacian(p) == (1 / dt) * fvc::div(phi) //Solving Pressure Poisson eqs
            );
            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
            U = U + dt * fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
