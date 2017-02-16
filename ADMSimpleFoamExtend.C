/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow
    Consistent formulation without time-step and relaxation dependence by Jasak

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "ActuatorDisk.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createDiskFields.H"

    simpleControl simple(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    //========================================================================
    //              Preliminary Field Loops
    //========================================================================
    
    PreCalcs Pre;
    Pre.readDataB2(mesh);
    Pre.setDiskID(mesh, diskID, Large_for, Large_back);
    Pre.setDatStop(mesh, diskID, DatStop);
    
    ActuatorDiskModel AcDisk;
    AcDisk.readDataB2(mesh);
    
    forAll(mesh.C(), cellI) {
      if(diskID[cellI] == 2) {
    DiskRad[cellI] = AcDisk.getRadius(mesh.C()[cellI]);
    DiskTheta[cellI] = AcDisk.getTheta(mesh.C()[cellI], DiskRad[cellI]);
    DiskGama[cellI] = AcDisk.getBladeAngle(DiskRad[cellI]);
    DiskChord[cellI] = AcDisk.getChordLength(DiskRad[cellI]);
      }
    }
    
    //Read in Airfoil profile data for NASA 0413-LS profile
    AcDisk.readCSV_NASA0413LS();
    
    Info << "\n!!Running ActuatorDiskModel_15022016!!\n" << endl;

    Info<< "\nStarting time loop\n" << endl;
    
    scalar relaxFactor;
    
    Istream& is1 = mesh.solutionDict().subDict("B2Fan").lookup("VolForceRelax");
    is1.format(IOstream::ASCII); 
    is1 >> relaxFactor;

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
#           include "TransportEqn.H"
#           include "UEqn.H"
#           include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
