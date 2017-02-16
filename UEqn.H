    // Solve the momentum equation
    vector U_up, U_down;
    vector VolumeForce_previous(0,0,0);
    scalar volumeForceDeltaX = 0;
    scalar volumeForceDeltaY = 0;
    scalar volumeForceDeltaZ = 0;
    scalar volumeForceDeltaCounter = 1;
    
    forAll(mesh.C(), cellI) { 
      if(diskID[cellI] == 2) { 
    //Store previous iteration value
    VolumeForce_previous.x() = VolumeForce[cellI].x();
    VolumeForce_previous.y() = VolumeForce[cellI].y();
    VolumeForce_previous.z() = VolumeForce[cellI].z();
    
    U_up.x() = UarbForward[cellI].x();
    U_up.y() = UarbForward[cellI].y();
    U_up.z() = UarbForward[cellI].z();
    
    U_down.x() = UarbBackward[cellI].x();
    U_down.y() = UarbBackward[cellI].y();
    U_down.z() = UarbBackward[cellI].z();
    
    scalar volume = mesh.V()[cellI];
    
    VolumeForce[cellI] = AcDisk.calcVolumeForce(mesh.C()[cellI], DiskRad[cellI], U_up, U_down, DiskTheta[cellI], DiskGama[cellI], DiskChord[cellI], relativeVelocity[cellI], beta[cellI], alpha[cellI], Cl[cellI], Cd[cellI], Torque[cellI], volume, AverageVelocity[cellI]);
        
    //Explicitly relax source terms
    VolumeForce[cellI].x() = VolumeForce_previous.x() + relaxFactor*(VolumeForce[cellI].x() - VolumeForce_previous.x());
    VolumeForce[cellI].y() = VolumeForce_previous.y() + relaxFactor*(VolumeForce[cellI].y() - VolumeForce_previous.y());
    VolumeForce[cellI].z() = VolumeForce_previous.z() + relaxFactor*(VolumeForce[cellI].z() - VolumeForce_previous.z());
    
    volumeForceDeltaX = VolumeForce[cellI].x() - VolumeForce_previous.x() + volumeForceDeltaX;
    volumeForceDeltaY = VolumeForce[cellI].y() - VolumeForce_previous.y() + volumeForceDeltaY;
    volumeForceDeltaZ = VolumeForce[cellI].z() - VolumeForce_previous.z() + volumeForceDeltaZ;
    volumeForceDeltaCounter = volumeForceDeltaCounter + 1;
        
    
    VolumeForce_previous.x() = 0;
    VolumeForce_previous.y() = 0;
    VolumeForce_previous.z() = 0;
      }
    }
    
        if(Pstream::master() == true) { 
      Sout << "\nVolumeForce Calculated " << nl << endl;
    }
    
    AcDisk.Post(mesh, Torque, VolumeForce, diskID, p, U, phi, pInlet, DiskTheta, mesh.V(), DiskRad);
    
    tmp<fvVectorMatrix> HUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    // Get under-relaxation factor
    const scalar UUrf = mesh.solutionDict().equationRelaxationFactor(U.name());

    // Momentum solution
    solve
    (
        relax(HUEqn(), UUrf)
     ==
        -fvc::grad(p)
    );
