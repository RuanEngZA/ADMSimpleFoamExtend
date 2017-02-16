Info << "Reading Fields Related to Actuator Disk" << endl;

	volScalarField BoundaryDummy
	(
	  IOobject
	  (
	    "BoundaryDummy",
	  runTime.timeName(),
	  mesh,
	  IOobject::MUST_READ,
	  IOobject::AUTO_WRITE
	  ),
	mesh
	);

    Info << "Reading diskID" << endl;
    volScalarField diskID
	(
	  IOobject
	  (
	    "diskID",
	  runTime.timeName(),
	  mesh,
	  IOobject::MUST_READ,
	  IOobject::AUTO_WRITE
	  ),
	mesh
	);
    
    Info << "Reading field UarbForward\n" << endl;
    volVectorField UarbForward
    (
        IOobject
        (
            "UarbForward",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UarbForward", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.0,0.0,1.0)), 
	BoundaryDummy.boundaryField().types()
    );
    
    volVectorField UarbForward_temp
    (
        IOobject
        (
            "UarbForward_temp",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UarbForward_temp", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.0,0.0,1.0)), 
	BoundaryDummy.boundaryField().types()
    );
    
    Info << "Reading field UarbBackward\n" << endl;
    volVectorField UarbBackward
    (
        IOobject
        (
            "UarbBackward",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UarbBackward", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.0,0.0,-1.0)), 
	BoundaryDummy.boundaryField().types()
    );
    
    Info << "Reading field Large\n" << endl;
    volScalarField Large_for
    (
      IOobject
      (
	"Large_for",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("Large_for", dimensionSet(0,0,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField Large_back
    (
      IOobject
      (
	"Large_back",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("Large_back", dimensionSet(0,0,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volVectorField UTransport_1
    (
      IOobject
      (
	"UTransport_1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("UTransport_1", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
    );
    
    volVectorField UTransport_3
    (
      IOobject
      (
	"UTransport_3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("UTransport_3", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
    );
    
    Info << "Reading field DatStop\n" << endl;
    volScalarField DatStop
    (
      IOobject
      (
	"DatStop",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("DatStop", dimensionSet(0,0,-1,0,0,0,0), 1.0), 
      BoundaryDummy.boundaryField().types()
    );
          
    dimensionedScalar DS
    (
      "DS", 
      dimensionSet(0, 2, -1, 0, 0, 0, 0), 
      1e-5
    ); 
    
    Info << "Reading surfaceScalarField phiFor\n"<<endl;
    surfaceScalarField phiFor
    (
	IOobject
	(
	    "phiFor",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	linearInterpolate(UarbForward) & mesh.Sf()
    );
    
    Info << "Reading surfaceScalarField phiBack\n"<<endl;
    surfaceScalarField phiBack
    (
	IOobject
	(
	    "phiBack",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	linearInterpolate(UarbBackward) & mesh.Sf()
    );
    
    volVectorField DiskRad 
    (
      IOobject
      (
	"DiskRad",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("DiskRad", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField DiskTheta 
    (
      IOobject
      (
	"DiskTheta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("DiskTheta", dimensionSet(0,0,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField DiskGama 
    (
      IOobject
      (
	"DiskGama",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("DiskGama", dimensionSet(0,0,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField DiskChord 
    (
      IOobject
      (
	"DiskChord",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("DiskChord", dimensionSet(0,0,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volVectorField VolumeForce 
    (
      IOobject
      (
	"VolumeForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("VolumeForce", dimensionSet(0,1,-2,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
    );
    
    volVectorField VolumeForceRelax 
    (
      IOobject
      (
	"VolumeForceRelax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("VolumeForceRelax", dimensionSet(0,1,-2,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
    );
    
    volVectorField relativeVelocity
    (
      IOobject
      (
	"relativeVelocity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("relativeVelocity", dimensionSet(0,1,-2,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField beta
    (
      IOobject
      (
	"beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("beta", dimensionSet(0,1,-2,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField alpha
    (
      IOobject
      (
	"alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("alpha", dimensionSet(0,1,-2,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField Cl
    (
      IOobject
      (
	"Cl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("Cl", dimensionSet(0,1,-2,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField Cd
    (
      IOobject
      (
	"Cd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("Cd", dimensionSet(0,1,-2,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField cellVolume
    (
      IOobject
      (
	"cellVolume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("cellVolume", dimensionSet(0,1,-2,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField Torque
    (
      IOobject
      (
	"Torque",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("Torque", dimensionSet(0,1,-2,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField pInlet
    (
      IOobject
      (
	"pInlet",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("pInlet", dimensionSet(0,1,-2,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volVectorField AverageVelocity
    (
      IOobject
      (
	"AverageVelocity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector("AverageVelocity", dimensionSet(0,1,-2,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
    );