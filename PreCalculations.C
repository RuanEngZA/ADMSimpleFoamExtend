/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ActuatorDiskModel_21122015

Description
    Steady-state solver for buoyant, turbulent flow of compressible fluids,
    including radiation, for ventilation and heat-transfer.

\*---------------------------------------------------------------------------*/

#include "ActuatorDisk.H"
#include "fvMesh.H"
#include "unitConversion.H"
#include <stdlib.h>

namespace Foam {
  defineTypeNameAndDebug(PreCalcs, 1);
  
  PreCalcs::PreCalcs(){}	//Default Constructor
  PreCalcs::~PreCalcs(){}	//Default Destructor
  
  void PreCalcs::readDataB2(const fvMesh &iMesh) { 
  
  Istream& is1 = iMesh.solutionDict().subDict("B2Fan").lookup("AcDiskLocation");
  is1.format(IOstream::ASCII); 
  is1 >> GenesisAcDiskLocation;
  
  Istream& is2 = iMesh.solutionDict().subDict("B2Fan").lookup("EPS");
  is2.format(IOstream::ASCII); 
  is2 >> EPS;
  
  Istream& is3 = iMesh.solutionDict().subDict("B2Fan").lookup("HubSize");
  is3.format(IOstream::ASCII); 
  is3 >> HubSize;
  
  Istream& is4 = iMesh.solutionDict().subDict("B2Fan").lookup("AcDiskWidth");
  is4.format(IOstream::ASCII); 
  is4 >> AcDiskWidth;
  
  Istream& is5 = iMesh.solutionDict().subDict("B2Fan").lookup("FlowDirection");
  is5.format(IOstream::ASCII); 
  is5 >> FlowDir;
  
  Istream& is6 = iMesh.solutionDict().subDict("B2Fan").lookup("FanDiameter");
  is6.format(IOstream::ASCII); 
  is6 >> FanDiameter;
  
  Istream& is7 = iMesh.solutionDict().subDict("B2Fan").lookup("ChordLength");
  is7.format(IOstream::ASCII); 
  is7 >> ChordLength;
  
  Istream& is8 = iMesh.solutionDict().subDict("B2Fan").lookup("DiskSize");
  is8.format(IOstream::ASCII); 
  is8 >> DiskSize;
  
  Istream& is9 = iMesh.solutionDict().subDict("B2Fan").lookup("UpDiskLocation");
  is9.format(IOstream::ASCII); 
  is9 >> GenesisUpDiskLocation;
  
  Istream& is10 = iMesh.solutionDict().subDict("B2Fan").lookup("DownDiskLocation");
  is10.format(IOstream::ASCII); 
  is10 >> GenesisDownDiskLocation;
 
  
  Info << "**B2Fan Properties** " << endl;
  Info << "Fan Diameter: " << FanDiameter << endl;
  Info << "Hub Size: " << HubSize << endl;
  Info << "Chord Length: " << ChordLength  << endl;
  Info << "Fan Speed: " << FanSpeed << endl;
  Info << "Numbers of Blade: " << BladeNo << endl;
  Info << "Blade Angle at Hub: " << GamaHub << "\n" << endl;
  
  Info << "**Actuator Disk Properties** " << endl;
  Info << "First Actuator Disk Location: " << GenesisAcDiskLocation << endl;
  Info << "First Upstream Disk Location: " << GenesisUpDiskLocation << endl;
  Info << "First Downstream Disk Location: " << GenesisDownDiskLocation << endl;
  Info << "Actuator Disk Width: " << AcDiskWidth << endl;
  Info << "Flow Direction: " << FlowDir << endl;
  Info << "Tolerance: " << EPS << "\n" << endl;
 
  }
  
  void PreCalcs::setDiskID(const fvMesh &iMesh, scalarField &idiskID, scalarField &iLarge_for, scalarField &iLarge_back) { 
  /*--------------------------------------------------------------------------
			      Set Disk ID
    Loop through all cells and set the volScalarField diskID to:
    
    2 for Actuator Disk
    1 for Upstream Disk
    3 for Downstream Disk
    
    --------------------------------------------------------------------------*/
    
  Info << "Setting diskID for Actuator Disk\n" << endl;
  //Check if the point is situated in the axial direction of the Actuator Disk
  vector AcDiskAxialSearch(0, 0, 0);
  AcDiskAxialSearch.z() = AcDiskWidth/2;
	 
  forAll(iMesh.C(), cellI) {     
    //Check if the cell is situated in the axial and radial direction of the disk
    if(AbsRangeCheckZ(iMesh.C()[cellI], GenesisAcDiskLocation, AcDiskAxialSearch) && radRangeCheckZ(iMesh.C()[cellI], GenesisAcDiskLocation, FanDiameter, HubSize)) { 		 			 	idiskID[cellI] = 2;
    }
	}
	
	//Determine up and downstream disks
	vector UpDiskCentreLoc(GenesisUpDiskLocation);
	vector DownDiskCentreLoc(GenesisDownDiskLocation);	    
	vector DiskAxialSearch (0, 0, 0);
	DiskAxialSearch.z() = DiskSize/2;

	forAll(iMesh.C(), cellI1) { 
	  if (AbsRangeCheckZ(iMesh.C()[cellI1], UpDiskCentreLoc, DiskAxialSearch) && radRangeCheckZ(iMesh.C()[cellI1], UpDiskCentreLoc, FanDiameter, HubSize)) { 
	    idiskID[cellI1] = 1;
	    iLarge_back[cellI1] = 1;
		  }
	      
	      if (AbsRangeCheckZ(iMesh.C()[cellI1], DownDiskCentreLoc, DiskAxialSearch) && radRangeCheckZ(iMesh.C()[cellI1], DownDiskCentreLoc, FanDiameter, HubSize)) { 
		idiskID[cellI1] = 3;
		iLarge_for[cellI1] = 1;
		  }
	      }
  }
  
  void PreCalcs::setDatStop(const fvMesh &iMesh, scalarField &idiskID, scalarField &iDatStop) { 
  /*--------------------------------------------------------------------------
   * 				Set Data Stop Vector
   * Sets the value of DatStop volScalarField to 1 outside of actuator disk 
   * data transport region and 0 within the transport region
   * ------------------------------------------------------------------------*/
  vector DatStopUpLoc(GenesisAcDiskLocation);
  vector DatStopAxial(0,0,0);
  DatStopAxial.z() = -(GenesisDownDiskLocation.z() - GenesisUpDiskLocation.z()) + EPS;
  
  forAll(iMesh.C(), cellI) {     
    if(AbsRangeCheckZ(iMesh.C()[cellI], DatStopUpLoc, DatStopAxial/2) && radRangeCheckZ(iMesh.C()[cellI], DatStopUpLoc, FanDiameter, HubSize)) { 
      iDatStop[cellI] = 0;
      }
    }  
  }
        
  bool PreCalcs::AbsRangeCheckZ(const vector &MeshC, const vector &Genesis, const vector &length) {  
    scalar z_up = Genesis.z() + length.z();
    scalar z_down = Genesis.z() - length.z();
    
    return (MeshC.z() >= z_down && MeshC.z() <= z_up);
  }
  
  bool PreCalcs::radRangeCheckZ(const vector &MeshC, const vector &Genesis, const scalar &fanDiameter, const scalar &hubSize) {   
    vector PointLine(MeshC - Genesis);
    scalar PointX = PointLine.x()*PointLine.x();
    scalar PointY = PointLine.y()*PointLine.y();
    scalar fanRad = (fanDiameter/2)*(fanDiameter/2);
    scalar hubRad = (hubSize/2)*(hubSize/2);
  
    return ((PointX + PointY) <= fanRad && (PointX + PointY) > hubRad) ;
  
  }
  
  
}