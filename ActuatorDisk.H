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
    ActuatorDiskModel_

Description
    Actuator Disk Model implemented on a simplified mesh

\*---------------------------------------------------------------------------*/

#ifndef ACDISKS
#define ACDISKS

#include "fvCFD.H"
#include "polyMesh.H"
#include <vector>

namespace Foam {
  
  class PreCalcs {
  public:
    
    TypeName("PreCalcs");
    PreCalcs();
    virtual ~PreCalcs();
    
    void readDataB2(const fvMesh &iMesh);
    void setDiskID(const fvMesh &iMesh, scalarField &idiskID, scalarField &iLarge_for, scalarField &iLarge_back);
    bool AbsRangeCheckZ(const vector &MeshC, const vector &Genesis, const vector &length);
    bool radRangeCheckZ(const vector &MeshC, const vector &Genesis, const scalar &fanDiameter, const scalar &hubSize);
    void setDatStop(const fvMesh &iMesh, scalarField &idiskID, scalarField &iDatStop);
      
  private:
    vector GenesisAcDiskLocation, FlowDir, GenesisUpDiskLocation, GenesisDownDiskLocation;
    scalar EPS, HubSize,  AcDiskWidth, FanDiameter, ChordLength, DiskSize, FanSpeed, BladeNo, GamaHub;
    
  };
  
  class ActuatorDiskModel {
  public: 
    
    TypeName("ActuatorDiskModel");
    ActuatorDiskModel();
    virtual ~ActuatorDiskModel();
    
    void readDataB2(const fvMesh &iMesh);
    vector getRadius(const vector &iMeshC);
    scalar getTheta(const vector &iMeshC, vector &DiskRad);
    scalar getBladeAngle(const vector &radiusVec);
    scalar getChordLength (const vector &radiusVec);
    
    void readCSV_NASA0413LS();
    scalar getLiftNASA0413(const scalar &alpha, const scalar &Re, scalar &ClField);
    scalar getDragNASA0413(const scalar &alpha, const scalar &Re, scalar &CdField);
    scalar getCl(const scalar &alpha, const std::vector<long double> &alpha_Arr, const std::vector<long double> &cl_Arr);
    scalar getCd(const scalar &alpha, const std::vector<long double> &alpha_Arr, const std::vector<long double> &cd_Arr);
    
    vector calcVolumeForce(const vector &iMeshC, const vector &radiusVec, const vector &U_up, const vector &U_down, const scalar &theta, const scalar &gamma, const scalar &chordLength, vector &relativeVelocity, scalar &betaField, scalar &alphaField, scalar &ClField, scalar &CdField, scalar &Torque, scalar &cellVolume, vector &AverageVelocity);
    
    void Post(const fvMesh &iMesh, const scalarField &Torque, const vectorField &volumeForce, const scalarField &diskID, const volScalarField &p, const volVectorField &U, const surfaceScalarField &phi, scalarField &pInlet, const scalarField &diskTheta, const scalarField &cellVolume, const vectorField &diskRadius);
    
    //Air Property Data (needs to be removed when moving to BuoyantSimpleFoam)
    const static scalar rho = 1.0857;
    const static scalar mu = 1.7948e-5;
    
    //Circle Values
    const static scalar Pi = 3.141592654;
    const static scalar deg2Rad = 3.141592654/180; 
    const static scalar rad2Deg = 180/3.141592654;
    
    std::vector<long double> alpha_2e6_Arr_0413, cl_2e6_Arr_0413, cd_2e6_Arr_0413;
    std::vector<long double> alpha_4e6_Arr_0413, cl_4e6_Arr_0413, cd_4e6_Arr_0413;
    std::vector<long double> alpha_6e6_Arr_0413, cl_6e6_Arr_0413, cd_6e6_Arr_0413;
    std::vector<long double> alpha_9e6_Arr_0413, cl_9e6_Arr_0413, cd_9e6_Arr_0413;
    
  private:
    
    vector GenesisAcDiskLocation, FlowDir, GenesisUpDiskLocation, GenesisDownDiskLocation;
    scalar EPS, HubSize,  AcDiskWidth, FanDiameter, ChordLength, DiskSize, FanSpeed, BladeNo, GamaHub;
    
  };
  
}//end namespace Foam
#endif
