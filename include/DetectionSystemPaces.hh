//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectionSystemPaces.hh,v 1.1 2012-11-14 14:00:00 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTIONSYSTEMPACES_HH
#define DETECTIONSYSTEMPACES_HH

#include "G4AssemblyVolume.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class DetectionSystemPaces {
public:
   DetectionSystemPaces();
   ~DetectionSystemPaces();

   G4int Build();   //G4SDManager* mySDman);
   G4int PlaceDetector(G4LogicalVolume* expHallLog, G4int ndet);

private:
   // Assembly volumes
   G4AssemblyVolume* fAssembly{nullptr};
   G4AssemblyVolume* fAssemblyDetector{nullptr};
   G4AssemblyVolume* fAssemblySilicon{nullptr};

   // Logical volumes
   G4LogicalVolume* fAluminumHemisphereLog{nullptr};
   G4LogicalVolume* fAluminumAnnulusTopLog{nullptr};
   G4LogicalVolume* fAluminumAnnulusBotLog{nullptr};
   G4LogicalVolume* fCanisterLog{nullptr};
   G4LogicalVolume* fSiliconBlockLog{nullptr};
   G4LogicalVolume* fSiliconDeadLayerLog{nullptr};
   G4LogicalVolume* fTeflonAnnulusTopLog{nullptr};
   G4LogicalVolume* fTeflonAnnulusBotLog{nullptr};
   G4LogicalVolume* fDelrinHemisphereLog{nullptr};

   G4double fCutClearance;

   // Solid parameters
   G4double fAluminumHemisphereInnerRadius;
   G4double fAluminumHemisphereOuterRadius;
   G4double fAluminumHemisphereBeamHoleRadius;
   G4double fAluminumHemisphereBeamHoleRimHeight;
   G4double fAluminumAnnulusTopInnerRadius;
   G4double fAluminumAnnulusTopOuterRadius;
   G4double fAluminumAnnulusTopThickness;
   G4double fAluminumAnnulusBotInnerRadius;
   G4double fAluminumAnnulusBotOuterRadius;
   G4double fAluminumAnnulusBotThickness;
   G4double fCanisterInnerRadius;
   G4double fCanisterOuterRadius;
   G4double fCanisterThickness;
   G4double fSiliconBlockRadius;
   G4double fSiliconBlockThickness;
   G4double fSiliconDeadLayerThickness;
   G4double fScrewRadius;
   G4double fScrewPlacementRadius;
   //G4double screwLength;
   G4double fTeflonAnnulusTopInnerRadius;
   G4double fTeflonAnnulusTopOuterRadius;
   G4double fTeflonAnnulusTopThickness;
   G4double fTeflonAnnulusBotInnerRadius;
   G4double fTeflonAnnulusBotOuterRadius;
   G4double fTeflonAnnulusBotThickness;
   G4double fDelrinHemisphereInnerRadius;
   G4double fDelrinHemisphereOuterRadius;
   G4double fDelrinHemisphereBeamHoleRadius;

   G4double fAluminumHemisphereDist;
   G4double fDelrinHemisphereDist;

   // distances from FRONT of detector facing the source, for assembly
   G4double fAluminumAnnulusTopDist;   // will be half the annulus thickness
   G4double fAluminumAnnulusBotDist;
   G4double fCanisterDist;
   G4double fSiliconBlockDist;
   G4double fSiliconDeadLayerFrontDist;
   G4double fSiliconDeadLayerBackDist;
   G4double fTeflonAnnulusTopDist;
   G4double fTeflonAnnulusBotDist;

   // detector placement and orientation scalars
	std::array<G4double, 5> fPacesPlacementDistance;
   std::array<G4double, 5> fPacesPlacementTheta;
   std::array<G4double, 5> fPacesPlacementPhi;
   std::array<G4double, 5> fPacesOrientationTheta;
   std::array<G4double, 5> fPacesOrientationPhi;

   // assembly volume sizes
   G4double fDetectorAssemblyRadius;
   G4double fDetectorAssemblyThickness;

   // Construction methods
   G4int AddAluminumHemisphere();
   G4int AddAluminumAnnulusTop();
   G4int AddAluminumAnnulusBot();
   G4int AddCanister();
   G4int AddSiliconBlock();
   G4int AddSiliconDeadLayer();
   G4int AddScrews() const;
   G4int AddTeflonAnnulusTop();
   G4int AddTeflonAnnulusBot();
   G4int AddDelrinHemisphere();
   G4int CombineAssemblySilicon();
   G4int CombineAssemblyDetector();
   G4int CombineAssembly();
};

#endif
