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
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTIONSYSTEMLANTHANUMBROMIDE_HH
#define DETECTIONSYSTEMLANTHANUMBROMIDE_HH

#include "G4SystemOfUnits.hh"   // new version geant4.10 requires units
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

const G4double inchtocm = 2.54;

class DetectionSystemLanthanumBromide {
public:
   DetectionSystemLanthanumBromide();
   ~DetectionSystemLanthanumBromide();

   G4int    Build();   //G4SDManager* mySDman);
   G4int    PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber, G4double radialpos);
   G4double GetDetectorLengthOfUnitsCM() const { return fDetectorLengthZ; };
   G4double GetCrystalRadius() const { return fCrystalOuterRadius; };
   G4double GetCrystalLength() const { return fCrystalLengthZ; };
   G4double GetR() const;
   G4double GetTheta(G4int i) const;
   G4double GetPhi(G4int i) const;
   G4double GetYaw(G4int i) const;
   G4double GetPitch(G4int i) const;
   G4double GetRoll(G4int i) const;
   G4double GetCrystalRadialPosition() const { return fPackingFrontLidThickness + fDiscFrontLidThickness + fSealFrontLidThickness + fCanFrontLidThickness + 12.5 * cm; };

private:
   // Logical volumes
   G4LogicalVolume* fDetectorVolumeLog{nullptr};
   G4LogicalVolume* fCrystalBlockLog{nullptr};

   G4LogicalVolume* fCanCylinderLog{nullptr};
   G4LogicalVolume* fCanFrontLidLog{nullptr};
   G4LogicalVolume* fCanBackLidLog{nullptr};

   G4LogicalVolume* fSealFrontLidLog{nullptr};
   G4LogicalVolume* fDiscFrontLidLog{nullptr};

   G4LogicalVolume* fPackingCylinderLog{nullptr};
   G4LogicalVolume* fPackingFrontLidLog{nullptr};

   // Assembly volumes
   G4AssemblyVolume* fAssembly{nullptr};

   //    SensitiveDetector* crystalBlockSD;

   G4int fCopyNumber;
   G4int fNumberOfDetectors;
   G4int fNumberOfSegments;

   G4String fCrystalMaterial;
   G4String fCanMaterial;
   G4String fSealMaterial;
   G4String fDiscMaterial;
   G4String fPackingMaterial;

   G4double fDetailViewEndAngle;
   G4double fCrystalLengthZ;
   G4double fCrystalInnerRadius;
   G4double fCrystalOuterRadius;
   G4double fPackingLengthZ;
   G4double fPackingInnerRadius;
   G4double fPackingOuterRadius;
   G4double fPackingLidInnerRadius;
   G4double fPackingLidOuterRadius;
   G4double fPackingFrontLidThickness;
   G4double fDiscLidInnerRadius;
   G4double fDiscLidOuterRadius;
   G4double fDiscFrontLidThickness;
   G4double fSealLidInnerRadius;
   G4double fSealLidOuterRadius;
   G4double fSealFrontLidThickness;
   G4double fCanLengthZ;
   G4double fCanInnerRadius;
   G4double fCanOuterRadius;
   G4double fCanLidInnerRadius;
   G4double fCanLidOuterRadius;
   G4double fCanFrontLidThickness;
   G4double fCanBackLidThickness;

   G4double fDetectorLengthZ;
	std::array<std::array<G4double, 5>, 8> fDetectorAngles;

   G4double fSetRadialPos;

   G4Tubs* BuildCrystal() const;
   G4Tubs* BuildAluminumCan() const;
   G4Tubs* BuildAluminumCanFrontLid() const;
   G4Tubs* BuildAluminumCanBackLid() const;
   G4Tubs* BuildPacking() const;
   G4Tubs* BuildPackingFrontLid() const;
   G4Tubs* BuildDiscFrontLid() const;
   G4Tubs* BuildSealFrontLid() const;

   G4int BuildSealVolume();
   G4int BuildDiscVolume();
   G4int BuildPackingVolume();
   G4int BuildAluminumCanVolume();
   G4int BuildCrystalVolume();
   G4int BuildOneDetector();

   G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi) const;
};

#endif
