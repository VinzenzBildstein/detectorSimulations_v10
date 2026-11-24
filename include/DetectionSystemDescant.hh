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

#ifndef DETECTIONSYSTEMDESCANT_HH
#define DETECTIONSYSTEMDESCANT_HH

#include "G4SystemOfUnits.hh"   // new version geant4.10 requires units
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SubtractionSolid.hh"
#include "G4Colour.hh"

class G4AssemblyVolume;

class DetectionSystemDescant {
public:
   explicit DetectionSystemDescant(G4bool leadShield);
   ~DetectionSystemDescant();

   G4int Build();
   G4int PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber);
   G4int PlaceDetectorAuxPorts(G4LogicalVolume* expHallLog, G4int detectorNumber, G4double radialpos);

	G4int PlaceDetector(G4LogicalVolume* expHallLog, const G4String& color, G4ThreeVector pos, const G4ThreeVector& rot);

private:
   // Logical volumes
   G4LogicalVolume* fBlueVolumeLog{nullptr};
   G4LogicalVolume* fGreenVolumeLog{nullptr};
   G4LogicalVolume* fRedVolumeLog{nullptr};
   G4LogicalVolume* fWhiteVolumeLog{nullptr};
   G4LogicalVolume* fYellowVolumeLog{nullptr};

   G4LogicalVolume* fBlueScintillatorVolumeLog{nullptr};
   G4LogicalVolume* fGreenScintillatorVolumeLog{nullptr};
   G4LogicalVolume* fRedScintillatorVolumeLog{nullptr};
   G4LogicalVolume* fWhiteScintillatorVolumeLog{nullptr};
   G4LogicalVolume* fYellowScintillatorVolumeLog{nullptr};

   G4LogicalVolume* fBlueLeadVolumeLog{nullptr};
   G4LogicalVolume* fGreenLeadVolumeLog{nullptr};
   G4LogicalVolume* fRedLeadVolumeLog{nullptr};
   G4LogicalVolume* fWhiteLeadVolumeLog{nullptr};
   G4LogicalVolume* fYellowLeadVolumeLog{nullptr};

   G4LogicalVolume* fQuartzWindow5inchLog{nullptr};
   G4LogicalVolume* fQuartzWindow3inchLog{nullptr};

   // Assembly volumes
   G4AssemblyVolume* fAssemblyBlue{nullptr};               // Contains all non-sensitive materials
   G4AssemblyVolume* fAssemblyBlueScintillator{nullptr};   // Contains all only sensitive materials, eg. the scintillation material
   G4AssemblyVolume* fAssemblyGreen{nullptr};
   G4AssemblyVolume* fAssemblyGreenScintillator{nullptr};
   G4AssemblyVolume* fAssemblyRed{nullptr};
   G4AssemblyVolume* fAssemblyRedScintillator{nullptr};
   G4AssemblyVolume* fAssemblyWhite{nullptr};
   G4AssemblyVolume* fAssemblyWhiteScintillator{nullptr};
   G4AssemblyVolume* fAssemblyYellow{nullptr};
   G4AssemblyVolume* fAssemblyYellowScintillator{nullptr};

   G4double fCanLength;
   G4double fCanThickness;
   G4double fCanThicknessFront;
   G4double fCanInnerRadius;
   G4double fLeadShieldThickness;
   G4double fRadialDistance;
   G4String fCanMaterial;
   G4String fLiquidMaterial;
   G4String fLeadMaterial;
   G4bool   fIncludeLead;

   G4double fCanBackThickness;

   // for making PMT cuts
   G4double fStartPhi;
   G4double fEndPhi;
   G4double fInnerRadiusWindow;
   G4double fOuterRadiusWindow5inch;
   G4double fOuterRadiusWindow3inch;
   G4double fHalfLengthZWindowCut;
   G4double fWhitePMTOffsetY;
   G4double fBluePMTOffsetX;
   G4double fRedPMTOffsetX;
   G4double fYellowGreenPMTOffsetX;
   G4double fYellowGreenPMTOffsetY;
   G4String fQuartzMaterial;
   G4double fOpticalWindowHalfThickness;

   // Saint Gobain data files, 6 points around the front face of the can, and 6 on the back
   // from Dan Brennan
	std::array<std::array<G4double, 3>, 12> fBlueDetector;
   std::array<std::array<G4double, 3>, 12> fGreenDetector;
   std::array<std::array<G4double, 3>, 12> fRedDetector;
   std::array<std::array<G4double, 3>, 12> fWhiteDetector;
   std::array<std::array<G4double, 3>, 12> fYellowDetector;

   G4double fTrimGreenY;
   G4double fTrimYellowY;
   G4double fTrimBlueX;

   // These are the angles of the cuts for each of the 6 sides of the detectors
	std::array<G4double, 6> fBluePhi;
   std::array<G4double, 6> fGreenPhi;
   std::array<G4double, 6> fRedPhi;
   std::array<G4double, 6> fWhitePhi;
   std::array<G4double, 6> fYellowPhi;

   // The Euler angles from James' MSc thesis which gives us the detector positions
   // Some of the angles for the green and yellow detectors are wrong in James' thesis,
   // note the +180 on a few angles.
	std::array<std::array<G4double, 3>, 15> fBlueAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 10> fGreenAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 15> fRedAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 20> fWhiteAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 10> fYellowAlphaBetaGamma;

   // for LaBr3 detector locations
	std::array<std::array<G4double, 5>, 8> fDetectorAngles;
   G4double fSetRadialPos;

   // The colours of the detectors
   G4Colour fBlueColour;
   G4Colour fGreenColour;
   G4Colour fRedColour;
   G4Colour fWhiteColour;
   G4Colour fYellowColour;
   G4Colour fLiquidColour;    // Scintillator colour
   G4Colour fGreyColour;      // stainless steel pmt tube colour
   G4Colour fMagentaColour;   // optical window colour
   G4Colour fBlackColour;

   G4bool fSurfCheck;

   G4int BuildCanVolume();
   G4int BuildDetectorVolume();

	G4SubtractionSolid* CanVolume(G4bool insideVol, G4double volumeLength, std::array<std::array<G4double, 3>, 12> detector, std::array<G4double, 6> detectorPhi);
   G4SubtractionSolid* CutVolumeOnFourPoints(G4int idx, G4bool insideVol, G4double volumeLength, std::array<G4double, 6> detectorPhi, G4SubtractionSolid* volume, const G4ThreeVector& frontP1, const G4ThreeVector& frontP2);

   G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

   G4ThreeVector SolveLineEquationX(const G4ThreeVector& p1, const G4ThreeVector& p2, G4double x);
   G4ThreeVector SolveLineEquationY(const G4ThreeVector& p1, const G4ThreeVector& p2, G4double y);
   G4ThreeVector SolveLineEquationZ(const G4ThreeVector& p1, const G4ThreeVector& p2, G4double z);
};

#endif
