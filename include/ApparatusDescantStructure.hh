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

#ifndef APPARATUSDESCANTSTRUCTURE_HH
#define APPARATUSDESCANTSTRUCTURE_HH

#include "G4SystemOfUnits.hh"   // new version geant4.10 requires units
//#include "G4IntersectionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SubtractionSolid.hh"
#include "G4Colour.hh"

class G4AssemblyVolume;
class ApparatusDescantStructure {
public:
   ApparatusDescantStructure();
   ~ApparatusDescantStructure();

   G4int Build();

   G4int PlaceDescantStructure(G4LogicalVolume* expHallLog);

private:
   // Logical volumes
   G4LogicalVolume* fDescantStructureLog{nullptr};

   // Assembly volumes
   G4AssemblyVolume* fAssemblyDescantStructure{nullptr};

   G4String fStructureMaterial;
   G4double fStructureInnerRadius{0.};
   G4double fStructureOuterRadius{0.};
   G4double fStartPhi{0.};
   G4double fDeltaPhi{0.};
   G4double fStartTheta{0.};
   G4double fDeltaTheta{0.};
   G4double fStructureRadialDistance{0.};

   //for subtraction polyhedra
   G4int    fNumSides{0};
   G4int    fNumZplanes{0};
	std::array<G4double, 2> fZPlanes;

   G4double fPolyRadiusRed{0.};
	std::array<G4double, 2> fInnerRadiusRed;
	std::array<G4double, 2> fOuterRadiusRed;

   G4double fPolyRadiusGreenYellow{0.};
	std::array<G4double, 2> fInnerRadiusGreenYellow;
	std::array<G4double, 2> fOuterRadiusGreenYellow;

   G4double fPolyRadiusWhite{0.};
	std::array<G4double, 2> fInnerRadiusWhite;
	std::array<G4double, 2> fOuterRadiusWhite;

   //for the subtraction box
   G4double fSubBoxX{0.};
   G4double fSubBoxY{0.};
   G4double fSubBoxZ{0.};
   G4double fMoveCut{0.};

   //for back polyhedra
   G4double fBackRadius{0.};
	std::array<G4double, 2> fInnerRadiusBack;
	std::array<G4double, 2> fOuterRadiusBack;
	std::array<G4double, 2> fZPlanesBack;

   // for subtraction cylinder
   G4double fMinRadius{0.};
   G4double fMaxRadius{0.};
   G4double fHalfLength{0.};
   G4double fSubtractionStartPhi{0.};
   G4double fSubtractionDeltaPhi{0.};

   // The Euler angles from James' MSc thesis which gives us the detector positions
   // Some of the angles for the green and yellow detectors are wrong in James' thesis,
   // note the +180 on a few angles.
	std::array<std::array<G4double, 3>, 15> fBlueAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 10> fGreenAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 15> fRedAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 20> fWhiteAlphaBetaGamma;
   std::array<std::array<G4double, 3>, 10> fYellowAlphaBetaGamma;

   // Descant Structure Colour
   G4Colour fGreyColour;

   G4int BuildDescantShell();
};

#endif
