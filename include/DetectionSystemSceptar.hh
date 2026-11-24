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

#ifndef DETECTIONSYSTEMSCEPTAR_HH
#define DETECTIONSYSTEMSCEPTAR_HH

#include "G4SubtractionSolid.hh"
#include "G4AssemblyVolume.hh"
#include "G4Trd.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectionSystemSceptar {
public:
   DetectionSystemSceptar();
   ~DetectionSystemSceptar();

   G4int Build();   //G4SDManager* mySDman);
   G4int PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber);

private:
   // Assembly volumes
   G4AssemblyVolume* fAssembly{nullptr};
   G4AssemblyVolume* fAssemblySquare{nullptr};
   G4AssemblyVolume* fAssemblyAngled{nullptr};
   G4AssemblyVolume* fAssemblySquareSD{nullptr};
   G4AssemblyVolume* fAssemblyAngledSD{nullptr};

   G4double fConvert;
   G4double fSquareScintillatorLength;
   G4double fSquareScintillatorWidth;
   G4double fSquareScintillatorThickness;
   G4double fAngledScintillatorLength;
   G4double fAngledScintillatorLongWidth;
   G4double fAngledScintillatorShortWidth;
   G4double fAngledScintillatorThickness;
   G4double fMylarThickness;
   G4double fScintGap;
   G4double fScintAngle1;
   G4double fScintAngle2;
   G4double fScintAngle3;
   G4double fScintAngle4;
   G4double fScintAngleMove;
   G4double fSquareScintRadialDistance;
   G4double fAngledScintRadialDistance;
   G4double fAngledScintMoveBack;
   G4double fDelrinInnerRadius;
   G4double fDelrinOuterRadius;
   G4double fDelrin2InnerRadius;
   G4double fDelrin2OuterRadius;
   G4double fHevimetInnerRadius;
   G4double fHevimetOuterRadius;
   G4double fDelrinHoleRadius;
   G4double fSeparateHemispheres;

   //##################################################################
   //### methods to construct all of the components of the detector ###
   //##################################################################
   G4int ConstructScintillator();
   G4int ConstructDelrinShell();
   G4int Construct2ndDelrinShell();
   G4int ConstructHevimetShell();

   //Logical Volumes
   //
   G4LogicalVolume* fSquareMylarLog{nullptr};
   G4LogicalVolume* fAngledMylarLog{nullptr};
   G4LogicalVolume* fSquareScintillatorLog{nullptr};
   G4LogicalVolume* fAngledScintillatorLog{nullptr};
   G4LogicalVolume* fDelrinShellLog{nullptr};
   G4LogicalVolume* fDelrinShell2Log{nullptr};
   G4LogicalVolume* fHevimetShellLog{nullptr};

   G4Trd*              SquareMylar() const;
   G4SubtractionSolid* SquareMylarWithCut() const;
   G4SubtractionSolid* AngledMylar() const;
   G4SubtractionSolid* AngledMylarWithCut() const;
   G4Trd*              SquareScintillator() const;
   G4SubtractionSolid* AngledScintillator() const;
   G4SubtractionSolid* DelrinShell() const;
   G4SubtractionSolid* DelrinShell2() const;
   G4SubtractionSolid* HevimetShell() const;
};

#endif
