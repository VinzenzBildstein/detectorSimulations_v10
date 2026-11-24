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
// $Id: DetectorMessenger.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTORMESSENGER_HH
#define DETECTORMESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger : public G4UImessenger {
public:
   explicit DetectorMessenger(DetectorConstruction*);
   ~DetectorMessenger();

   void SetNewValue(G4UIcommand*, G4String) override;

private:
   DetectorConstruction* fDetector{nullptr};

   G4UIdirectory*             fDetDir{nullptr};
   G4UIdirectory*             fAppDir{nullptr};
   G4UIdirectory*             fWorldDir{nullptr};
   G4UIdirectory*             fDetSysDir{nullptr};
   G4UIcmdWithAString*        fWorldMaterialCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fWorldDimensionsCmd{nullptr};
   G4UIcmdWithABool*          fWorldVisCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fWorldMagneticFieldCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fWorldStepLimit{nullptr};
   G4UIcmdWithoutParameter*   fUpdateCmd{nullptr};

   // Generic Target Apparatus
   G4UIcmdWithAString*        fGenericTargetCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fGenericTargetDimensionsCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fGenericTargetPositionCmd{nullptr};

   G4UIcmdWithAString* fLayeredTargetAddCmd{nullptr};

   // Spice Target Apparatus

   //field box
   G4UIcmdWithAString*        fFieldBoxMaterialCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fFieldBoxDimensionsCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fFieldBoxPositionCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fFieldBoxMagneticFieldCmd{nullptr};

   //field for SPICE
   G4UIcmdWithAString* fTabMagneticFieldCmd{nullptr};   ///19/7

   G4UIcmdWithAString*        fAddApparatusSpiceTargetChamberCmd{nullptr};
   G4UIcmdWithoutParameter*   fAddApparatus8piVacuumChamberCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fAddApparatus8piVacuumChamberAuxMatShellCmd{nullptr};
   G4UIcmdWithAnInteger*      fAddApparatusGriffinStructureCmd{nullptr};

   G4UIcmdWithAString*        fAddBoxMatCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fAddBoxThicknessCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fAddBoxInnerDimensionsCmd{nullptr};
   G4UIcmdWith3Vector*        fAddBoxColourCmd{nullptr};
   G4UIcmdWithoutParameter*   fAddBoxCmd{nullptr};

   G4UIcmdWithAString*        fAddGridMatCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fAddGridSizeCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fAddGridDimensionsCmd{nullptr};
   G4UIcmdWith3Vector*        fAddGridColourCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fAddGridPosOffsetCmd{nullptr};
   G4UIcmdWithoutParameter*   fAddGridCmd{nullptr};

   // Detection Systems
   G4UIcmdWithAnInteger* fAddDetectionSystemGammaTrackingCmd{nullptr};
   G4UIcmdWithAnInteger* fAddDetectionSystemSodiumIodideCmd{nullptr};
   G4UIcmdWith3Vector*   fAddDetectionSystemLanthanumBromideCmd{nullptr};
   G4UIcmdWith3Vector*   fAddDetectionSystemAncillaryBGOCmd{nullptr};

   G4UIcmdWithAnInteger* fAddDetectionSystem8piCmd{nullptr};
   G4UIcmdWithAnInteger* fAddDetectionSystem8piDetectorCmd{nullptr};

   G4UIcmdWithAnInteger* fAddDetectionSystemDescantCmd{nullptr};
   G4UIcmdWith3Vector*   fAddDetectionSystemDescantAuxPortsCmd{nullptr};

   G4UIcmdWithoutParameter*   fAddApparatusDescantStructureCmd{nullptr};
   G4UIcmdWithAString*        fSetDetectionSystemDescantColorCmd{nullptr};
   G4UIcmdWith3Vector*        fSetDetectionSystemDescantRotationCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fAddDetectionSystemDescantCartCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fAddDetectionSystemDescantSpherCmd{nullptr};

   G4UIcmdWith3Vector* fAddDetectionSystemTestcanCmd{nullptr};

   G4UIcmdWithAnInteger* fAddDetectionSystemSceptarCmd{nullptr};
   G4UIcmdWithAnInteger* fAddDetectionSystemGriffinForwardCmd{nullptr};
   G4UIcmdWithAnInteger* fAddDetectionSystemGriffinForwardDetectorCmd{nullptr};
   G4UIcmdWithAnInteger* fAddDetectionSystemGriffinBackCmd{nullptr};
   G4UIcmdWithAnInteger* fAddDetectionSystemGriffinBackDetectorCmd{nullptr};
   G4UIcmdWithoutParameter* fAddDetectionSystemSpiceCmd{nullptr};
   G4UIcmdWithAString*      fAddDetectionSystemTrificCmd{nullptr};
   G4UIcmdWithAString*      fAddMagnetSystemSpiceCmd{nullptr};
   G4UIcmdWithAnInteger*    fAddDetectionSystemPacesCmd{nullptr};
   G4UIcmdWithAnInteger*    fAddDetectionSystemGriffinHevimetCmd{nullptr};

   G4UIcmdWithAnInteger*      fAddDetectionSystemGriffinCustomDetectorCmd{nullptr};
   G4UIcmdWithAnInteger*      fAddDetectionSystemGriffinCustomCmd{nullptr};
   G4UIcmdWithAnInteger*      fAddDetectionSystemGriffinShieldSelectCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fAddDetectionSystemGriffinSetRadialDistanceCmd{nullptr};
   G4UIcmdWithAnInteger*      fAddDetectionSystemGriffinSetExtensionSuppLocationCmd{nullptr};
   G4UIcmdWith3Vector*        fAddDetectionSystemGriffinSetPositionCmd{nullptr};
   G4UIcmdWith3Vector*        fAddDetectionSystemGriffinSetDeadLayerCmd{nullptr};
   G4UIcmdWithABool*          fUseTIGRESSPositionsCmd{nullptr};

   G4UIcmdWithADoubleAndUnit* fTrifWinCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fTrifDegCmd{nullptr};
   G4UIcmdWithAString*        fTrifDMatCmd{nullptr};
   G4UIcmdWithABool*          fTrifAluCmd{nullptr};
   G4UIcmdWithABool*          fTrifFlatCmd{nullptr};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
