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
// $Id: DetectorMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADouble.hh"

#include <string>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
   : fDetector(Det), fDetSysDir(new G4UIdirectory("/DetSys/")),
   fDetDir(new G4UIdirectory("/DetSys/det/")),
   fAppDir(new G4UIdirectory("/DetSys/app/")),
   fWorldDir(new G4UIdirectory("/DetSys/world/")),
   fWorldMaterialCmd(new G4UIcmdWithAString("/DetSys/world/material", this)),
   fWorldDimensionsCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/world/dimensions", this)),
   fWorldVisCmd(new G4UIcmdWithABool("/DetSys/world/vis", this)),
   fWorldMagneticFieldCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/world/magneticField", this)),
   fWorldStepLimit(new G4UIcmdWithADoubleAndUnit("/DetSys/world/StepLimit", this)),
   fGenericTargetCmd(new G4UIcmdWithAString("/DetSys/app/genericTarget", this)),
   fGenericTargetDimensionsCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/app/genericTargetDimensions", this)),
   fGenericTargetPositionCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/app/genericTargetPosition", this)),
   fLayeredTargetAddCmd(new G4UIcmdWithAString("/DetSys/app/LayeredTargetAddLayer", this)),
   fFieldBoxMaterialCmd(new G4UIcmdWithAString("/DetSys/app/fieldBoxMaterial", this)),
   fFieldBoxDimensionsCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxDimensions", this)),
   fFieldBoxPositionCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxPosition", this)),
   fFieldBoxMagneticFieldCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxMagneticField", this)),
   fTabMagneticFieldCmd(new G4UIcmdWithAString("/DetSys/world/TabMagneticField", this)),   ///19/7
   fAddBoxMatCmd(new G4UIcmdWithAString("/DetSys/det/boxMat", this)),
   fAddBoxThicknessCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/det/boxThickness", this)),
   fAddBoxInnerDimensionsCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/det/boxInnerDimensions", this)),
   fAddBoxColourCmd(new G4UIcmdWith3Vector("/DetSys/det/boxColour", this)),
   fAddBoxCmd(new G4UIcmdWithoutParameter("/DetSys/det/addBox", this)),
   fAddGridMatCmd(new G4UIcmdWithAString("/DetSys/det/gridMat", this)),
   fAddGridSizeCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/det/gridSize", this)),
   fAddGridDimensionsCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/det/gridDimensions", this)),
   fAddGridColourCmd(new G4UIcmdWith3Vector("/DetSys/det/gridColour", this)),
   fAddGridCmd(new G4UIcmdWithoutParameter("/DetSys/det/addGrid", this)),
   fAddGridPosOffsetCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/det/gridPosOffset", this)),
   fAddApparatusSpiceTargetChamberCmd(new G4UIcmdWithAString("/DetSys/app/addSpiceTargetChamber", this)),
   fAddApparatus8piVacuumChamberCmd(new G4UIcmdWithoutParameter("/DetSys/app/add8piVacuumChamber", this)),
   fAddApparatus8piVacuumChamberAuxMatShellCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/app/add8piVacuumChamberAuxMatShell", this)),
   fAddApparatusGriffinStructureCmd(new G4UIcmdWithAnInteger("/DetSys/app/addGriffinStructure", this)),
   fUpdateCmd(new G4UIcmdWithoutParameter("/DetSys/det/update", this)),
   fAddDetectionSystemGammaTrackingCmd(new G4UIcmdWithAnInteger("/DetSys/det/addGammaTracking", this)),
   fAddDetectionSystemSodiumIodideCmd(new G4UIcmdWithAnInteger("/DetSys/det/addSodiumIodide", this)),
   fAddDetectionSystemLanthanumBromideCmd(new G4UIcmdWith3Vector("/DetSys/det/addLanthanumBromide", this)),
   fAddDetectionSystemAncillaryBGOCmd(new G4UIcmdWith3Vector("/DetSys/det/addAncillaryBGO", this)),
   fAddDetectionSystem8piCmd(new G4UIcmdWithAnInteger("/DetSys/det/add8pi", this)),
   fAddDetectionSystem8piDetectorCmd(new G4UIcmdWithAnInteger("/DetSys/det/add8piDetector", this)),
   fAddDetectionSystemDescantCmd(new G4UIcmdWithAnInteger("/DetSys/det/addDescant", this)),
   fAddDetectionSystemDescantAuxPortsCmd(new G4UIcmdWith3Vector("/DetSys/det/addDescantAuxPorts", this)),
   fAddApparatusDescantStructureCmd(new G4UIcmdWithoutParameter("/DetSys/det/addDescantStructure", this)),
   fAddDetectionSystemTestcanCmd(new G4UIcmdWith3Vector("/DetSys/det/addTestcan", this)),
   fSetDetectionSystemDescantColorCmd(new G4UIcmdWithAString("/DetSys/det/setDescantColor", this)),
   fSetDetectionSystemDescantRotationCmd(new G4UIcmdWith3Vector("/DetSys/det/setDescantRotation", this)),
   fAddDetectionSystemDescantCartCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/det/addDescantCart", this)),
   fAddDetectionSystemDescantSpherCmd(new G4UIcmdWith3VectorAndUnit("/DetSys/det/addDescantSpher", this)),
   fAddDetectionSystemGriffinForwardCmd(new G4UIcmdWithAnInteger("/DetSys/det/addGriffinForward", this)),
   fAddDetectionSystemGriffinForwardDetectorCmd(new G4UIcmdWithAnInteger("/DetSys/det/addGriffinForwardDetector", this)),
   fAddDetectionSystemGriffinBackCmd(new G4UIcmdWithAnInteger("/DetSys/det/addGriffinBack", this)),
   fAddDetectionSystemGriffinBackDetectorCmd(new G4UIcmdWithAnInteger("/DetSys/det/addGriffinBackDetector", this)),
   fAddDetectionSystemGriffinCustomDetectorCmd(new G4UIcmdWithAnInteger("/DetSys/det/addGriffinCustomDetector", this)),
   fAddDetectionSystemGriffinCustomCmd(new G4UIcmdWithAnInteger("/DetSys/det/addGriffinCustom", this)),
   fAddDetectionSystemGriffinShieldSelectCmd(new G4UIcmdWithAnInteger("/DetSys/det/SetCustomShieldsPresent", this)),
   fAddDetectionSystemGriffinSetRadialDistanceCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/det/SetCustomRadialDistance", this)),
   fAddDetectionSystemGriffinSetExtensionSuppLocationCmd(new G4UIcmdWithAnInteger("/DetSys/det/SetCustomExtensionSuppressorLocation", this)),
   fAddDetectionSystemGriffinSetPositionCmd(new G4UIcmdWith3Vector("/DetSys/det/SetCustomPosition", this)),
   fAddDetectionSystemGriffinSetDeadLayerCmd(new G4UIcmdWith3Vector("/DetSys/det/SetCustomDeadLayer", this)),
   fAddDetectionSystemGriffinHevimetCmd(new G4UIcmdWithAnInteger("/DetSys/det/includeGriffinHevimet", this)),
   fAddDetectionSystemSceptarCmd(new G4UIcmdWithAnInteger("/DetSys/det/addSceptar", this)),
   fAddDetectionSystemSpiceCmd(new G4UIcmdWithoutParameter("/DetSys/det/addSpiceDetector", this)),
   fAddDetectionSystemTrificCmd(new G4UIcmdWithAString("/DetSys/det/addTrificDetector", this)),
   fAddDetectionSystemPacesCmd(new G4UIcmdWithAnInteger("/DetSys/det/addPaces", this)),
   fUseTIGRESSPositionsCmd(new G4UIcmdWithABool("/DetSys/det/UseTIGRESSPositions", this)),
   fTrifFlatCmd(new G4UIcmdWithABool("/DetSys/det/TrificFlatWindow", this)),
   fTrifAluCmd(new G4UIcmdWithABool("/DetSys/det/TrificAluWindow", this)),
   fTrifWinCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/det/TrificWindowThickness", this)),
   fTrifDegCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/det/TrificDegraderThickness", this)),
   fTrifDMatCmd(new G4UIcmdWithAString("/DetSys/det/TrificDegraderMaterial", this))
{
   fDetSysDir->SetGuidance("UI commands of this simulation");

   fDetDir->SetGuidance("detector control");

   fAppDir->SetGuidance("apparatus control");

   fWorldDir->SetGuidance("world control");

   fWorldMaterialCmd->SetGuidance("Select material for the world.");
   fWorldMaterialCmd->SetParameterName("choice", false);
   fWorldMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fWorldDimensionsCmd->SetGuidance("Set world dimensions - x y z unit.");
   fWorldDimensionsCmd->SetUnitCategory("Length");
   fWorldDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fWorldVisCmd->SetGuidance("Set the visulization of the world");
   fWorldVisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fWorldMagneticFieldCmd->SetGuidance("Set world magnetic field - x y z unit.");
   fWorldMagneticFieldCmd->SetUnitCategory("Magnetic flux density");
   fWorldMagneticFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fWorldStepLimit->SetGuidance("Set user step limit for the world volume.");
   fWorldStepLimit->SetUnitCategory("Length");
   fWorldStepLimit->AvailableForStates(G4State_PreInit, G4State_Idle);

   fGenericTargetCmd->SetGuidance("Select material of the target.");
   fGenericTargetCmd->SetParameterName("choice", false);
   fGenericTargetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fGenericTargetDimensionsCmd->SetGuidance("Set target dimensions - x y z unit.");
   fGenericTargetDimensionsCmd->SetUnitCategory("Length");
   fGenericTargetDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fGenericTargetPositionCmd->SetGuidance("Set target position - x y z unit.");
   fGenericTargetPositionCmd->SetUnitCategory("Length");
   fGenericTargetPositionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fLayeredTargetAddCmd->SetGuidance("Add layer to target - material mg/cm3 mg/cm2.");
   fLayeredTargetAddCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fFieldBoxMaterialCmd->SetGuidance("Select material of the target.");
   fFieldBoxMaterialCmd->SetParameterName("choice", false);
   fFieldBoxMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fFieldBoxDimensionsCmd->SetGuidance("Set target dimensions - x y z unit.");
   fFieldBoxDimensionsCmd->SetUnitCategory("Length");
   fFieldBoxDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fFieldBoxPositionCmd->SetGuidance("Set target position - x y z unit.");
   fFieldBoxPositionCmd->SetUnitCategory("Length");
   fFieldBoxPositionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fFieldBoxMagneticFieldCmd->SetGuidance("Set magnetic field - x y z unit.");
   fFieldBoxMagneticFieldCmd->SetUnitCategory("Magnetic flux density");
   fFieldBoxMagneticFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fTabMagneticFieldCmd->SetGuidance("Set tabulated magnetic field.");
   fTabMagneticFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   // Box Stuff
   fAddBoxMatCmd->SetGuidance("Set box material.");
   fAddBoxMatCmd->SetParameterName("choice", false);
   fAddBoxMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddBoxThicknessCmd->SetGuidance("Set box thickness.");
   fAddBoxThicknessCmd->SetUnitCategory("Length");
   fAddBoxThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddBoxInnerDimensionsCmd->SetGuidance("Set box inner dimensions.");
   fAddBoxInnerDimensionsCmd->SetUnitCategory("Length");
   fAddBoxInnerDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddBoxColourCmd->SetGuidance("Set box colour.");
   fAddBoxColourCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddBoxCmd->SetGuidance("Add a box.");
   fAddBoxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   // Grid Stuff
   fAddGridMatCmd->SetGuidance("Set grid material.");
   fAddGridMatCmd->SetParameterName("choice", false);
   fAddGridMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddGridSizeCmd->SetGuidance("Set grid size.");
   fAddGridSizeCmd->SetUnitCategory("Length");
   fAddGridSizeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddGridDimensionsCmd->SetGuidance("Set grid dimensions.");
   fAddGridDimensionsCmd->SetUnitCategory("Length");
   fAddGridDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddGridColourCmd->SetGuidance("Set grid colour.");
   fAddGridColourCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddGridCmd->SetGuidance("Add a grid.");
   fAddGridCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddGridPosOffsetCmd->SetGuidance("Set grid offset.");
   fAddGridPosOffsetCmd->SetUnitCategory("Length");
   fAddGridPosOffsetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddApparatusSpiceTargetChamberCmd->SetGuidance("Add SPICE target chamber.");
   fAddApparatusSpiceTargetChamberCmd->AvailableForStates(G4State_Idle);

   fAddApparatus8piVacuumChamberCmd->SetGuidance("Add 8pi vacuum chamber.");
   fAddApparatus8piVacuumChamberCmd->AvailableForStates(G4State_Idle);

   fAddApparatus8piVacuumChamberAuxMatShellCmd->SetGuidance("Add AuxMat shell around 8pi vacuum chamber");
   fAddApparatus8piVacuumChamberAuxMatShellCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddApparatusGriffinStructureCmd->SetGuidance("Add Griffin Structure");
   fAddApparatusGriffinStructureCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fUpdateCmd->SetGuidance("Update geometry.");
   fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
   fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
   fUpdateCmd->AvailableForStates(G4State_Idle);

   fAddDetectionSystemGammaTrackingCmd->SetGuidance("Add Detection System GammaTracking");
   fAddDetectionSystemGammaTrackingCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemSodiumIodideCmd->SetGuidance("Add Detection System SodiumIodide");
   fAddDetectionSystemSodiumIodideCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemLanthanumBromideCmd->SetGuidance("Add Detection System LanthanumBromide - number of dets, radius in cm, empty");
   fAddDetectionSystemLanthanumBromideCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemAncillaryBGOCmd->SetGuidance("Add Detection System AncillaryBGO - number of dets, radius in cm, hevimet yes(1) no(0)");
   fAddDetectionSystemAncillaryBGOCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystem8piCmd->SetGuidance("Add Detection System 8pi");
   fAddDetectionSystem8piCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystem8piDetectorCmd->SetGuidance("Add 8pi Detector");
   fAddDetectionSystem8piDetectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   fAddDetectionSystemDescantCmd->SetGuidance("Add Detection System DESCANT");
   fAddDetectionSystemDescantCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   fAddDetectionSystemDescantAuxPortsCmd->SetGuidance("Add 8 DESCANT detectors in the auxillary LaBr3 detector locations");
   fAddDetectionSystemDescantAuxPortsCmd->SetGuidance("/DetSys/det/addDescantAuxPorts _nDet_ _radialPos_ _leadShield_");
   fAddDetectionSystemDescantAuxPortsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddApparatusDescantStructureCmd->SetGuidance("Add DESCANT structure");
   fAddApparatusDescantStructureCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemTestcanCmd->SetGuidance("Add Testcan Detection System");
   fAddDetectionSystemTestcanCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fSetDetectionSystemDescantColorCmd->SetGuidance("Set color of next descant detector to be added via addDescantCart or addDescantSpher");
   fSetDetectionSystemDescantColorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
   fSetDetectionSystemDescantColorCmd->SetCandidates("white blue red green yellow");

   fSetDetectionSystemDescantRotationCmd->SetGuidance("Set rotation of next descant detector to be added via addDescantCart or addDescantSpher (alhpa beta gamma)");
   fSetDetectionSystemDescantRotationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemDescantCartCmd->SetGuidance("Add single DESCANT detector at provided cartesian coordinates");
   fAddDetectionSystemDescantCartCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemDescantSpherCmd->SetGuidance("Add single DESCANT detector at provided spherical coordinates");
   fAddDetectionSystemDescantSpherCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinForwardCmd->SetGuidance("Add Detection System GriffinForward");
   fAddDetectionSystemGriffinForwardCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinForwardDetectorCmd->SetGuidance("Add GriffinForward Detector");
   fAddDetectionSystemGriffinForwardDetectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinBackCmd->SetGuidance("Add Detection System GriffinBack");
   fAddDetectionSystemGriffinBackCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinBackDetectorCmd->SetGuidance("Add GriffinBack Detector");
   fAddDetectionSystemGriffinBackDetectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinCustomDetectorCmd->SetGuidance("Adds a detector using the paramaters specified");
   fAddDetectionSystemGriffinCustomDetectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinCustomCmd->SetGuidance("Adds a detection system using the paramaters specified");
   fAddDetectionSystemGriffinCustomCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   //////// Commands that are required for addGriffinCustom
   fAddDetectionSystemGriffinShieldSelectCmd->SetGuidance("Selects whether or not the detector suppressors are included");
   fAddDetectionSystemGriffinShieldSelectCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinSetRadialDistanceCmd->SetGuidance("Selects the radial distance for the detector from the origin");
   fAddDetectionSystemGriffinSetRadialDistanceCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinSetExtensionSuppLocationCmd->SetGuidance("Selects a position for the extension suppressors. Either forward (0) or back (1).");
   fAddDetectionSystemGriffinSetExtensionSuppLocationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinSetPositionCmd->SetGuidance("Sets the position and number for the detector placed in the next call to addGriffinCustom.");
   fAddDetectionSystemGriffinSetPositionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemGriffinSetDeadLayerCmd->SetGuidance("Sets the dead layer for the specified crystal.");
   fAddDetectionSystemGriffinSetDeadLayerCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
   ////////

   fAddDetectionSystemGriffinHevimetCmd->SetGuidance("Includes the Hevimet for a Griffin detector.");
   fAddDetectionSystemGriffinHevimetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemSceptarCmd->SetGuidance("Add Detection System Sceptar");
   fAddDetectionSystemSceptarCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemSpiceCmd->SetGuidance("Add Detection System Spice");
   fAddDetectionSystemSpiceCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemTrificCmd->SetGuidance("Add Detection System Trific");
   fAddDetectionSystemTrificCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAddDetectionSystemPacesCmd->SetGuidance("Add Detection System Paces");
   fAddDetectionSystemPacesCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fUseTIGRESSPositionsCmd->SetGuidance("Use TIGRESS detector positions rather than GRIFFIN");
   fUseTIGRESSPositionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fTrifFlatCmd->SetGuidance("Set the window to be flat for future trific instances");
   fTrifFlatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fTrifAluCmd->SetGuidance("Set the window to be aluminsed for future trific instances");
   fTrifAluCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fTrifWinCmd->SetGuidance("Set the window thickness for future trific instances");
   fTrifWinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fTrifDegCmd->SetGuidance("Set the thickness of trific upstream degrader");
   fTrifDegCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fTrifDMatCmd->SetGuidance("Select material of trific upstream degrader.");
   fTrifDMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
   delete fDetDir;
   delete fAppDir;
   delete fWorldDir;
   delete fWorldMaterialCmd;
   delete fWorldDimensionsCmd;
   delete fWorldVisCmd;
   delete fWorldMagneticFieldCmd;
   delete fWorldStepLimit;
   delete fDetSysDir;
   delete fUpdateCmd;
   delete fGenericTargetCmd;
   delete fGenericTargetDimensionsCmd;
   delete fGenericTargetPositionCmd;
   delete fLayeredTargetAddCmd;

   delete fFieldBoxMaterialCmd;
   delete fFieldBoxDimensionsCmd;
   delete fFieldBoxPositionCmd;
   delete fFieldBoxMagneticFieldCmd;
   delete fTabMagneticFieldCmd;
   delete fAddBoxMatCmd;
   delete fAddBoxThicknessCmd;
   delete fAddBoxInnerDimensionsCmd;
   delete fAddBoxColourCmd;
   delete fAddBoxCmd;
   delete fAddGridMatCmd;
   delete fAddGridSizeCmd;
   delete fAddGridDimensionsCmd;
   delete fAddGridColourCmd;
   delete fAddGridPosOffsetCmd;
   delete fAddGridCmd;
   delete fAddApparatusSpiceTargetChamberCmd;
   delete fAddDetectionSystemGammaTrackingCmd;
   delete fAddApparatus8piVacuumChamberCmd;
   delete fAddApparatus8piVacuumChamberAuxMatShellCmd;
   delete fAddApparatusGriffinStructureCmd;

   delete fAddDetectionSystemSodiumIodideCmd;
   delete fAddDetectionSystemLanthanumBromideCmd;
   delete fAddDetectionSystemAncillaryBGOCmd;

   delete fAddDetectionSystem8piCmd;
   delete fAddDetectionSystem8piDetectorCmd;
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   delete fAddDetectionSystemDescantCmd;
   delete fAddDetectionSystemDescantAuxPortsCmd;
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   delete fAddApparatusDescantStructureCmd;

   delete fAddDetectionSystemTestcanCmd;

   delete fSetDetectionSystemDescantColorCmd;
   delete fSetDetectionSystemDescantRotationCmd;
   delete fAddDetectionSystemDescantCartCmd;
   delete fAddDetectionSystemDescantSpherCmd;

   delete fAddDetectionSystemSceptarCmd;
   delete fAddDetectionSystemSpiceCmd;
   delete fAddDetectionSystemTrificCmd;
   delete fAddDetectionSystemPacesCmd;

   delete fAddDetectionSystemGriffinForwardCmd;
   delete fAddDetectionSystemGriffinForwardDetectorCmd;
   delete fAddDetectionSystemGriffinBackCmd;
   delete fAddDetectionSystemGriffinBackDetectorCmd;
   delete fAddDetectionSystemGriffinCustomDetectorCmd;
   delete fAddDetectionSystemGriffinCustomCmd;
   delete fAddDetectionSystemGriffinHevimetCmd;

   delete fAddDetectionSystemGriffinShieldSelectCmd;
   delete fAddDetectionSystemGriffinSetRadialDistanceCmd;
   delete fAddDetectionSystemGriffinSetExtensionSuppLocationCmd;
   delete fAddDetectionSystemGriffinSetPositionCmd;
   delete fAddDetectionSystemGriffinSetDeadLayerCmd;

   delete fUseTIGRESSPositionsCmd;

   delete fTrifWinCmd;
   delete fTrifDegCmd;
   delete fTrifAluCmd;
   delete fTrifFlatCmd;
   delete fTrifDMatCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
   if(command == fWorldMaterialCmd) {
      fDetector->SetWorldMaterial(newValue);
      return;
   }
   if(command == fWorldDimensionsCmd) {
      fDetector->SetWorldDimensions(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fWorldVisCmd) {
      fDetector->SetWorldVis(G4UIcmdWithABool::GetNewBoolValue(newValue));
   }
   if(command == fUpdateCmd) {
      fDetector->UpdateGeometry();
   }
   if(command == fWorldStepLimit) {
      fDetector->SetWorldStepLimit(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }
   if(command == fGenericTargetCmd) {
      fDetector->SetGenericTargetMaterial(newValue);
   }
   if(command == fGenericTargetDimensionsCmd) {
      fDetector->SetGenericTargetDimensions(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fGenericTargetPositionCmd) {
      fDetector->SetGenericTargetPosition(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fLayeredTargetAddCmd) {
      G4String           Material;
      G4double           areal = 0.;
      std::istringstream is(newValue);   ///string
      is >> Material >> areal;
      fDetector->LayeredTargetAdd(Material, areal);
   }
   if(command == fTabMagneticFieldCmd) {
      G4String           PathAndTableName;
      G4double           z_offset   = 0.;
      G4double           z_rotation = 0.;
      std::istringstream is(newValue);   ///string
      is >> PathAndTableName >> z_offset >> z_rotation;
      fDetector->SetTabMagneticField(PathAndTableName, z_offset, z_rotation);   // z in mm, angle in degree
   }
   if(command == fAddGridMatCmd) {
      fDetector->SetGridMat(newValue);
   }
   if(command == fAddGridSizeCmd) {
      fDetector->SetGridSize(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }
   if(command == fAddGridDimensionsCmd) {
      fDetector->SetGridDimensions(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddGridColourCmd) {
      fDetector->SetGridColour(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddGridPosOffsetCmd) {
      fDetector->SetGridPosOffset(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddGridCmd) {
      fDetector->AddGrid();
   }
   if(command == fAddApparatusSpiceTargetChamberCmd) {
      fDetector->AddApparatusSpiceTargetChamber(newValue);
   }
   if(command == fAddApparatus8piVacuumChamberCmd) {
      fDetector->AddApparatus8piVacuumChamber();
   }
   if(command == fAddApparatus8piVacuumChamberAuxMatShellCmd) {
      fDetector->AddApparatus8piVacuumChamberAuxMatShell(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }
   if(command == fAddApparatusGriffinStructureCmd) {
      fDetector->AddApparatusGriffinStructure(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemLanthanumBromideCmd) {
      fDetector->AddDetectionSystemLanthanumBromide(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddDetectionSystemAncillaryBGOCmd) {
      fDetector->AddDetectionSystemAncillaryBGO(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddDetectionSystem8piCmd) {
      fDetector->AddDetectionSystem8pi(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystem8piDetectorCmd) {
      fDetector->AddDetectionSystem8piDetector(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if(command == fAddDetectionSystemDescantCmd) {
      fDetector->AddDetectionSystemDescant(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemDescantAuxPortsCmd) {
      fDetector->AddDetectionSystemDescantAuxPorts(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if(command == fAddApparatusDescantStructureCmd) {
      fDetector->AddApparatusDescantStructure();
   }

   if(command == fSetDetectionSystemDescantColorCmd) {
      fDetector->SetDetectionSystemDescantColor(newValue);
   }
   if(command == fSetDetectionSystemDescantRotationCmd) {
      fDetector->SetDetectionSystemDescantRotation(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddDetectionSystemDescantCartCmd) {
      fDetector->AddDetectionSystemDescantCart(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddDetectionSystemDescantSpherCmd) {
      fDetector->AddDetectionSystemDescantSpher(G4UIcmdWith3VectorAndUnit::GetNew3VectorRawValue(newValue), G4UIcmdWithADoubleAndUnit::GetNewUnitValue(newValue));
   }
   if(command == fAddDetectionSystemTestcanCmd) {
      fDetector->AddDetectionSystemTestcan(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }

   if(command == fAddDetectionSystemSceptarCmd) {
      fDetector->AddDetectionSystemSceptar(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinForwardCmd) {
      fDetector->AddDetectionSystemGriffinForward(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinForwardDetectorCmd) {
      fDetector->AddDetectionSystemGriffinForwardDetector(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinBackCmd) {
      fDetector->AddDetectionSystemGriffinBack(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinBackDetectorCmd) {
      fDetector->AddDetectionSystemGriffinBackDetector(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinCustomDetectorCmd) {
      fDetector->AddDetectionSystemGriffinCustomDetector(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinCustomCmd) {
      fDetector->AddDetectionSystemGriffinCustom(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinHevimetCmd) {
      fDetector->AddDetectionSystemGriffinHevimet(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinShieldSelectCmd) {
      fDetector->AddDetectionSystemGriffinShieldSelect(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinSetRadialDistanceCmd) {
      fDetector->AddDetectionSystemGriffinSetRadialDistance(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinSetExtensionSuppLocationCmd) {
      fDetector->AddDetectionSystemGriffinSetExtensionSuppLocation(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinSetPositionCmd) {
      fDetector->AddDetectionSystemGriffinSetPosition(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddDetectionSystemGriffinSetDeadLayerCmd) {
      fDetector->AddDetectionSystemGriffinSetDeadLayer(G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
   }
   if(command == fAddDetectionSystemSpiceCmd) {
      fDetector->AddDetectionSystemSpice();
   }
   if(command == fAddDetectionSystemTrificCmd) {
      //Done as a string because Torr isnt a default unit
      G4double           torr = 0.;
      std::istringstream is(newValue);   ///string
      is >> torr;
      if(torr > 1) { fDetector->AddDetectionSystemTrific(torr);
		} else { fDetector->AddDetectionSystemTrific(760.0); }
   }
   if(command == fAddDetectionSystemPacesCmd) {
      fDetector->AddDetectionSystemPaces(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
   }
   if(command == fUseTIGRESSPositionsCmd) {
      fDetector->UseTIGRESSPositions(G4UIcmdWithABool::GetNewBoolValue(newValue));
   }
   if(command == fTrifWinCmd) {
      fDetector->SetTrifWindowThickness(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }
   if(command == fTrifDegCmd) {
      fDetector->SetTrifDegraderThickness(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }
   if(command == fTrifAluCmd) {
      fDetector->SetTrifAluminised(G4UIcmdWithABool::GetNewBoolValue(newValue));
   }
   if(command == fTrifDMatCmd) {
      fDetector->SetTrifDegraderMat(newValue);
   }
   if(command == fTrifFlatCmd) {
      fDetector->SetTrifFlatWindow(G4UIcmdWithABool::GetNewBoolValue(newValue));
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
