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
// $Id: PhysicsListMessenger.cc 68030 2013-03-13 13:51:27Z gcosmo $
//
/// \file radioactivedecay/rdecay02/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
   : fPPhysicsList(pPhys),
     fPhysDir(new G4UIdirectory("/DetSys/phys/")),
     fGammaCutCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/phys/setGCut", this)),
     fElectCutCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/phys/setECut", this)),
     fProtoCutCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/phys/setPCut", this)),
     fAllCutCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/phys/setCuts", this)),
     fECutCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/phys/TargetCuts", this)),
     fMCutCmd(new G4UIcmdWithADoubleAndUnit("/DetSys/phys/DetectorCuts", this)),
     fPListCmd(new G4UIcmdWithAString("/DetSys/phys/SelectPhysics", this)),
     fConstructOpCmd(new G4UIcmdWithABool("/DetSys/phys/ConstructOpticalPhysics", this)),
     fSpiceStepperCmd(new G4UIcmdWithABool("/DetSys/phys/SpiceStepper", this))
{
   fPhysDir->SetGuidance("physics control.");

   fGammaCutCmd->SetGuidance("Set gamma cut.");
   fGammaCutCmd->SetParameterName("Gcut", false);
   fGammaCutCmd->SetUnitCategory("Length");
   fGammaCutCmd->SetRange("Gcut>0.0");
   fGammaCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fElectCutCmd->SetGuidance("Set electron cut.");
   fElectCutCmd->SetParameterName("Ecut", false);
   fElectCutCmd->SetUnitCategory("Length");
   fElectCutCmd->SetRange("Ecut>0.0");
   fElectCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fProtoCutCmd->SetGuidance("Set positron cut.");
   fProtoCutCmd->SetParameterName("Pcut", false);
   fProtoCutCmd->SetUnitCategory("Length");
   fProtoCutCmd->SetRange("Pcut>0.0");
   fProtoCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fAllCutCmd->SetGuidance("Set cut for all.");
   fAllCutCmd->SetParameterName("cut", false);
   fAllCutCmd->SetUnitCategory("Length");
   fAllCutCmd->SetRange("cut>0.0");
   fAllCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fECutCmd->SetGuidance("Set cuts for the target");
   fECutCmd->SetParameterName("Ecut", false);
   fECutCmd->SetUnitCategory("Length");
   fECutCmd->SetRange("Ecut>0.0");
   fECutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fMCutCmd->SetGuidance("Set cuts for the Detector");
   fMCutCmd->SetParameterName("Ecut", false);
   fMCutCmd->SetUnitCategory("Length");
   fMCutCmd->SetRange("Ecut>0.0");
   fMCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fPListCmd->SetGuidance("Select modula physics list.");
   fPListCmd->SetParameterName("PList", false);
   fPListCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fConstructOpCmd->SetGuidance("Choose to build optical physics models");
   fConstructOpCmd->SetDefaultValue(false);
   fConstructOpCmd->AvailableForStates(G4State_Idle);   // Joey had removed the pre init state, not sure why

   fSpiceStepperCmd->SetGuidance("Choose to invoke SPICE stepper");
   fSpiceStepperCmd->SetDefaultValue(false);
   fSpiceStepperCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
   delete fPhysDir;
   delete fGammaCutCmd;
   delete fElectCutCmd;
   delete fProtoCutCmd;
   delete fAllCutCmd;
   delete fPListCmd;
   delete fECutCmd;
   delete fMCutCmd;
   delete fConstructOpCmd;
   delete fSpiceStepperCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                       G4String     newValue)
{
   if(command == fGammaCutCmd) {
      fPPhysicsList->SetCutForGamma(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }

   if(command == fElectCutCmd) {
      fPPhysicsList->SetCutForElectron(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }

   if(command == fProtoCutCmd) {
      fPPhysicsList->SetCutForPositron(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }

   if(command == fAllCutCmd) {
      G4double cut = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue);
      fPPhysicsList->SetCutForGamma(cut);
      fPPhysicsList->SetCutForElectron(cut);
      fPPhysicsList->SetCutForPositron(cut);
   }

   if(command == fECutCmd) {
      fPPhysicsList->SetTargetCut(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }

   if(command == fMCutCmd) {
      fPPhysicsList->SetDetectorCut(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
   }

   if(command == fPListCmd) {
      fPPhysicsList->SelectPhysicsList(newValue);
   }

   if(command == fConstructOpCmd) {
      fPPhysicsList->ConstructOp(G4UIcmdWithABool::GetNewBoolValue(newValue));
   }

   if(command == fSpiceStepperCmd) {
      fPPhysicsList->SpiceStepper(G4UIcmdWithABool::GetNewBoolValue(newValue));
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
