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
// $Id: PrimaryGeneratorMessenger.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PRIMARYGENERATORMESSENGER_HH
#define PRIMARYGENERATORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
class PrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger : public G4UImessenger {
public:
   explicit PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
   virtual ~PrimaryGeneratorMessenger();

   void SetNewValue(G4UIcommand*, G4String) override;

private:
   PrimaryGeneratorAction* fAction{nullptr};

   G4UIcmdWithAnInteger*      fNumberOfDecayingLaBrDetectorsCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fEfficiencyEnergyCmd{nullptr};
   G4UIcmdWith3Vector*        fEfficiencyDirectionCmd{nullptr};
   G4UIcmdWith3VectorAndUnit* fEfficiencyPositionCmd{nullptr};
   G4UIcmdWithAString*        fEfficiencyParticleCmd{nullptr};
   G4UIcmdWith3Vector*        fEfficiencyPolarizationCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fConeRadiusCmd{nullptr};   //SPICE cone radius
   G4UIcmdWithADoubleAndUnit* fConeZValueCmd{nullptr};   //SPICE cone with values
   G4UIcmdWithADoubleAndUnit* fConeRValueCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fConeAngleCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fConeMinAngleCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fBeamSpotSigmaCmd{nullptr};
   G4UIcmdWithAnInteger*      fBeamDistroCmd{nullptr};
   G4UIcmdWithAString*        fBeamFileCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fKentuckyEnergyCmd{nullptr};
   G4UIcmdWithAString*        fKentuckyReactionCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fMinimumPhiCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fMaximumPhiCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fMinimumThetaCmd{nullptr};
   G4UIcmdWithADoubleAndUnit* fMaximumThetaCmd{nullptr};
   G4UIcmdWithAnInteger*      fVerbosityCmd{nullptr};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
