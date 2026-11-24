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
/// \file analysis/shared/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
// $Id: EventAction.hh 68015 2013-03-13 13:27:27Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "HistoManager.hh"

#include "DetectorConstruction.hh"   // for DetectorProperties

class RunAction;
class HistoManager;

static const int MAXSTEPS    = 1000;
static const int MAXHITS     = 100;
static const int NUMSTEPVARS = 15;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction {
public:
   EventAction(RunAction*, HistoManager*);
   ~EventAction() = default;

   void BeginOfEventAction(const G4Event*) override;
   void EndOfEventAction(const G4Event*) override;

   G4int GetEventNumber() const { return fEvtNb; };

   void AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& targetZ, const G4double& kinEnergy);
   void AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& targetZ);

   G4bool SpiceTest();

private:
   RunAction*    fRunAction;
   HistoManager* fHistoManager;

   G4int fPrintModulo;
   G4int fEvtNb;

   void ClearVariables();

   //Applying a resolution to SPCIE energies if desired
	std::array<G4double, 10000> fAmp;   //Bin amp (effective y)
   std::array<G4double, 10000> fAmpx;
   G4double fAmpTot;

   // Tracking info
	std::array<std::array<G4int, MAXHITS>, NUMSTEPVARS>              fHitTrackerI;
	std::array<std::array<G4double, MAXHITS>, NUMSTEPVARS>           fHitTrackerD;
   G4int              fNumberOfHits;
	std::array<DetectorProperties, MAXHITS> fProperties;

	std::array<std::array<G4int, MAXHITS>, NUMSTEPVARS>    fStepTrackerI;
	std::array<std::array<G4double, MAXHITS>, NUMSTEPVARS> fStepTrackerD;
   G4int    fNumberOfSteps;

   // Energy deposit in detection systems
	std::array<std::array<G4double, 12>, 10> fSpiceEnergyDet;
	std::array<std::array<G4double, 12>, 10> fSpiceTrackDet;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
