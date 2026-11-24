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
/// \file analysis/shared/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.hh 77256 2013-11-22 10:10:23Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include <unordered_map>

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
//class DetectionSystemGammaTracking;
class DetectionSystemGriffin;
class DetectionSystem8pi;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DetectionSystemDescant;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ApparatusDescantStructure;
class ApparatusLayeredTarget;
class DetectionSystemDescant;

class DetectionSystemSceptar;
class DetectionSystemSpice;
class DetectionSystemTrific;
class DetectionSystemPaces;
class DetectionSystemSodiumIodide;
class DetectionSystemLanthanumBromide;
class DetectionSystemBox;
class DetectionSystemAncillaryBGO;

//class MagneticField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct DetectorProperties {
   G4int systemID;
   G4int detectorNumber;
   G4int crystalNumber;
   void  Clear()
   {
      systemID       = 0;
      detectorNumber = 0;
      crystalNumber  = 0;
   }
};

bool operator==(const DetectorProperties& lhs, const DetectorProperties& rhs);

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
   DetectorConstruction();
   ~DetectorConstruction();

   G4int GriffinDetectorsMap(G4int i) const
   {
      if(i < 0 || i > 15) { return -1; }
      return fGriffinDetectorsMap[i];
   }

   void PassEfficiencyPosition(G4ThreeVector num) { fDetEffPosition = std::move(num); }

   void SetWorldMaterial(G4String);
   void SetWorldDimensions(const G4ThreeVector&);
   void SetWorldVis(G4bool);
   void SetWorldStepLimit(G4double);

   //Generic Target
   void SetGenericTargetMaterial(G4String);
   void SetGenericTargetDimensions(const G4ThreeVector&);
   void SetGenericTargetPosition(const G4ThreeVector&);
   void SetGenericTarget();

   G4double LayeredTargetLayerStart(int);

   void LayeredTargetAdd(const G4String&, G4double);

   void SetTabMagneticField(const G4String&, G4double, G4double);
   // Grid Functions
   void SetGridMat(G4String input) { fGridMat = std::move(input); };
   void SetGridSize(G4double input) { fGridSize = input; };
   void SetGridDimensions(G4ThreeVector input) { fGridDimensions = std::move(input); };
   void SetGridColour(G4ThreeVector input) { fGridColour = std::move(input); };
   void SetGridPosOffset(G4ThreeVector input) { fGridOffset = std::move(input); };
   void AddGrid();
   void AddApparatusSpiceTargetChamber(G4String);
   void AddApparatus8piVacuumChamber();
   void AddApparatus8piVacuumChamberAuxMatShell(G4double thickness);
   void AddApparatusGriffinStructure(G4int selector);

	// Trific setter functions
   void SetTrifWindowThickness(G4double input) { fTrifWindowThickness = input; }
   void SetTrifDegraderThickness(G4double input) { fTrifDegraderThickness = input; }
   void SetTrifDegraderMat(G4String input) { fTrifDegraderMat = std::move(input); }
   void SetTrifAluminised(G4bool input) { fTrifAluminised = input; }
   void SetTrifFlatWindow(G4bool input) { fTrifFlatWindow = input; }

   G4double GetWorldSizeX() const { return fWorldSizeX; };
   G4double GetWorldSizeY() const { return fWorldSizeY; };
   G4double GetWorldSizeZ() const { return fWorldSizeZ; };

   const G4VPhysicalVolume* GetphysiWorld() const { return fPhysiWorld; };

   G4VPhysicalVolume* Construct() override;

   void UpdateGeometry();

   void AddDetectionSystemSodiumIodide(G4int ndet);
   void AddDetectionSystemLanthanumBromide(const G4ThreeVector& input);
   void AddDetectionSystemAncillaryBGO(const G4ThreeVector& input);

   void AddDetectionSystem8pi(G4int ndet);

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   void AddDetectionSystemDescant(G4int ndet);
   void AddDetectionSystemDescantAuxPorts(const G4ThreeVector& input);

   void SetDetectionSystemDescantRotation(const G4ThreeVector& input);
   void SetDetectionSystemDescantColor(G4String input);
   void AddDetectionSystemDescantCart(const G4ThreeVector& input);
   void AddDetectionSystemDescantSpher(const G4ThreeVector& input, G4double unit);
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   void AddApparatusDescantStructure();

   void AddDetectionSystemTestcan(const G4ThreeVector& input);

   void AddDetectionSystem8piDetector(G4int ndet);
   void AddDetectionSystemGriffinForward(G4int ndet);
   void AddDetectionSystemGriffinForwardDetector(G4int ndet);
   void AddDetectionSystemGriffinBack(G4int ndet);
   void AddDetectionSystemGriffinBackDetector(G4int ndet);
   //void AddDetectionSystemGriffinPositionConfig(G4ThreeVector input);
   void AddDetectionSystemGriffinHevimet(G4int input);
   void AddDetectionSystemGriffinCustom(G4int ndet);
   void AddDetectionSystemGriffinCustomDetector(G4int ndet = 0);
   void AddDetectionSystemGriffinShieldSelect(G4int ShieldSelect);
   void AddDetectionSystemGriffinSetRadialDistance(G4double detectorDist);
   void AddDetectionSystemGriffinSetExtensionSuppLocation(G4int detectorPos);
   void AddDetectionSystemGriffinSetPosition(const G4ThreeVector& params);
   void AddDetectionSystemGriffinSetDeadLayer(const G4ThreeVector& params);

   void AddDetectionSystemSceptar(G4int ndet);
   void AddDetectionSystemPaces(G4int ndet);
   void AddDetectionSystemSpice();
   void AddDetectionSystemTrific(G4double);

   G4double GetLanthanumBromideCrystalRadius();
   G4double GetLanthanumBromideCrystalLength();
   G4double GetLanthanumBromideR();
   G4double GetLanthanumBromideTheta(G4int i);
   G4double GetLanthanumBromidePhi(G4int i);
   G4double GetLanthanumBromideYaw(G4int i);
   G4double GetLanthanumBromidePitch(G4int i);
   G4double GetLanthanumBromideRoll(G4int i);
   G4double GetLanthanumBromideCrystalRadialPosition();

   G4bool GridCell() const { return fGridCell; }
   G4bool Griffin() const { return fGriffin; }
   G4bool LaBr() const { return fLaBr; }
   G4bool AncBgo() const { return fAncBgo; }
   G4bool NaI() const { return fNaI; }
   G4bool Sceptar() const { return fSceptar; }
   G4bool EightPi() const { return fEightPi; }
   G4bool Spice() const { return fSpice; }
   G4bool Paces() const { return fPaces; }
   G4bool Descant() const { return fDescant; }
   G4bool Testcan() const { return fTestcan; }

   void UseTIGRESSPositions(G4bool input) { fUseTigressPositions = input; }

   bool               HasProperties(G4VPhysicalVolume* vol) const { return fPropertiesMap.find(vol) != fPropertiesMap.end(); }
   DetectorProperties GetProperties(G4VPhysicalVolume* vol) const { return fPropertiesMap.at(vol); }
   void               SetProperties();
   void               SetPropertiesRecursive(G4LogicalVolume* vol);

   void Print();
   void PrintRecursive(G4LogicalVolume* vol);

private:
   bool               CheckVolumeName(const G4String& volumeName);
   DetectorProperties ParseVolumeName(const G4String& volumeName);

   G4int fGriffinDetectorsMapIndex;
	std::array<G4int, 16> fGriffinDetectorsMap;

   G4double fWorldSizeX;
   G4double fWorldSizeY;
   G4double fWorldSizeZ;
   G4bool   fWorldVis;
   G4bool   fBuiltDetectors;
   G4double fGriffinFwdBackPosition;
   G4int    fDetectorShieldSelect;
   G4double fDetectorRadialDistance;
   G4int    fExtensionSuppressorLocation;
   G4int    fCustomDetectorNumber;
   G4int    fCustomDetectorPosition;
   G4int    fCustomDetectorVal;
   G4int    fHevimetSelector;
   G4bool   fUseTigressPositions;

   // Box
   G4String      fBoxMat;
   G4double      fBoxThickness;
   G4ThreeVector fBoxInnerDimensions;
   G4ThreeVector fBoxColour;

   G4Box*             fSolidWorld{nullptr};   //pointer to the solid World
   G4LogicalVolume*   fLogicWorld{nullptr};   //pointer to the logical World
   G4VPhysicalVolume* fPhysiWorld{nullptr};   //pointer to the physical World

   // Grid
   G4String      fGridMat;
   G4double      fGridSize;
   G4ThreeVector fGridDimensions;
   G4ThreeVector fGridColour;
   G4ThreeVector fGridOffset;

   void DefineSuppressedParameters();
   void DefineMaterials();

   std::array<std::array<G4double, 5>, 20> fCoords;
   G4bool                                  fSetGenericTargetMaterial;
   G4bool                                  fSetGenericTargetDimensions;
   G4bool                                  fSetGenericTargetPosition;
   G4String                                fGenericTargetMaterial;
   G4ThreeVector                           fGenericTargetDimensions;
   G4ThreeVector                           fGenericTargetPosition;

   G4bool        fSetFieldBoxMaterial;
   G4bool        fSetFieldBoxDimensions;
   G4bool        fSetFieldBoxPosition;
   G4bool        fSetFieldBoxMagneticField;
   G4String      fFieldBoxMaterial;
   G4ThreeVector fFieldBoxDimensions;
   G4ThreeVector fFieldBoxPosition;
   G4ThreeVector fFieldBoxMagneticField;

   G4String fMatWorldName;

   ApparatusLayeredTarget* fApparatusLayeredTarget;
   DetectorMessenger*      fDetectorMessenger;

   std::array<std::array<G4double, 4>, 16> fGriffinDeadLayer;

   G4ThreeVector fDescantRotation;
   G4String      fDescantColor;

   G4ThreeVector fDetEffPosition;

   //booleans which control which histograms are created (these are set by the detector construction)
   G4bool fGridCell;
   G4bool fGriffin;
   G4bool fLaBr;
   G4bool fAncBgo;
   G4bool fNaI;
   G4bool fSceptar;
   G4bool fEightPi;
   G4bool fDescant;
   G4bool fTestcan;
   G4bool fSpice;
   G4bool fPaces;

   G4double fTrifWindowThickness;
   G4double fTrifDegraderThickness;
   G4String fTrifDegraderMat;
   G4bool   fTrifAluminised;
   G4bool   fTrifFlatWindow;

   //unordered maps which hold properties of the physical volumes created
   std::unordered_map<G4VPhysicalVolume*, DetectorProperties> fPropertiesMap;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
