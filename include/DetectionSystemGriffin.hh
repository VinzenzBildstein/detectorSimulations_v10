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

#ifndef DETECTIONSYSTEMGRIFFIN_HH
#define DETECTIONSYSTEMGRIFFIN_HH

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4AssemblyVolume.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Colour.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectionSystemGriffin {
public:
   DetectionSystemGriffin(G4int sel, G4int suppSwitch, G4double detRad, G4int hevimetSel);
   ~DetectionSystemGriffin();

   void Build();
   // For detector specific dead layers
   void BuildDeadLayerSpecificCrystal(G4int det);

   void     BuildEverythingButCrystals(G4int det);
   G4double GetCrystalDistanceFromOrigin() const { return fCrystalDistFromOrigin; }

   G4double TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
   G4double TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
   G4double TransZ(G4double x, G4double z, G4double theta);

   // For detector specific dead layers
	G4int PlaceDeadLayerSpecificCrystal(G4LogicalVolume* expHallLog, G4int detectorNumber, G4int positionNumber, G4bool posTigress);
   G4int PlaceEverythingButCrystals(G4LogicalVolume* expHallLog, G4int detectorNumber, G4int positionNumber, G4bool posTigress);

   void SetDeadLayer(G4int detNum, G4int cryNum, G4double deadLayer) { fGriffinDeadLayers[detNum][cryNum] = deadLayer; }

private:
   G4String fSdName0;
   G4String fSdName1;
   G4String fSdName2;
   G4String fSdName3;
   G4String fSdName4;
   G4String fSdName5;
   G4String fColNameGe;
   G4String fColNameLeftCasing;
   G4String fColNameRightCasing;
   G4String fColNameLeftExtension;
   G4String fColNameRightExtension;
   G4String fColNameBackPlug;

   G4bool fIncludeExtensionSuppressors;
   G4bool fIncludeSideSuppressors;
   G4bool fIncludeBackSuppressors;

   G4String fBackSuppressorMaterial;
   G4String fBGOMaterial;

   G4RotationMatrix* fRotateNull;
   G4ThreeVector     fMoveNull;

   G4double fCutClearance;
   G4double fExtraCutLength;

   G4bool fSurfCheck;

	std::array<std::array<double, 5>, 20> fCoords;

   //Jun 21, 2005, Epapr7.80: modification to the prototype suppressor shield,
   //since the hevimet collimators they provided were too big, and therefore the
   //extensions do not come far enough forward
   G4double fExtensionAccidentalBackShift;

   G4String fHevimetChoice;

   //Settings for the suppressor shield
   G4double fThicknessDouble;
   G4double fSideThicknessDouble;
   G4double fExtensionThicknessDouble;

   G4double fRadialDistance;
   G4int    fNumberOfDetectors;
   G4double fOriginToCrystalDistance;
   G4bool   fDeadLayerIncludeFlag;
   G4double fInnerDeadLayerThickness;
   G4double fOuterDeadLayerThickness;

   //Cold Finger
   G4double fColdFingerOuterShellRadius;
   G4double fColdFingerShellThickness;
   G4double fColdFingerShellLength;
   G4double fRearPlateThickness;
   G4double fColdFingerEndPlateThickness;
   G4double fColdFingerRadius;
   G4double fColdFingerSpace;
   G4double fColdFingerLength;
   G4double fCoolantLength;
   G4double fCoolantRadius;
   G4double fCoolantThickness;

   //New cooling structures (Jan 2005)
   G4double fExtraBlockThickness;
   G4double fExtraBlockDistanceFromBackPlate;
   G4double fExtraBlockInnerDiameter;
   G4double fTrianglePostsDistanceFromCrystals;
   G4double fTrianglePostStartingDepth;
   G4double fFetAirHoleRadius;
   G4double fCoolingSideBlockThickness;
   G4double fCoolingSideBlockWidth;
   G4double fCoolingSideBlockHorizontalDepth;
   G4double fStructureMatColdFingerThickness;
   G4double fStructureMatColdFingerRadius;
   G4double fCoolingBarThickness;
   G4double fCoolingBarWidth;

   G4int fSuppressorDesign;
   G4int fSuppressorPositionSelector;
   G4int fHevimetSelector;

   G4double fForwardInnerRadius;
   G4double fBackInnerRadius;

   G4double fSuppressorForwardRadius;
   G4double fSuppressorBackRadius;

   // For the optimization of the depth segmentation
   G4double fDepthSegmentationAdjustment;

   //the germanium detector's variables (for one "leaf" of the clover)
   G4double fGermaniumOuterRadius;
   G4double fGermaniumHoleRadius;
   G4double fGermaniumWidth;
   G4double fGermaniumLength;
   G4double fGermaniumHoleDistFromFace;
   G4double fGermaniumDistFromCanFace;
   G4double fGermaniumBentLength;
   G4double fGermaniumShift;   //the amount by which the two sides adjacent to the other
   //germanium crystals are cut closer to the center

   G4double fGermaniumSeparation;   //the space between the quarter detectors
   G4double fInterCrystalElectrodeMatThickness;
   G4double fElectrodeMatStartingDepth;
   G4double fGermaniumCornerConeEndLength;   //the ending length of the cones
   //at the corners of each leaf

   //basic germanium detector's values
   G4double fDetectorBlockLength;
   G4double fDetectorBlockHeight;
   G4double fDetectorBlockTrapezoidalHeight;

   G4double fDetectorTotalWidth;
   G4double fDetectorTotalLength;
   G4double fBentEndAngle;
   G4double fBentEndLength;
   G4double fCanFaceThickness;
   G4double fCanSideThickness;

   G4double fHevimetTipThickness;
   G4double fHevimetTipAngle;

   //Values for the BGO
   G4double fSuppressorCutExtension;
   G4double fSuppressorShellThickness;

   G4double fBackBGOThickness;
   G4double fBGOChoppedTip;
   G4double fBGOMovableSpace;
   G4double fSideSuppressorBackOffset;
   G4double fSideBGOThickness;
   G4double fBGOCanSeperation;
   G4double fSideBGOLength;
   G4double fSideSuppressorLength;
   G4double fBGOTrapLength;
   G4double fSuppressorExtensionThickness;
   G4double fSuppressorExtensionLength;
   G4double fSuppressorExtensionAngle;

   G4double fSuppressorExtensionLengthDet;

   //Values for the HeavyMet
   G4double fHeavyMetThickness;
   G4double fHeavyMetInsideAngle;

   G4double fAirBoxFrontWidth;
   G4double fAirBoxFrontLength;
   G4double fAirBoxBackLength;

   G4double fAirBoxBackLengthDet;
   G4double fAirBoxFrontLengthDet;
   G4double fAirBoxFrontWidthDet;

   G4double fShift;
   G4double fSuppShift;

   G4int fCopyNumber;
   G4int fCopyNumberTwo;
   G4int fGermaniumCopyNumber;
   G4int fLeftSuppressorSideCopyNumber;
   G4int fRightSuppressorSideCopyNumber;
   G4int fLeftSuppressorExtensionCopyNumber;
   G4int fRightSuppressorExtensionCopyNumber;
   G4int fBackSuppressorCopyNumber;

   G4double fRhombiDiameter;
   G4double fNewRhombiRadius;

   G4double fNewRhombiRadiusDet;

   G4double fDetectorPositionShift;
   G4double fAppliedBackShift;

   G4int fGermaniumSelector;
   G4int fCanSelector;
   G4int fBGOSelector;
   G4int fColdFingerSelector;

   //Redacted parameters/////////////
   G4double fDetectorPlacementCxn;
   G4double fTrianglePostDim;
   G4double fSuppressorExtRightX;
   G4double fSuppressorExtRightY;
   G4double fSuppressorExtRightZ;
   G4double fSuppressorExtLeftX;
   G4double fSuppressorExtLeftY;
   G4double fSuppressorExtLeftZ;
   G4double fWedgeDim;
   G4double fQuarterDetectorCxn;
   G4double fQuarterDetectorCxnB;
   G4String fElectrodeMaterial;
   G4String fStructureMaterial;

   // Assembly volumes
   G4AssemblyVolume* fAssembly{nullptr};
   G4AssemblyVolume* fGermaniumAssembly{nullptr};
   G4AssemblyVolume* fLeftSuppressorCasingAssembly{nullptr};
   G4AssemblyVolume* fRightSuppressorCasingAssembly{nullptr};
   G4AssemblyVolume* fLeftSuppressorExtensionAssembly{nullptr};
   G4AssemblyVolume* fRightSuppressorExtensionAssembly{nullptr};
   G4AssemblyVolume* fSuppressorBackAssembly{nullptr};
   // For detector specific dead layers
	std::array<G4AssemblyVolume*, 4> fAssemblyCry;
   std::array<G4AssemblyVolume*, 4> fGermaniumAssemblyCry;
   std::array<G4AssemblyVolume*, 4> fLeftSuppressorCasingAssemblyCry;
   std::array<G4AssemblyVolume*, 4> fRightSuppressorCasingAssemblyCry;
   std::array<G4AssemblyVolume*, 4> fLeftSuppressorExtensionAssemblyCry;
   std::array<G4AssemblyVolume*, 4> fRightSuppressorExtensionAssemblyCry;
   std::array<G4AssemblyVolume*, 4> fSuppressorBackAssemblyCry;
   G4AssemblyVolume* fSuppressorShellAssembly{nullptr};
   G4AssemblyVolume* fBackAndSideSuppressorShellAssembly{nullptr};
   G4AssemblyVolume* fExtensionSuppressorShellAssembly{nullptr};
   G4AssemblyVolume* fHevimetAssembly{nullptr};

   // Logical volumes
   G4LogicalVolume* fAirBoxLog{nullptr};

   // methods to construct all of the components of the detector
   void ConstructNewSuppressorCasingWithShells(G4int det);
   void BuildelectrodeMatElectrodes();
   void ConstructComplexDetectorBlockWithDeadLayer();

   // For detector specific dead layers
   void ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer(G4int det, G4int cry);
   void ConstructNewSuppressorCasingDetectorSpecificDeadLayer(G4int det, G4int cry);
   void ConstructDetector();
   void ConstructComplexDetectorBlock();
   void ConstructColdFinger();
   void ConstructNewHeavyMet();

   //LogicalVolumes used in ConstructBasicDetectorBlock
   G4LogicalVolume* fGermaniumBlockLog{nullptr};

   //Logical Volumes used in ConstructComplexDetectorBlock:
   G4LogicalVolume* fGermaniumBlock1Log{nullptr};
   G4LogicalVolume* fGermaniumHoleLog{nullptr};
   G4LogicalVolume* fInnerDeadLayerLog{nullptr};
   G4LogicalVolume* fInnerDeadLayerCapLog{nullptr};
   G4LogicalVolume* fOuterDeadLayerLog{nullptr};

   G4LogicalVolume* fInterCrystalElectrodeMatBackLog{nullptr};
   G4LogicalVolume* fInterCrystalElectrodeMatFrontLog{nullptr};

   //Logical Volumes used in ConstructBGOCasing:
   G4LogicalVolume* fBackBGOLog{nullptr};
   G4LogicalVolume* fBGOCasingLog{nullptr};

   //Logical Volumes used in ConstructNewSuppressorCasing:
   G4LogicalVolume* fBackQuarterSuppressorShellLog{nullptr};
   G4LogicalVolume* fRightSuppressorShellLog{nullptr};
   G4LogicalVolume* fLeftSuppressorShellLog{nullptr};
   G4LogicalVolume* fRightSuppressorShellExtensionLog{nullptr};
   G4LogicalVolume* fLeftSuppressorShellExtensionLog{nullptr};

   G4LogicalVolume* fCapForRightSuppressorLog{nullptr};

   G4LogicalVolume* fBackQuarterSuppressorLog{nullptr};
   G4LogicalVolume* fRightSuppressorLog{nullptr};
   G4LogicalVolume* fLeftSuppressorLog{nullptr};
   G4LogicalVolume* fRightSuppressorExtensionLog{nullptr};
   G4LogicalVolume* fLeftSuppressorExtensionLog{nullptr};

   //Logical Volumes used in ConstructDetector:
   G4LogicalVolume* fFrontFaceLog{nullptr};
   G4LogicalVolume* fRightBentPieceLog{nullptr};
   G4LogicalVolume* fLeftBentPieceLog{nullptr};
   G4LogicalVolume* fTopBentPieceLog{nullptr};
   G4LogicalVolume* fBottomBentPieceLog{nullptr};
   G4LogicalVolume* fRightWedgeLog{nullptr};
   G4LogicalVolume* fLeftWedgeLog{nullptr};
   G4LogicalVolume* fTopWedgeLog{nullptr};
   G4LogicalVolume* fBottomWedgeLog{nullptr};
   G4LogicalVolume* fUpperRightConeLog{nullptr};
   G4LogicalVolume* fLowerRightConeLog{nullptr};
   G4LogicalVolume* fUpperLeftConeLog{nullptr};
   G4LogicalVolume* fLowerLeftConeLog{nullptr};
   G4LogicalVolume* fUpperRightTubeLog{nullptr};
   G4LogicalVolume* fLowerRightTubeLog{nullptr};
   G4LogicalVolume* fUpperLeftTubeLog{nullptr};
   G4LogicalVolume* fLowerLeftTubeLog{nullptr};
   G4LogicalVolume* fRightSidePanelLog{nullptr};
   G4LogicalVolume* fLeftSidePanelLog{nullptr};
   G4LogicalVolume* fTopSidePanelLog{nullptr};
   G4LogicalVolume* fBottomSidePanelLog{nullptr};
   G4LogicalVolume* fRearPlateLog{nullptr};
   G4LogicalVolume* fFingerShellLog{nullptr};
   G4LogicalVolume* fTankLog{nullptr};
   G4LogicalVolume* fTankLid1Log{nullptr};
   G4LogicalVolume* fTankLid2Log{nullptr};
   G4LogicalVolume* fTankLiquidLog{nullptr};

   //Logical Volumes used in ConstructColdFinger:
   G4LogicalVolume* fEndPlateLog{nullptr};
   G4LogicalVolume* fFingerLog{nullptr};
   G4LogicalVolume* fExtraColdBlockLog{nullptr};
   G4LogicalVolume* fTrianglePostLog{nullptr};
   G4LogicalVolume* fFetAirHoleLog{nullptr};
   G4LogicalVolume* fCoolingBarLog{nullptr};
   G4LogicalVolume* fCoolingSideBlockLog{nullptr};
   G4LogicalVolume* fStructureMatColdFingerLog{nullptr};

   //Logical Volumes used in ConstructNewHeavyMet:
   G4LogicalVolume* fHevimetLog{nullptr};

   //internal methods for ConstructCan()
   G4Box*              SquareFrontFace() const;
   G4Trap*             CornerWedge() const;
   G4Para*             BentSidePiece() const;
   G4Box*              OtherBentSidePiece() const;
   G4Cons*             RoundedEndEdge() const;
   G4Tubs*             CornerTube() const;
   G4Box*              SidePanel() const;
   G4SubtractionSolid* RearPlate() const;
   G4Tubs*             ColdFingerShell() const;
   G4Tubs*             LiquidNitrogenTank() const;
   G4Tubs*             LiquidNitrogenTankLid() const;
   G4Tubs*             LiquidNitrogen() const;

   //internal methods for ConstructBasicDetectorBlock()
   G4Box* RectangularSegment() const;
   G4Trd* TrapezoidalSegment() const;

   //internal methods for ConstructComplexDetectorBlock()
   G4SubtractionSolid* QuarterDetector() const;
   // For detector specific dead layers
   G4SubtractionSolid* QuarterSpecificDeadLayerDetector(G4int det, G4int cry);

   //internal methods for ConstructComplexDetectorBlockWithPlastic()
   G4UnionSolid* InterCrystalelectrodeMatBack() const;
   G4UnionSolid* InterCrystalelectrodeMatFront() const;

   //internal methods for ConstructColdFinger()
   G4Tubs*             AirHole() const;
   G4Tubs*             AirHoleCut() const;
   G4Box*              EndPlate() const;
   G4Tubs*             Finger() const;
   G4SubtractionSolid* ExtraColdBlock();
   G4Trd*              TrianglePost() const;
   G4Box*              CoolingBar() const;
   G4Box*              CoolingSideBlock() const;
   G4Tubs*             StructureMatColdFinger() const;

   //internal methods for ConstructBGOCasing()
   G4SubtractionSolid* BackBGO();
   G4SubtractionSolid* BGOCasing();
   G4SubtractionSolid* FrontBGO();
   G4Trd*              SideBGO();

   //internal methods for ConstructNewSuppressorCasing()
   G4SubtractionSolid* BackSuppressorQuarter() const;
   G4SubtractionSolid* FrontSlantSuppressor(const G4String& sidePosition, G4bool choppingSuppressor) const;
   G4SubtractionSolid* SideSuppressorExtension(const G4String& sidePosition, G4bool choppingSuppressor) const;
   G4Trap*             SideSuppressorExtensionUncut() const;

   //internal methods for New SuppressorCasingWithShells
   G4SubtractionSolid* ShellForBackSuppressorQuarter() const;
   G4SubtractionSolid* ShellForFrontSlantSuppressor(const G4String& sidePosition);
   G4SubtractionSolid* ShellForSuppressorExtension(const G4String& sidePosition);

   //internal methods for ConstructNewHeavyMet()
   G4SubtractionSolid* NewHeavyMet() const;

   G4String fCrystalMaterial;
   G4String fCanMaterial;
   G4String fVacuumMaterial;
   G4double fCrystalLengthX;
   G4double fCrystalLengthY;
   G4double fCrystalLengthZ;
   G4double fCrystalInnerRadius;
   G4double fCrystalOuterRadius;
   G4double fCanThickness;
   G4double fCanInnerRadius;
   G4double fCanLidInnerRadius;
   G4double fCanLidOuterRadius;
   G4double fCanFrontLidThickness;
   G4double fCanBackLidThickness;
   G4double fCanFaceDistFromOrigin;
   G4double fCrystalDistFromCanFace;
   G4double fCrystalDistFromCanBack;
   G4double fCanLengthZ;
   G4double fCrystalDistFromOrigin;

   // For detector specific dead layers
	std::array<std::array<G4double, 4>, 16> fGriffinDeadLayers;
   std::array<G4Colour, 4> fGriffinCrystalColours;
   std::array<G4Colour, 4> fGriffinDeadLayerColours;

   // internal methods
   void BuildOneDetector(G4int det);
   //    void PlaceDetector(G4int detectorNumber);
};

#endif
