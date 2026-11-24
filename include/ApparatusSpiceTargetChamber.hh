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

#ifndef APPARATUSSPICETARGETCHAMBER_HH
#define APPARATUSSPICETARGETCHAMBER_HH

#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4String.hh"

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

// define custom colours for visualisation
#define SN_COL 1.0, 1.0, 1.0
#define CU_COL 1.0, 1.0, 0.0
#define PEEK_COL 0.5, 0.5, 0.0
#define KAPTON_COL 0.2, 0.7, 0.1
#define PB_COL 0.6, 0.1, 0.1
#define NDFEB_COL 0.7, 0.3, 0.3
#define DELRIN_COL 0.0, 0.0, 1.0
#ifndef AL_COL
#define AL_COL 0.5, 0.5, 0.5
#endif

#define TARGET_CHAMBER_COL 0.15, 0.15, 0.15
#define ELECTROBOX_COL 0.15, 0.15, 0.15

class ApparatusSpiceTargetChamber {
public:
   explicit ApparatusSpiceTargetChamber(G4String);
   ~ApparatusSpiceTargetChamber();

   void Build(G4LogicalVolume*);

   static G4double TargetChamberDistanceFromTarget() { return fTargetChamberDistanceFromTarget; }

private:
   ////////////////////////////////////////////
   // Logical Volumes used in ApparatusSpiceTargetChamber
   ////////////////////////////////////////////
   G4LogicalVolume* fTargetChamberFrontRingLog{nullptr};
   G4LogicalVolume* fTargetChamberFrontConeLog{nullptr};
   G4LogicalVolume* fTargetChamberSphereLog{nullptr};
   G4LogicalVolume* fTargetChamberCylinderDownLog{nullptr};
   G4LogicalVolume* fTargetWheelLog{nullptr};
   G4LogicalVolume* fFrameLogical{nullptr};
   G4LogicalVolume* fCollimLogical{nullptr};
   G4LogicalVolume* fExtLogical{nullptr};
   G4LogicalVolume* fRodLogical{nullptr};
   G4LogicalVolume* fFirstGearLog{nullptr};
   G4LogicalVolume* fSecondGearLog{nullptr};
   G4LogicalVolume* fThirdGearLog{nullptr};
   G4LogicalVolume* fGearPlateOneLog{nullptr};
   G4LogicalVolume* fGearPlateTwoLog{nullptr};
   G4LogicalVolume* fGearStickLog{nullptr};
   G4LogicalVolume* fTargetMountPlateLog{nullptr};
   G4LogicalVolume* fTargetMountRodLog{nullptr};
   G4LogicalVolume* fBiasPlateLog{nullptr};
   G4LogicalVolume* fPhotonShieldLayerOneLog{nullptr};
   G4LogicalVolume* fPhotonShieldLayerTwoLog{nullptr};
   G4LogicalVolume* fPhotonShieldLayerThreeLog{nullptr};
   G4LogicalVolume* fPsTargetClampLog{nullptr};
   G4LogicalVolume* fPsDetectorClampLog{nullptr};
   G4LogicalVolume* fShieldBoltLog{nullptr};
   G4LogicalVolume* fDetectorBoltLog{nullptr};
   G4LogicalVolume* fMagnetLog{nullptr};
   G4LogicalVolume* fMagnetCoverLog{nullptr};
   G4LogicalVolume* fMclampChamberLog{nullptr};
   G4LogicalVolume* fElectroBoxLog{nullptr};
   G4LogicalVolume* fShieldCoverLog{nullptr};
   G4LogicalVolume* fColdFingerLog{nullptr};
   G4LogicalVolume* fS3CaseLogical{nullptr};
   G4LogicalVolume* fConicalCollimatorLog{nullptr};   //11/8
   G4LogicalVolume* fXRayInsertLog{nullptr};

   ////////////////////////////////////////////
   // Physical Volumes used in ApparatusSpiceTargetChamber
   ////////////////////////////////////////////
   G4VPhysicalVolume* fTargetChamberFrontRingPhys{nullptr};
   G4VPhysicalVolume* fTargetChamberFrontConePhys{nullptr};
   G4VPhysicalVolume* fTargetChamberSpherePhys{nullptr};
   G4VPhysicalVolume* fTargetChamberCylinderDownPhys{nullptr};
   G4VPhysicalVolume* fTargetWheelPhys{nullptr};
   G4VPhysicalVolume* fFramePhysical{nullptr};
   G4VPhysicalVolume* fCollimPhysical{nullptr};
   G4VPhysicalVolume* fExtPhysical{nullptr};
   G4VPhysicalVolume* fRodPhysical{nullptr};
   G4VPhysicalVolume* fFirstGearPhys{nullptr};
   G4VPhysicalVolume* fSecondGearPhys{nullptr};
   G4VPhysicalVolume* fThirdGearPhys{nullptr};
   G4VPhysicalVolume* fGearPlateOnePhys{nullptr};
   G4VPhysicalVolume* fGearPlateTwoPhys{nullptr};
   G4VPhysicalVolume* fGearStickPhys{nullptr};
   G4VPhysicalVolume* fTargetMountPlatePhys{nullptr};
   G4VPhysicalVolume* fTargetMountRodPhys{nullptr};
   G4VPhysicalVolume* fBiasPlatePhys{nullptr};
   G4VPhysicalVolume* fPhotonShieldLayerOnePhys{nullptr};
   G4VPhysicalVolume* fPhotonShieldLayerTwoPhys{nullptr};
   G4VPhysicalVolume* fPhotonShieldLayerThreePhys{nullptr};
   G4VPhysicalVolume* fPsTargetClampPhys{nullptr};
   G4VPhysicalVolume* fPsDetectorClampPhys{nullptr};
   G4VPhysicalVolume* fTargetBoltPhys{nullptr};
   G4VPhysicalVolume* fDetectorBoltPhys{nullptr};
   G4VPhysicalVolume* fMagnetPhys{nullptr};
   G4VPhysicalVolume* fMagnetCoverPhys{nullptr};
   G4VPhysicalVolume* fMclampChamberPhys{nullptr};
   G4VPhysicalVolume* fMclampShieldPhys{nullptr};
   G4VPhysicalVolume* fElectroBoxPhys{nullptr};
   G4VPhysicalVolume* fShieldCoverPhys{nullptr};
   G4VPhysicalVolume* fColdFingerPhys{nullptr};
   G4VPhysicalVolume* fS3CasePhysical{nullptr};
   G4VPhysicalVolume* fConicalCollimatorPhys{nullptr};   //11/8
   G4VPhysicalVolume* fXRayInsertPhys{nullptr};

   ////////////////////////////////////////////
   // Properties used in ApparatusSpiceTargetChamber
   ////////////////////////////////////////////

   //-------------------------
   // Materials:
   //-------------------------
   G4String fMagnetMaterial;                   //
   G4String fTargetChamberMaterial;            //
   G4String fDownstreamConeMaterial;           //
   G4String fMountLead;                        //doesn't exist
   G4String fPhotonShieldLayerOneMaterial;     //
   G4String fPhotonShieldLayerTwoMaterial;     //
   G4String fPhotonShieldLayerThreeMaterial;   //
   G4String fPsClampMaterial;                  //
   G4String fMagnetChamberClampMaterial;       //
   G4String fTargetWheelMaterial;              //
   G4String fGearStickMaterial;                //
   G4String fBiasPlateMaterial;                //
   G4String fTargetWheelGearMaterialA;         //
   G4String fTargetWheelGearMaterialB;         //
   G4String fGearPlateMaterial;                //
   G4String fTargetMountPlateMaterial;         //
   G4String fElectroBoxMaterial;               //
   G4String fPSBoltMaterial;                   //
   G4String fShieldCoverMaterial;              //
   G4String fMagnetCoverMaterial;              //
   G4String fColdFingerMaterial;               //
   G4String fS3CableCaseMaterial;              //
   G4String fConicalCollimatorMaterial;        //11/8
   G4String fXRayInsertMaterial;               //

   //-----------------------------
   // Dimensions of Target Chamber
   //-----------------------------
   static G4double fTargetChamberDistanceFromTarget;
	std::array<G4int, 6> fTargetWheelTargetPos = {{ 0, 0, 1, 0, 1, 0 }};
   G4double fTargetWheelFaceZ{0. * CLHEP::mm};

   G4double fTargetChamberCylinderInnerRadius;
   G4double fTargetChamberCylinderOuterRadius;
   G4double fTargetChamberCylinderLength;
   G4double fTargetChamberCylinderFlatRadius;
   G4double fTargetChamberMELTabRad;
   G4double fTargetChamberMELTabWidth;
   G4double fTargetChamberMagnetIndent;
   G4double fTargetChamberCylinderLipInnerRadius;
   G4double fTargetChamberCylinderLipThickness;
   G4double fTargetChamberCylinderDetOuterRadius;
   G4double fTargetChamberCylinderDetThickness;
   G4double fTargetChamberCylinderHoleLength;
   G4double fTargetChamberCylinderHoleBreadth;
   G4double fTargetChamberCylinderHoleDistanceX;
   G4double fTargetChamberCylinderHoleDistanceY;

   G4double fTargetChamberSphereInnerRadius;
   G4double fTargetChamberSphereOuterRadius;
   G4double fTargetChamberSphereCentre;
   G4double fTargetChamberSphereCuttingOuterRadius;
   G4double fTargetChamberSphereCutOffPoint;
   G4double fFrontRingInnerRadius;
   G4double fFrontRingOuterRadius;
   G4double fFrontRingLength;
   G4double fFrontRingCylinderOuterRadius;
   G4double fFrontRingCylinderLength;

   G4double fConeFrontInnerRad;
   G4double fConeFrontOuterRad;
   G4double fConeLength;

   // --------------------------
   // Dimensions of Target Wheel
   // --------------------------
   G4double fTargetWheelRotation;

   G4double fTargetWheelRadius;
   G4double fTargetWheelThickness;
   G4double fTargetWheelIndentRadius;
   G4double fTargetWheelIndentDepth;

   G4double fTargetOffsetRadius;

   G4double fTargetWheelHoleSpaceAng;
   G4double fTargetWheelTargetHoleRad;
   G4double fTargetWheelExtInnerRad1;
   G4double fTargetWheelExtInnerRad2;
   G4double fTargetWheelExtBaseRad;
   G4double fTargetWheelExtLength;
   G4double fTargetWheelExtTargLip;
   G4double fTargetWheelExtBaseLip;

   // Mount Plate
   G4double fTargetMountPlateRadius;
   G4double fTargetMountPlateThickness;
   G4double fTargetMountPlateCutDepth;
   G4double fTargetMountPlateZOffset;
   G4double fTargetMountPlateHoleRad;

   G4double fMountPlateRodRadius;
   G4double fMountPlateRodPlaceR;
   G4double fMountPlateRodLength;

   G4double fMountPlateRodZoffset;

   //Collimator
   G4double fTargColRadius1;
   G4double fTargColRadius2;
   G4double fTargColHoleHalfSep;
   G4double fTargColPlateAng;
   G4double fTargColPlateThick;
   G4double fTargColPlateFlat;
   G4double fTargColPlateRad;
   G4double fTargetWheelOffset;

   // Gears
   G4double fFirstGearRadius;
   G4double fFirstGearPlaneOffset;
   G4double fSecondGearRadius;
   G4double fSecondGearPlaneOffset;
   G4double fThirdGearRadius;
   G4double fThirdGearPlaneOffset;
   G4double fAllGearThickness;
   G4double fAllGearZOffset;
   // Gear Plates
   G4double fGearPlateOneRadius;
   G4double fGearPlateTwoRadius;
   G4double fGearPlateThickness;
   // Gear Stick
   G4double fGearStickLength;
   G4double fGearStickRadius;
   // Targets
   G4double fBunnyTargetThickness;
   G4double fBunnyRodLength;

   G4double fEvapTargetThickness;
   G4double fEvapTargetRadius;
   G4double fEvapTargetInner;
   G4double fEvapTargetTabWidth;
   G4double fEvapTargetTabExtra;
   G4double fEvapTargetTabIndent;
   G4double fEvapTargetScrewRad;
   G4double fEvapTargetScrewHeight;

   //----------------------------
   // Dimensions of Photon Shield
   //----------------------------
   G4double fPhotonShieldFrontRadius;
   G4double fPhotonShieldBackRadius;
   G4double fPhotonShieldInnerRadius;
   G4double fPhotonShieldLength;
   G4double fPhotonShieldBackFacePos;

   G4double fPhotonShieldLayerOneThickness;
   G4double fPhotonShieldLayerTwoThickness;
   G4double fPhotonShieldLayerThreeThickness;

   // ----------------------------------
   // Dimensions of Photon Shield Clamps
   // ----------------------------------
   G4double fPsTargetClampRadius;
   G4double fPsTargetClampThickness;
   G4double fPsTargetClampExtension;
   G4double fPsTargetClampExtWidth;
   G4double fPsDetClampRadius;
   G4double fPsDetClampThickness;

   // ---------------------------------
   // Dimensions of Photon Shield Bolts
   // ---------------------------------
   G4double fPsBoltHoleRadius;
   G4double fPsBoltPlaceRadius;

   G4double fPsBoltLength;
   G4double fPsBoltRadius;
   G4double fPsBoltBackShield;
   G4double fPsNutLength;
   G4double fPsNutRad;

   //---------------------
   // Dimensions of Magnet
   //---------------------
   G4double fPlateOneThickness;
   G4double fPlateOneLength;
   G4double fPlateOneHeight;
   G4double fPlateOneLowerHeight;

   G4double fPlateOneEdgeX;
   G4double fPlateTwoThickness;
   G4double fPlateTwoLength;
   G4double fPlateTwoHeight;

   G4double fMagDistanceFromTarget;
   G4double fCuttingBoxAngle;

   // ------------------------------------
   // Dimensions of Magnet Clamp (Chamber)
   // ------------------------------------
   G4double fMclampChamberThickness;
   G4double fMclampChamberLength;
   G4double fMclampChamberHeight;
   G4double fMagnetRadiusClampOffset;

   // ------------------------------------------
   // Dimensions of Magnet Clamp (Photon Shield)
   // ------------------------------------------
   G4double fPsClampThicknessTop;
   G4double fPsClampThicknessBottom;
   G4double fPsClampThicknessCut;
   G4double fPsClampThicknessInner;
   G4double fPsClampThicknessOuter;
   G4double fPsClampThicknessFace;
   G4double fPsClampTabHeight;

   // -------------------------
   // Dimensions of ElectroBox
   // -------------------------
   G4double fElectroboxOuterBoxLength;
   G4double fElectroboxInnerBoxLength;
   G4double fElectroboxMidLength;
   G4double fElectroboxLipLength;
   G4double fElectroboxLipInnerRadius;
   G4double fElectroboxZOffset;

   // -------------------------
   // Dimensions of ColdFinger
   // -------------------------
   G4double fColdfingerThickness;
   G4double fColdfingerLength;
   G4double fColdfingerWidth;
   G4double fColdfingerZOffset;
   G4double fColdfingerHoleRadius;

   // ------------------------------
   // Dimensions of Kapton Coating
   // ------------------------------
   G4double fShieldCoatingThickness;

   // ---------------------------
   // individual offsets for visualisation
   // ---------------------------
   G4double fFrontDomeOffset;
   G4double fMiddleRingOffset;
   G4double fBackDetAndAlBoxOffset;

   // -----------------------
   // Dimensions of Beam Pipe
   // -----------------------
   G4double fPipeInnerRadius;
   G4double fPipeOuterRadius;
   G4double fPipeZLength;
   G4double fPipeZOffset;

   // --------------------------------
   // Dimensions of Conical Collimator
   // --------------------------------
   G4double fConicalInnerRad1;
   G4double fConicalInnerRad2;
   G4double fConicalOuterRad1;
   G4double fConicalOuterRad2;
   G4double fConicalLength1;
   G4double fConicalLength2;
   G4double fConicalLength3;
   G4double fConicalTabRad;
   G4double fConicalTabHeight;
   G4double fConicalTabWidth;

   // --------------------------
   // Dimensions of X-ray Insert
   // --------------------------
   G4double fInsertCylOuterRadius;
   G4double fInsertCylInnerRadius;
   G4double fInsertCylZLength;
   G4double fInsertEdgeZLength;

   G4double fInsertHoleRadius;

   G4double fInsertConeZLength;
   G4double fInsertConeR1;
   G4double fInsertConeR2;
   G4double fInsertEdgeAngle;

   // --------------------------
   // Dimensions of Source
   // --------------------------

   G4double fEffSourceHolderHeight;
   G4double fEffSourceHolderRad;
   G4double fEffSourceHolderInner;
   G4double fEffSourceHoldPocket;
   G4double fEffSourceHoldPocketRad;
   G4double fEffSourceScrewRad;
   G4double fEffSourceScrewHeight;

   G4double fMEezagSourceHeight;
   G4double fMEezagSourceRad;
   G4double fMEezagSourceWall;
   G4double fMEezagSourceRim;
   G4double fMEezagSourceZ;
   G4double fMEezagSourceActiveZ;
   G4double fMEezagSourceAcrylicLength;

   G4double fBracketSideLength;
   G4double fBracketBackLength;
   G4double fBracketDepth;
   G4double fBracketBackWidth;
   G4double fBracketSideWidth;
   G4double fBracketoffset;
   G4double fBracketAngle;
   G4double fBracketScrewRad;
   G4double fBracketScrewHeight;

   // --------------------------
   // Construction option bools
   // --------------------------

   bool fBuildMEL{false};
   bool fEmpty{false};
   bool fLight{false};
   bool fBunny{false};
   bool fEff{false};
   bool fSource{false};
   bool fCol{false};

   //////////////////////////////////////////////////////
   // internal methods and functions in ApparatusSpiceTargetChamber::Build()
   //////////////////////////////////////////////////////

   // methods
   void BuildTargetChamberFrontRing();
   void BuildTargetChamberSphere();
   void BuildTargetChamberCylinderDownstream(bool);
   void BuildTargetWheel();
   void BuildTargetWheelSimple();
   void BuildTargetWheelGears();
   void BuildTargetWheelGearPlates();
   void BuildBunnyTarget();
   void BuildEvapTarget();
   void BuildCollimator();
   void BuildGearStick();
   void BuildBiasPlate();
   void BuildTargetMountPlate();
   void BuildPhotonShield();
   void BuildPhotonShieldClamps();
   void BuildPhotonShieldClampBolts();
   void BuildCollectorMagnetandCover();
   void BuildMagnetCovering(G4VSolid*);
   void BuildMagnetClampChamber();
   void BuildElectroBox();
   void BuildShieldCovering();
   void BuildColdFinger();
   void BuildS3CableHolder();
   void BuildConicalCollimator();
   void BuildXrayInsert();

   void BuildPlaceEfficiencySource(G4LogicalVolume*);
   void BuildMEezagSource(G4LogicalVolume*) const;
   void BuildPlaceSource(G4LogicalVolume*);

   void PlaceTargetChamberFrontRing(G4LogicalVolume*);
   void PlaceTargetChamberSphere(G4LogicalVolume*);
   void PlaceTargetChamberCylinderDownstream(G4LogicalVolume*);
   void PlaceTargetWheel(G4LogicalVolume*);
   void PlaceTargetWheelGears(G4LogicalVolume*);
   void PlaceTargetWheelGearPlates(G4LogicalVolume*);
   void PlaceBunnyTargets(G4LogicalVolume*);
   void PlaceBunnyTarget(G4int, G4LogicalVolume*);
   void PlaceEvapTargets(G4LogicalVolume*);
   void PlaceEvapTarget(G4int, G4LogicalVolume*);
   void PlaceCollimator(G4LogicalVolume*);
   void PlaceGearStick(G4LogicalVolume*);
   void PlaceTargetMountPlate(G4LogicalVolume*);
   void PlaceBiasPlate(G4LogicalVolume*);
   void PlacePhotonShield(G4LogicalVolume*);
   void PlacePhotonShieldClamps(G4LogicalVolume*);
   void PlacePhotonShieldClampBolts(G4int, G4LogicalVolume*);
   void PlaceMagnetAndCover(G4int, G4LogicalVolume*);
   void PlaceElectroBox(G4LogicalVolume*);
   void PlaceShieldCovering(G4LogicalVolume*);
   void PlaceColdFinger(G4LogicalVolume*);
   void PlaceS3CableHolder(G4LogicalVolume*);
   void PlaceConicalCollimator(G4LogicalVolume*);
   void PlaceXrayInsert(G4LogicalVolume*);
};
#endif
