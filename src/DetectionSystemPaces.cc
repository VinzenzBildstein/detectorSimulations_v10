#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemPaces.hh"

#include "G4SystemOfUnits.hh"   // new version geant4.10 requires units

DetectionSystemPaces::DetectionSystemPaces()
   /*Measurement status*/
   //~   estimate
   //OK  confirmed measurement
   // Cut clearance
   : fCutClearance(0.01 * mm),
   // Aluminum hemisphere
   fAluminumHemisphereInnerRadius((82.50 / 2) * mm),   //OK
   fAluminumHemisphereOuterRadius((88.86 / 2) * mm),   //OK
   fAluminumHemisphereBeamHoleRadius((34.50 / 2) * mm),   //~
   fAluminumHemisphereBeamHoleRimHeight((3.0) * mm),         //~
   // Annulus that sits in front of silicon crystal
   fAluminumAnnulusTopInnerRadius((15.86 / 2) * mm),   //OK
   fAluminumAnnulusTopOuterRadius((31.71 / 2) * mm),   //OK
   fAluminumAnnulusTopThickness((2.37) * mm),        //OK
   // Annulus that sits in back of silicon crystal
   fAluminumAnnulusBotInnerRadius((15.86 / 2) * mm),   //~
   fAluminumAnnulusBotOuterRadius((31.71 / 2) * mm),   //~
   fAluminumAnnulusBotThickness((4.50) * mm),        //~
   // Canister holding crystal
   fCanisterInnerRadius((18.50 / 2) * mm),   //~
   fCanisterOuterRadius((21.00 / 2) * mm),   //~
   fCanisterThickness((5.00) * mm),        //~
   // Silicon crystal detector
   fSiliconBlockRadius((15.86 / 2) * mm),   //~
   fSiliconBlockThickness((4.90) * mm),        //~
   fSiliconDeadLayerThickness((0.004) * mm),       //~ Nominal Structure Stopping Power of Window - Ortec, 2u Si equivalent
   // Screws
   fScrewRadius((1.00 / 2) * mm),    //~
   fScrewPlacementRadius((22.50 / 2) * mm),   //~
   // Teflon bases
   fTeflonAnnulusTopInnerRadius((21.00 / 2) * mm),   //~
   fTeflonAnnulusTopOuterRadius((32.00 / 2) * mm),   //~
   fTeflonAnnulusTopThickness((2.22) * mm),        //OK
   fTeflonAnnulusBotInnerRadius((31.65 / 2) * mm),   //~
   fTeflonAnnulusBotOuterRadius((34.50 / 2) * mm),   //~
   fTeflonAnnulusBotThickness((4.00) * mm),        //~
   // delrin hemisphere
   fDelrinHemisphereInnerRadius((84.4) * mm),        //OK
   fDelrinHemisphereOuterRadius((89.4) * mm),        //OK
   fDelrinHemisphereBeamHoleRadius((20.00 / 2) * mm),   //OK
   // placement distances
   fAluminumHemisphereDist((0.) * mm),
   fDelrinHemisphereDist((0.) * mm),
   fAluminumAnnulusTopDist((fAluminumAnnulusTopThickness / 2)),
   fCanisterDist((fAluminumAnnulusTopThickness) + (fCanisterThickness / 2)),
   fSiliconBlockDist((fAluminumAnnulusTopThickness) + (fSiliconBlockThickness / 2)),   // Takes dead layer into account
   fSiliconDeadLayerFrontDist((fAluminumAnnulusTopThickness) + (fSiliconDeadLayerThickness / 2)),
   fSiliconDeadLayerBackDist((fAluminumAnnulusTopThickness) + (fSiliconBlockThickness) - (fSiliconDeadLayerThickness / 2)),
   fTeflonAnnulusTopDist((fAluminumAnnulusTopThickness) + (fSiliconBlockThickness) - (fTeflonAnnulusTopThickness / 2)),
   fTeflonAnnulusBotDist((fTeflonAnnulusTopDist) + (fTeflonAnnulusTopThickness / 2) + (fTeflonAnnulusBotThickness / 2)),
   fAluminumAnnulusBotDist((fTeflonAnnulusTopDist) + (fTeflonAnnulusTopThickness / 2) + (fAluminumAnnulusBotThickness / 2)),
   // assembly volume sizes
   fDetectorAssemblyRadius((34.50 / 2) * mm),   //?
   fDetectorAssemblyThickness((15.0) * mm)        //?
{
   // Placement and Orientation scalars
   /*pacesPlacementDistance[5] = {  30.00*mm ,  30.00*mm ,  30.00*mm ,  30.00*mm ,  30.00*mm };
	  pacesPlacementPhi[5] =    {   0.0*deg ,  72.0*deg , 144.0*deg , 216.0*deg , 288.0*deg };
	  pacesPlacementTheta[5] =      {  60.0*deg ,  60.0*deg ,  60.0*deg ,  60.0*deg ,  60.0*deg };
	  pacesOrientationPhi[5] =  {   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg };
	  pacesOrientationTheta[5] =    {   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg };*/

   // distance is deduced from Alumium hemisphere assembly
   fPacesPlacementDistance[0] = 29.50 * mm;   // 29.5*mm Deduced from Alumium hemisphere
   fPacesPlacementDistance[1] = 29.50 * mm;
   fPacesPlacementDistance[2] = 29.50 * mm;
   fPacesPlacementDistance[3] = 29.50 * mm;
   fPacesPlacementDistance[4] = 29.50 * mm;

   // phi angle from measurement
   fPacesPlacementPhi[0] = 0.0000 * deg;      //OK
   fPacesPlacementPhi[1] = -68.7575 * deg;    //OK
   fPacesPlacementPhi[2] = -138.4082 * deg;   //OK
   fPacesPlacementPhi[3] = -212.1280 * deg;   //OK
   fPacesPlacementPhi[4] = -287.4275 * deg;   //OK

   // theta angle is guess estimate
   fPacesPlacementTheta[0] = (180 - 60.0) * deg;   //~
   fPacesPlacementTheta[1] = (180 - 60.0) * deg;   //~
   fPacesPlacementTheta[2] = (180 - 60.0) * deg;   //~
   fPacesPlacementTheta[3] = (180 - 60.0) * deg;   //~
   fPacesPlacementTheta[4] = (180 - 60.0) * deg;   //~

   fPacesOrientationPhi[0] = 0.0 * deg;
   fPacesOrientationPhi[1] = 0.0 * deg;
   fPacesOrientationPhi[2] = 0.0 * deg;
   fPacesOrientationPhi[3] = 0.0 * deg;
   fPacesOrientationPhi[4] = 0.0 * deg;

   fPacesOrientationTheta[0] = 0.1780 * deg;    //OK
   fPacesOrientationTheta[1] = 0.8275 * deg;    //OK
   fPacesOrientationTheta[2] = -0.2579 * deg;   //OK
   fPacesOrientationTheta[3] = 0.2992 * deg;    //OK
   fPacesOrientationTheta[4] = 0.1934 * deg;    //OK
}

DetectionSystemPaces::~DetectionSystemPaces()
{
   // logical volumes
   delete fAluminumHemisphereLog;
   delete fAluminumAnnulusTopLog;
   delete fAluminumAnnulusBotLog;
   delete fCanisterLog;
   delete fSiliconBlockLog;
   delete fSiliconDeadLayerLog;
   delete fTeflonAnnulusTopLog;
   delete fTeflonAnnulusBotLog;
   delete fDelrinHemisphereLog;
}

G4int DetectionSystemPaces::Build()
{
   // Build assembly volumes
   fAssembly         = new G4AssemblyVolume();
   fAssemblyDetector = new G4AssemblyVolume();
   fAssemblySilicon  = new G4AssemblyVolume();

   // Add silicon assembly
   AddSiliconBlock();
   AddSiliconDeadLayer();

   // Add detector assembly
   AddCanister();
   AddAluminumAnnulusTop();
   AddTeflonAnnulusTop();
   AddTeflonAnnulusBot();
   AddAluminumAnnulusBot();
   AddScrews();

   // Add total assembly
   AddAluminumHemisphere();

   // Assemble the ensemble
   CombineAssemblySilicon();
   CombineAssemblyDetector();
   CombineAssembly();

   return 1;
}   //end ::Build

G4int DetectionSystemPaces::PlaceDetector(G4LogicalVolume* expHallLog, G4int ndet)
{
   //place detectors
   if(ndet > 5 || ndet < 0) { ndet = 5; }
   G4ThreeVector     Ta;
   for(G4int i = 0; i < ndet; i++) {
      G4double dDist  = fPacesPlacementDistance[i];
      G4double dPhi   = fPacesPlacementPhi[i];
      G4double dTheta = fPacesPlacementTheta[i];
      auto*    Ra     = new G4RotationMatrix;
      Ta.setX(dDist * cos(dPhi) * sin(dTheta));
      Ta.setY(dDist * sin(dPhi) * sin(dTheta));
      Ta.setZ(dDist * 1.0 * cos(dTheta));

      G4double      oriPhi   = dPhi + fPacesOrientationPhi[i] + M_PI / 2;   //plus 90 deg
      G4double      oriTheta = dTheta + fPacesOrientationTheta[i];
      G4ThreeVector yprimeaxis(cos(oriPhi), sin(oriPhi), 0);
      Ra->set(yprimeaxis, oriTheta);

      fAssemblySilicon->MakeImprint(expHallLog, Ta, Ra, i + 1);
      fAssemblyDetector->MakeImprint(expHallLog, Ta, Ra, i * 20);
   }

   //place hemisphere
   auto* R0 = new G4RotationMatrix;
   G4ThreeVector     T0;
   fAssembly->MakeImprint(expHallLog, T0, R0, 5 * 40);

   return 1;
}   //end ::PlaceDetector

G4int DetectionSystemPaces::CombineAssemblySilicon()
{
   auto* Ra = new G4RotationMatrix;
   G4ThreeVector     Ta;
   Ta.setZ(fSiliconBlockDist);
   fAssemblySilicon->AddPlacedVolume(fSiliconBlockLog, Ta, Ra);

   return 1;
}   //end ::CombineAssemblySilicon

G4int DetectionSystemPaces::CombineAssemblyDetector()
{
   auto* Ra = new G4RotationMatrix;
   G4ThreeVector     Ta;
   Ta.setZ(fSiliconDeadLayerFrontDist);
   fAssemblyDetector->AddPlacedVolume(fSiliconDeadLayerLog, Ta, Ra);
   Ta.setZ(fSiliconDeadLayerBackDist);
   fAssemblyDetector->AddPlacedVolume(fSiliconDeadLayerLog, Ta, Ra);
   Ta.setZ(fCanisterDist);
   fAssemblyDetector->AddPlacedVolume(fCanisterLog, Ta, Ra);
   Ta.setZ(fAluminumAnnulusTopDist);
   fAssemblyDetector->AddPlacedVolume(fAluminumAnnulusTopLog, Ta, Ra);
   Ta.setZ(fAluminumAnnulusBotDist);
   fAssemblyDetector->AddPlacedVolume(fAluminumAnnulusBotLog, Ta, Ra);
   Ta.setZ(fTeflonAnnulusTopDist);
   fAssemblyDetector->AddPlacedVolume(fTeflonAnnulusTopLog, Ta, Ra);
   Ta.setZ(fTeflonAnnulusBotDist);
   fAssemblyDetector->AddPlacedVolume(fTeflonAnnulusBotLog, Ta, Ra);

   return 1;
}   //end ::CombineAssemblyDetector

G4int DetectionSystemPaces::CombineAssembly()
{
   //place aluminum hemisphere
   auto* Ra = new G4RotationMatrix;
   G4ThreeVector     Ta;
   Ta.setZ(fAluminumHemisphereDist);
   fAssembly->AddPlacedVolume(fAluminumHemisphereLog, Ta, Ra);

   return 1;
}   //end ::CombineAssembly

//Add silicon detector
G4int DetectionSystemPaces::AddSiliconBlock()
{
   //material
   G4Material* material = G4Material::GetMaterial("Silicon");
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Silicon"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* siliconBlockVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
   siliconBlockVisAtt->SetVisibility(true);

   //silicon block (without dead layer, smaller than actual)
   G4double deadLayerCut = fSiliconDeadLayerThickness + fCutClearance;
   G4double innerRadius  = 0. * mm;
   G4double outerRadius  = fSiliconBlockRadius;
   G4double halfLengthZ  = fSiliconBlockThickness / 2.0 - deadLayerCut;

   //primitive volume
   auto* siliconBlock = new G4Tubs("siliconBlock", innerRadius, outerRadius, halfLengthZ, 0 * M_PI, 2 * M_PI);

   //logical volume
   if(fSiliconBlockLog == nullptr) {
      fSiliconBlockLog = new G4LogicalVolume(siliconBlock, material, "pacesSiliconBlockLog", nullptr, nullptr, nullptr);   // Renamed from "siliconBlockLog" to "pacesSiliconBlockLog"
      fSiliconBlockLog->SetVisAttributes(siliconBlockVisAtt);
   }

   return 1;
}   //end ::end AddSiliconBlock

//Add detector dead layer
G4int DetectionSystemPaces::AddSiliconDeadLayer()
{
   //material
   G4Material* material = G4Material::GetMaterial("Silicon");   // "Nominal Structure Stopping Power of Window-Ortec, 2u Si"
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Silicon"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* siliconDeadLayerVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
   siliconDeadLayerVisAtt->SetVisibility(true);

   //dead layer
   G4double innerRadius = 0. * mm;
   G4double outerRadius = fSiliconBlockRadius;
   G4double halfLengthZ = fSiliconDeadLayerThickness / 2.0;

   //primitive volume
   auto* siliconDeadLayer = new G4Tubs("siliconDeadLayer", innerRadius, outerRadius, halfLengthZ, 0 * M_PI, 2 * M_PI);

   //logical volume
   if(fSiliconDeadLayerLog == nullptr) {
      fSiliconDeadLayerLog = new G4LogicalVolume(siliconDeadLayer, material, "siliconDeadLayerLog", nullptr, nullptr, nullptr);
      fSiliconDeadLayerLog->SetVisAttributes(siliconDeadLayerVisAtt);
   }

   return 1;
}   //end ::end AddSiliconDeadLayer

//Add detector canister
G4int DetectionSystemPaces::AddCanister()
{
   //material
   G4Material* material = G4Material::GetMaterial("Silicon");   // would be Silicon, Zhimin
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Silicon"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* canisterVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
   canisterVisAtt->SetVisibility(true);

   //main annulus
   G4double innerRadius = fCanisterInnerRadius + fCutClearance;
   G4double outerRadius = fCanisterOuterRadius;
   G4double halfLengthZ = fCanisterThickness / 2.0;

   //primitive volume
   auto* canister = new G4Tubs("canister", innerRadius, outerRadius, halfLengthZ, 0 * M_PI, 2 * M_PI);

   //logical volume
   if(fCanisterLog == nullptr) {
      fCanisterLog = new G4LogicalVolume(canister, material, "canisterLog", nullptr, nullptr, nullptr);
      fCanisterLog->SetVisAttributes(canisterVisAtt);
   }

   return 1;
}   //end ::end AddCanister

//Add annulus that is on top of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusTop()
{
   //material
   G4Material* material = G4Material::GetMaterial("Aluminum");
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Aluminum"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* aluminumAnnulusTopVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.0));
   aluminumAnnulusTopVisAtt->SetVisibility(true);

   //main annulus
   G4double innerRadius = fAluminumAnnulusTopInnerRadius;
   G4double outerRadius = fAluminumAnnulusTopOuterRadius;
   G4double halfLengthZ = fAluminumAnnulusTopThickness / 2.0;
   //screw cuts
   G4double cutInnerRadius = 0. * mm;
   G4double cutOuterRadius = fScrewRadius;
   G4double cutHalfLengthZ = 20. * mm;

   //primitive volume
   auto* alAnToWithoutCuts = new G4Tubs("alAnToWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0 * M_PI, 2 * M_PI);
   auto* cutScrewHole      = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0 * M_PI, 2 * M_PI);

   //cut out screw holes
   G4double      cI = 0.;
   G4double      cR  = fScrewPlacementRadius;
   G4double      cDt = 45. * deg;
   G4ThreeVector moveCut;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnTo0 = new G4SubtractionSolid("alAnTo0", alAnToWithoutCuts, cutScrewHole, nullptr, moveCut);
   cI            = 1;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnTo1 = new G4SubtractionSolid("alAnTo1", alAnTo0, cutScrewHole, nullptr, moveCut);
   cI            = 2;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnTo2 = new G4SubtractionSolid("alAnTo2", alAnTo1, cutScrewHole, nullptr, moveCut);
   cI            = 3;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnTo3 = new G4SubtractionSolid("alAnTo3", alAnTo2, cutScrewHole, nullptr, moveCut);
   cI            = 4;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnTo4 = new G4SubtractionSolid("alAnTo4", alAnTo3, cutScrewHole, nullptr, moveCut);
   cI            = 5;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnTo5 = new G4SubtractionSolid("alAnTo5", alAnTo4, cutScrewHole, nullptr, moveCut);
   cI            = 6;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnTo6 = new G4SubtractionSolid("alAnTo6", alAnTo5, cutScrewHole, nullptr, moveCut);
   //final cut
   cI = 7;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* aluminumAnnulusTop = new G4SubtractionSolid("aluminumAnnulusTop", alAnTo6, cutScrewHole, nullptr, moveCut);

   //logical volume
   if(fAluminumAnnulusTopLog == nullptr) {
      fAluminumAnnulusTopLog = new G4LogicalVolume(aluminumAnnulusTop, material, "aluminumAnnulusTopLog", nullptr, nullptr, nullptr);
      fAluminumAnnulusTopLog->SetVisAttributes(aluminumAnnulusTopVisAtt);
   }

   return 1;
}   //end ::end AddAluminumAnnulusTop

//Add annulus that is on bottom of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusBot()
{
   //material
   G4Material* material = G4Material::GetMaterial("Aluminum");
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Aluminum"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* aluminumAnnulusBotVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.0));
   aluminumAnnulusBotVisAtt->SetVisibility(true);

   //main annulus
   G4double innerRadius = fAluminumAnnulusBotInnerRadius;
   G4double outerRadius = fAluminumAnnulusBotOuterRadius;
   G4double halfLengthZ = fAluminumAnnulusBotThickness / 2.0;
   //screw cuts
   G4double cutInnerRadius = 0. * mm;
   G4double cutOuterRadius = fScrewRadius;
   G4double cutHalfLengthZ = 20. * mm;

   //primitive volume
   auto* alAnBoWithoutCuts = new G4Tubs("alAnBoWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0 * M_PI, 2 * M_PI);
   auto* cutScrewHole      = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0 * M_PI, 2 * M_PI);

   //cut out screw holes
   G4double      cI  = 0.;
   G4double      cR  = fScrewPlacementRadius;
   G4double      cDt = 45. * deg;
   G4ThreeVector moveCut;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnBo0 = new G4SubtractionSolid("alAnBo0", alAnBoWithoutCuts, cutScrewHole, nullptr, moveCut);
   cI            = 1;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnBo1 = new G4SubtractionSolid("alAnBo1", alAnBo0, cutScrewHole, nullptr, moveCut);
   cI            = 2;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnBo2 = new G4SubtractionSolid("alAnBo2", alAnBo1, cutScrewHole, nullptr, moveCut);
   cI            = 3;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnBo3 = new G4SubtractionSolid("alAnBo3", alAnBo2, cutScrewHole, nullptr, moveCut);
   cI            = 4;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnBo4 = new G4SubtractionSolid("alAnBo4", alAnBo3, cutScrewHole, nullptr, moveCut);
   cI            = 5;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnBo5 = new G4SubtractionSolid("alAnBo5", alAnBo4, cutScrewHole, nullptr, moveCut);
   cI            = 6;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* alAnBo6 = new G4SubtractionSolid("alAnBo6", alAnBo5, cutScrewHole, nullptr, moveCut);
   //final cut
   cI = 7;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* aluminumAnnulusBot = new G4SubtractionSolid("aluminumAnnulusBot", alAnBo6, cutScrewHole, nullptr, moveCut);

   //logical volume
   if(fAluminumAnnulusBotLog == nullptr) {
      fAluminumAnnulusBotLog = new G4LogicalVolume(aluminumAnnulusBot, material, "aluminumAnnulusBotLog", nullptr, nullptr, nullptr);
      fAluminumAnnulusBotLog->SetVisAttributes(aluminumAnnulusBotVisAtt);
   }

   return 1;
}   //end ::end AddAluminumAnnulusBot

//Add annulus that is on side of silicon detector
G4int DetectionSystemPaces::AddTeflonAnnulusTop()
{
   //material
   G4Material* material = G4Material::GetMaterial("Teflon");
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Teflon"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* teflonAnnulusTopVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
   teflonAnnulusTopVisAtt->SetVisibility(true);

   //main annulus
   G4double innerRadius = fTeflonAnnulusTopInnerRadius;
   G4double outerRadius = fTeflonAnnulusTopOuterRadius;
   G4double halfLengthZ = fTeflonAnnulusTopThickness / 2.0;
   //screw cuts
   G4double cutInnerRadius = 0. * mm;
   G4double cutOuterRadius = fScrewRadius;
   G4double cutHalfLengthZ = 20. * mm;

   //primitive volume
   auto* teAnToWithoutCuts = new G4Tubs("teAnToWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0 * M_PI, 2 * M_PI);
   auto* cutScrewHole      = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0 * M_PI, 2 * M_PI);

   //cut out screw holes
   G4double      cI  = 0.;
   G4double      cR  = fScrewPlacementRadius;
   G4double      cDt = 45. * deg;
   G4ThreeVector moveCut;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnTo0 = new G4SubtractionSolid("teAnTo0", teAnToWithoutCuts, cutScrewHole, nullptr, moveCut);
   cI            = 1;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnTo1 = new G4SubtractionSolid("teAnTo1", teAnTo0, cutScrewHole, nullptr, moveCut);
   cI            = 2;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnTo2 = new G4SubtractionSolid("teAnTo2", teAnTo1, cutScrewHole, nullptr, moveCut);
   cI            = 3;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnTo3 = new G4SubtractionSolid("teAnTo3", teAnTo2, cutScrewHole, nullptr, moveCut);
   cI            = 4;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnTo4 = new G4SubtractionSolid("teAnTo4", teAnTo3, cutScrewHole, nullptr, moveCut);
   cI            = 5;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnTo5 = new G4SubtractionSolid("teAnTo5", teAnTo4, cutScrewHole, nullptr, moveCut);
   cI            = 6;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnTo6 = new G4SubtractionSolid("teAnTo6", teAnTo5, cutScrewHole, nullptr, moveCut);
   //final cut
   cI = 7;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teflonAnnulusTop = new G4SubtractionSolid("teflonAnnulusTop", teAnTo6, cutScrewHole, nullptr, moveCut);

   //logical volume
   if(fTeflonAnnulusTopLog == nullptr) {
      fTeflonAnnulusTopLog = new G4LogicalVolume(teflonAnnulusTop, material, "teflonAnnulusTopLog", nullptr, nullptr, nullptr);
      fTeflonAnnulusTopLog->SetVisAttributes(teflonAnnulusTopVisAtt);
   }

   return 1;
}   //end ::end AddTeflonAnnulusTop

//Add annulus that is on bottom of silicon detector
G4int DetectionSystemPaces::AddTeflonAnnulusBot()
{
   //material
   G4Material* material = G4Material::GetMaterial("Teflon");
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Teflon"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* teflonAnnulusBotVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
   teflonAnnulusBotVisAtt->SetVisibility(true);

   //main annulus
   G4double innerRadius = fTeflonAnnulusBotInnerRadius;
   G4double outerRadius = fTeflonAnnulusBotOuterRadius;
   G4double halfLengthZ = fTeflonAnnulusBotThickness / 2.0;
   //screw cuts
   G4double cutInnerRadius = 0. * mm;
   G4double cutOuterRadius = fScrewRadius;
   G4double cutHalfLengthZ = 20. * mm;

   //primitive volume
   auto* teAnBoWithoutCuts = new G4Tubs("teAnBoWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0 * M_PI, 2 * M_PI);
   auto* cutScrewHole      = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0 * M_PI, 2 * M_PI);

   //cut out screw holes
   G4double      cI  = 0.;
   G4double      cR  = fScrewPlacementRadius;
   G4double      cDt = 45. * deg;
   G4ThreeVector moveCut;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnBo0 = new G4SubtractionSolid("teAnBo0", teAnBoWithoutCuts, cutScrewHole, nullptr, moveCut);
   cI            = 1;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnBo1 = new G4SubtractionSolid("teAnBo1", teAnBo0, cutScrewHole, nullptr, moveCut);
   cI            = 2;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnBo2 = new G4SubtractionSolid("teAnBo2", teAnBo1, cutScrewHole, nullptr, moveCut);
   cI            = 3;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnBo3 = new G4SubtractionSolid("teAnBo3", teAnBo2, cutScrewHole, nullptr, moveCut);
   cI            = 4;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnBo4 = new G4SubtractionSolid("teAnBo4", teAnBo3, cutScrewHole, nullptr, moveCut);
   cI            = 5;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnBo5 = new G4SubtractionSolid("teAnBo5", teAnBo4, cutScrewHole, nullptr, moveCut);
   cI            = 6;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teAnBo6 = new G4SubtractionSolid("teAnBo6", teAnBo5, cutScrewHole, nullptr, moveCut);
   //final cut
   cI = 7;
   moveCut.setX(cR * cos(cI * cDt));
   moveCut.setY(cR * sin(cI * cDt));
   auto* teflonAnnulusBot = new G4SubtractionSolid("teflonAnnulusBot", teAnBo6, cutScrewHole, nullptr, moveCut);

   //logical volume
   if(fTeflonAnnulusBotLog == nullptr) {
      fTeflonAnnulusBotLog = new G4LogicalVolume(teflonAnnulusBot, material, "teflonAnnulusBotLog", nullptr, nullptr, nullptr);
      fTeflonAnnulusBotLog->SetVisAttributes(teflonAnnulusBotVisAtt);
   }

   return 1;
}   //end ::end AddTeflonAnnulusBot

//Add the aluminum hemisphere, housing for PACES
G4int DetectionSystemPaces::AddAluminumHemisphere()
{
   //material
   G4Material* material = G4Material::GetMaterial("Aluminum");
   if(material == nullptr) {
      G4cout << " ----> Material "
             << "Aluminum"
             << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   //vis attributes
   auto* aluminumHemisphereVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
   aluminumHemisphereVisAtt->SetVisibility(true);

   //main hemisphere shell
   G4double innerRadius = fAluminumHemisphereInnerRadius;
   G4double outerRadius = fAluminumHemisphereOuterRadius;
   G4double thickness   = std::abs(outerRadius - innerRadius);
   G4double startPhi    = M_PI / 2;
   G4double endPhi      = M_PI - startPhi;
   //beam hole rim
   G4double rimInnerRadius = std::abs(innerRadius - fAluminumHemisphereBeamHoleRimHeight);
   G4double rimOuterRadius = innerRadius;
   G4double rimEndPhi      = asin((fAluminumHemisphereBeamHoleRadius + thickness) / innerRadius);
   G4double rimStartPhi    = M_PI - rimEndPhi;
   //beam hole cut
   G4double beamHoleInnerRadius = 0. * mm;
   G4double beamHoleOuterRadius = fAluminumHemisphereBeamHoleRadius;
   G4double beamHoleHalfLengthZ = outerRadius;
   //cuts where detectors go
   G4double detectorInnerRadius = 0. * mm;
   G4double detectorOuterRadius = fDetectorAssemblyRadius + fCutClearance;
   G4double detectorHalfLengthZ = fDetectorAssemblyThickness / 2.0;

   //primitive volumes
   auto* alHeShell   = new G4Sphere("alHeShell", innerRadius, outerRadius, 0 * M_PI, 2 * M_PI, startPhi, endPhi);
   auto* alHeRim     = new G4Sphere("alHeRim", rimInnerRadius, rimOuterRadius, 0 * M_PI, 2 * M_PI, rimStartPhi, rimEndPhi);
   auto*   cutBeamHole = new G4Tubs("cutBeamHole", beamHoleInnerRadius, beamHoleOuterRadius, beamHoleHalfLengthZ, 0 * M_PI, 2 * M_PI);
   auto*   cutDetector = new G4Tubs("cutDetector", detectorInnerRadius, detectorOuterRadius, detectorHalfLengthZ, 0 * M_PI, 2 * M_PI);

   //add rim and beam hole
   auto*   R0 = new G4RotationMatrix;
   G4ThreeVector       T0;
   auto*       alHeShellRim     = new G4UnionSolid("alHeShellRim", alHeShell, alHeRim, R0, T0);
   auto* alHeShellRimHole = new G4SubtractionSolid("alHeShellRimHole", alHeShellRim, cutBeamHole, R0, T0);
   //cut out detectors

   G4ThreeVector moveCut;
   G4ThreeVector yprimeaxis;
   auto*         rotateCut = new G4RotationMatrix;
   G4int         dI        = 0;
   //one
   G4double dDist  = fPacesPlacementDistance[dI] + detectorHalfLengthZ;
   G4double dPhi   = fPacesPlacementPhi[dI];
   G4double dTheta = fPacesPlacementTheta[dI];
   moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
   moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
   moveCut.setZ(dDist * 1.0 * cos(dTheta));
   G4double oriPhi   = dPhi + fPacesOrientationPhi[dI] - M_PI / 2;
   G4double oriTheta = dTheta + fPacesOrientationTheta[dI];
   yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
   rotateCut->set(yprimeaxis, oriTheta);
   auto* alHe0 = new G4SubtractionSolid("alHe0", alHeShellRimHole, cutDetector, rotateCut, moveCut);
   //two
   dI     = 1;
   dDist  = fPacesPlacementDistance[dI] + detectorHalfLengthZ;
   dPhi   = fPacesPlacementPhi[dI];
   dTheta = fPacesPlacementTheta[dI];
   moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
   moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
   moveCut.setZ(dDist * 1.0 * cos(dTheta));
   oriPhi   = dPhi + fPacesOrientationPhi[dI] - M_PI / 2;
   oriTheta = dTheta + fPacesOrientationTheta[dI];
   yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
   rotateCut->set(yprimeaxis, oriTheta);
   auto* alHe1 = new G4SubtractionSolid("alHe1", alHe0, cutDetector, rotateCut, moveCut);
   //three
   dI     = 2;
   dDist  = fPacesPlacementDistance[dI] + detectorHalfLengthZ;
   dPhi   = fPacesPlacementPhi[dI];
   dTheta = fPacesPlacementTheta[dI];
   moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
   moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
   moveCut.setZ(dDist * 1.0 * cos(dTheta));
   oriPhi   = dPhi + fPacesOrientationPhi[dI] - M_PI / 2;
   oriTheta = dTheta + fPacesOrientationTheta[dI];
   yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
   rotateCut->set(yprimeaxis, oriTheta);
   auto* alHe2 = new G4SubtractionSolid("alHe2", alHe1, cutDetector, rotateCut, moveCut);
   //four
   dI     = 3;
   dDist  = fPacesPlacementDistance[dI] + detectorHalfLengthZ;
   dPhi   = fPacesPlacementPhi[dI];
   dTheta = fPacesPlacementTheta[dI];
   moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
   moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
   moveCut.setZ(dDist * 1.0 * cos(dTheta));
   oriPhi   = dPhi + fPacesOrientationPhi[dI] - M_PI / 2;
   oriTheta = dTheta + fPacesOrientationTheta[dI];
   yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
   rotateCut->set(yprimeaxis, oriTheta);
   auto* alHe3 = new G4SubtractionSolid("alHe3", alHe2, cutDetector, rotateCut, moveCut);
   //five
   dI     = 4;
   dDist  = fPacesPlacementDistance[dI] + detectorHalfLengthZ;
   dPhi   = fPacesPlacementPhi[dI];
   dTheta = fPacesPlacementTheta[dI];
   moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
   moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
   moveCut.setZ(dDist * 1.0 * cos(dTheta));
   oriPhi   = dPhi + fPacesOrientationPhi[dI] - M_PI / 2;
   oriTheta = dTheta + fPacesOrientationTheta[dI];
   yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
   rotateCut->set(yprimeaxis, oriTheta);
   auto* aluminumHemisphere = new G4SubtractionSolid("aluminumHemisphere", alHe3, cutDetector, rotateCut, moveCut);
   //  //end

   //logical volume
   if(fAluminumHemisphereLog == nullptr) {
      fAluminumHemisphereLog = new G4LogicalVolume(aluminumHemisphere, material, "aluminumHemisphereLog", nullptr, nullptr, nullptr);
      fAluminumHemisphereLog->SetVisAttributes(aluminumHemisphereVisAtt);
   }

   return 1;
}   //end ::end AddAluminumHemisphere

//Add screws
G4int DetectionSystemPaces::AddScrews() const
{
   return 1;
}   //end ::end AddScrews
