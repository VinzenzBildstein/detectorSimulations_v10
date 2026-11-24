#include <cmath>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemGriffin.hh"

#include "G4SystemOfUnits.hh"   // new version geant4.10 requires units

///////////////////////////////////////////////////////////////////////
// The ::DetectionSystemGriffin constructor instatiates all the
// Logical and Physical Volumes used in the detector geometery (found
// in our local files as it contains NDA-protected info), and the
// ::~DetectionSystemGriffin destructor deletes them from the stack
// when they go out of scope
///////////////////////////////////////////////////////////////////////

DetectionSystemGriffin::~DetectionSystemGriffin()
{
   // LogicalVolumes in ConstructNewHeavyMet
   delete fHevimetLog;

   // LogicalVolumes in ConstructNewSuppressorCasing
   delete fLeftSuppressorExtensionLog;
   delete fRightSuppressorExtensionLog;
   delete fLeftSuppressorLog;
   delete fRightSuppressorLog;
   delete fBackQuarterSuppressorLog;

   delete fBackQuarterSuppressorShellLog;
   delete fRightSuppressorShellLog;
   delete fLeftSuppressorShellLog;
   delete fRightSuppressorShellExtensionLog;
   delete fLeftSuppressorShellExtensionLog;
   // LogicalVolumes in ConstructBGOCasing

   delete fBGOCasingLog;
   delete fBackBGOLog;

   // LogicalVolumes in ConstructColdFinger
   delete fCoolingSideBlockLog;
   delete fStructureMatColdFingerLog;
   delete fCoolingBarLog;
   delete fFetAirHoleLog;
   delete fTrianglePostLog;
   delete fExtraColdBlockLog;
   delete fFingerLog;
   delete fEndPlateLog;

   // LogicalVolumes in ConstructComplexDetectorBlock
   delete fGermaniumHoleLog;
   delete fGermaniumBlock1Log;
   delete fInnerDeadLayerLog;
   delete fInterCrystalElectrodeMatBackLog;
   delete fInterCrystalElectrodeMatFrontLog;

   // LogicalVolumes in ConstructBasicDetectorBlock
   delete fGermaniumBlockLog;

   // LogicalVolumes in ConstructDetector
   delete fTankLog;
   delete fTankLid1Log;
   delete fTankLid2Log;
   delete fFingerShellLog;
   delete fRearPlateLog;
   delete fBottomSidePanelLog;
   delete fTopSidePanelLog;
   delete fLeftSidePanelLog;
   delete fRightSidePanelLog;
   delete fLowerLeftTubeLog;
   delete fUpperLeftTubeLog;
   delete fLowerRightTubeLog;
   delete fUpperRightTubeLog;
   delete fLowerLeftConeLog;
   delete fUpperLeftConeLog;
   delete fLowerRightConeLog;
   delete fUpperRightConeLog;
   delete fBottomWedgeLog;
   delete fTopWedgeLog;
   delete fLeftWedgeLog;
   delete fRightWedgeLog;
   delete fBottomBentPieceLog;
   delete fTopBentPieceLog;
   delete fLeftBentPieceLog;
   delete fRightBentPieceLog;
   delete fFrontFaceLog;
}   // end ::~DetectionSystemGriffin

///////////////////////////////////////////////////////////////////////
// ConstructDetectionSystemGriffin builds the DetectionSystemGriffin
// at the origin
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::Build()
{

   fSurfCheck = true;
   G4int det  = 1;

   BuildOneDetector(det);

}   //end ::Build

/// These functions describe a rotation of coordinates x,y,z by theta around the y-axis and then by phi around the z-axis
G4double DetectionSystemGriffin::TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi)
{
   return (x * cos(theta) + z * sin(theta)) * cos(phi) - y * sin(phi);
}

G4double DetectionSystemGriffin::TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi)
{
   return (x * cos(theta) + z * sin(theta)) * sin(phi) + y * cos(phi);
}

G4double DetectionSystemGriffin::TransZ(G4double x, G4double z, G4double theta)
{
   return -x * sin(theta) + z * cos(theta);
}

G4int DetectionSystemGriffin::PlaceEverythingButCrystals(G4LogicalVolume* expHallLog, G4int detectorNumber, G4int positionNumber, G4bool posTigress)
{
   if(expHallLog == nullptr) {
      G4cerr << __PRETTY_FUNCTION__ << ": expHallLog == nullptr!" << std::endl; // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      exit(1);
   }

   G4double theta = fCoords[positionNumber][0] * deg;
   G4double phi   = fCoords[positionNumber][1] * deg;
   G4double alpha = fCoords[positionNumber][2] * deg;   // yaw
   G4double beta  = fCoords[positionNumber][3] * deg;   // pitch
   G4double gamma = fCoords[positionNumber][4] * deg;   // roll

   //In GRIFFIN upstream is now downstream relative to TIGRESS. However, we want to maintain the same numbering scheme as TIGRESS. Net result is that the lampshade angles change by 45 degrees.
   if(posTigress && (positionNumber < 4 || positionNumber > 11)) {
      phi   = fCoords[positionNumber][1] * deg - 45.0 * deg;
      gamma = fCoords[positionNumber][4] * deg - 45.0 * deg;
   }

   auto* rotate = new G4RotationMatrix;   // rotation matrix corresponding to direction vector
   rotate->rotateX(M_PI / 2.0);
   rotate->rotateX(alpha);
   rotate->rotateY(beta);
   rotate->rotateZ(gamma);

   // positioning
   G4double distFromOrigin    = fAirBoxBackLength / 2.0 + fAirBoxFrontLength + fNewRhombiRadius;
   G4double distFromOriginDet = fAirBoxBackLengthDet / 2.0 + fAirBoxFrontLengthDet + fNewRhombiRadiusDet;

   G4double x = 0;
   G4double y = 0;
   G4double z = distFromOriginDet;

   G4ThreeVector move(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));
   // general note: MakeImprint first rotates, then translates
   fAssembly->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);
   fBackAndSideSuppressorShellAssembly->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);

   z = distFromOrigin;   // This isolates the motion of the extension suppressors from the rest of the suppressors.

   move = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));

   fExtensionSuppressorShellAssembly->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);
   fHevimetAssembly->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);

   fCopyNumber = fBackSuppressorCopyNumber + detectorNumber * 4;

   // Back Suppressors
   if(fIncludeBackSuppressors) {
      G4double x0 = (fDetectorTotalWidth / 4.0);
      G4double y0 = (fDetectorTotalWidth / 4.0);
      G4double z0 = (fBackBGOThickness - fCanFaceThickness) / 2.0 + fSuppressorShellThickness + fDetectorTotalLength + fBGOCanSeperation + fShift + fAppliedBackShift + distFromOriginDet;

		std::array<G4RotationMatrix*, 4> rotateBackQuarterSuppressor;
      std::array<G4ThreeVector, 4>     moveBackQuarterSuppressor;

      for(int i = 0; i < 4; i++) {
         rotateBackQuarterSuppressor[i] = new G4RotationMatrix;
         rotateBackQuarterSuppressor[i]->rotateX(M_PI / 2.0 - M_PI / 2.0 * i);
         rotateBackQuarterSuppressor[i]->rotateX(alpha);
         rotateBackQuarterSuppressor[i]->rotateY(beta);
         rotateBackQuarterSuppressor[i]->rotateZ(gamma);

         x = -x0 * pow(-1, floor((i + 1) / 2));
         y = -y0 * pow(-1, floor((i + 2) / 2));
         z = z0;

         moveBackQuarterSuppressor[i] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));

         fSuppressorBackAssembly->MakeImprint(expHallLog, moveBackQuarterSuppressor[i], rotateBackQuarterSuppressor[i], fCopyNumber++, fSurfCheck);
      }
   }

   // Now Side Suppressors
   /////////////////////////////////////////////////////////////////////
   // Note : Left and Right are read from the BACK of the detector
   // Suppressors 1 and 2 cover germanium 1
   // Suppressors 3 and 4 cover germanium 2
   // Suppressors 5 and 6 cover germanium 3
   // Suppressors 7 and 8 cover germanium 4
   /////////////////////////////////////////////////////////////////////
   fCopyNumber    = fRightSuppressorSideCopyNumber + detectorNumber * 4;
   fCopyNumberTwo = fLeftSuppressorSideCopyNumber + detectorNumber * 4;

   // Replacement for sideSuppressorLength
   G4double shellSideSuppressorLength = fSideSuppressorLength + (fSuppressorShellThickness * 2.0);

   // Replacement for suppressorExtensionLength
   G4double shellSuppressorExtensionLength = fSuppressorExtensionLength + (fSuppressorShellThickness * 2.0) * (1.0 / tan(fBentEndAngle) - tan(fBentEndAngle));

   // Replacement for suppressorExtensionAngle: must totally recalculate
   G4double shellSuppressorExtensionAngle = atan(((fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0) / tan(fBentEndAngle) - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * sin(fBentEndAngle)) * tan(fBentEndAngle) - (fSuppressorForwardRadius + fHevimetTipThickness) * sin(fBentEndAngle)) / (shellSuppressorExtensionLength));

   // these two parameters are for shifting the extensions back and out when in their BACK position
   G4double extensionBackShift = fAirBoxFrontLength - (fHevimetTipThickness + shellSuppressorExtensionLength + (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * tan(fBentEndAngle)) * cos(fBentEndAngle);

   G4double extensionRadialShift = extensionBackShift * tan(fBentEndAngle);

   if(fIncludeSideSuppressors) {
		std::array<G4RotationMatrix*, 8> rotateSideSuppressor;
      std::array<G4ThreeVector, 8>     moveInnerSuppressor;

      G4double x0 = fSideBGOThickness / 2.0 + fSuppressorShellThickness + fDetectorTotalWidth / 2.0 + fBGOCanSeperation;
      G4double y0 = (fDetectorTotalWidth / 2.0 + fBGOCanSeperation + fSideBGOThickness / 2.0) / 2.0;
      G4double z0 = (shellSideSuppressorLength / 2.0 - fCanFaceThickness / 2.0 + fBentEndLength + (fBGOCanSeperation + fBGOChoppedTip) / tan(fBentEndAngle) + fShift + fAppliedBackShift - fSuppressorShellThickness / 2.0 + distFromOriginDet);

      for(int i = 0; i < 4; i++) {
         rotateSideSuppressor[2 * i] = new G4RotationMatrix;
         rotateSideSuppressor[2 * i]->rotateZ(M_PI / 2.0);
         rotateSideSuppressor[2 * i]->rotateY(M_PI / 2.0);
         rotateSideSuppressor[2 * i]->rotateX(M_PI / 2.0 - M_PI / 2.0 * i);
         rotateSideSuppressor[2 * i]->rotateX(alpha);
         rotateSideSuppressor[2 * i]->rotateY(beta);
         rotateSideSuppressor[2 * i]->rotateZ(gamma);

         x = x0 * cos((i - 1) * M_PI / 2.0) + y0 * sin((i - 1) * M_PI / 2);
         y = -x0 * sin((i - 1) * M_PI / 2.0) + y0 * cos((i - 1) * M_PI / 2);
         z = z0;

         moveInnerSuppressor[i * 2] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));

         fRightSuppressorCasingAssembly->MakeImprint(expHallLog, moveInnerSuppressor[2 * i], rotateSideSuppressor[2 * i], fCopyNumber++, fSurfCheck);

         rotateSideSuppressor[2 * i + 1] = new G4RotationMatrix;
         rotateSideSuppressor[2 * i + 1]->rotateY(-M_PI / 2.0);
         rotateSideSuppressor[2 * i + 1]->rotateX(-M_PI / 2.0 * (i + 2));
         rotateSideSuppressor[2 * i + 1]->rotateX(alpha);
         rotateSideSuppressor[2 * i + 1]->rotateY(beta);
         rotateSideSuppressor[2 * i + 1]->rotateZ(gamma);

         x = x0 * sin(i * M_PI / 2.0) + y0 * cos(i * M_PI / 2);
         y = x0 * cos(i * M_PI / 2.0) - y0 * sin(i * M_PI / 2);
         z = z0;

         moveInnerSuppressor[i * 2 + 1] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));

         fLeftSuppressorCasingAssembly->MakeImprint(expHallLog, moveInnerSuppressor[2 * i + 1], rotateSideSuppressor[2 * i + 1], fCopyNumberTwo++, fSurfCheck);
      }
   }

   // now we add the side pieces of suppressor that extend out in front of the can when it's in the back position
   fCopyNumber    = fRightSuppressorExtensionCopyNumber + detectorNumber * 4;
   fCopyNumberTwo = fLeftSuppressorExtensionCopyNumber + detectorNumber * 4;

	std::array<G4RotationMatrix*, 8> rotateExtension;
   std::array<G4ThreeVector, 8>     moveInnerExtension;

   G4double x0 = -((fSuppressorExtensionThickness / 2.0 + fSuppressorShellThickness) / cos(fBentEndAngle) + (shellSuppressorExtensionLength / 2.0 - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * tan(fBentEndAngle) / 2.0) * sin(fBentEndAngle) - (fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0) / tan(fBentEndAngle) - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * sin(fBentEndAngle)) * tan(fBentEndAngle));

   G4double y0 = (shellSuppressorExtensionLength * tan(shellSuppressorExtensionAngle) / 2.0 + (fSuppressorForwardRadius + fHevimetTipThickness) * sin(fBentEndAngle) - fBGOCanSeperation * 2.0 - fDetectorPlacementCxn) / 2.0;

   G4double z0 = -fCanFaceThickness / 2.0 - (shellSuppressorExtensionLength / 2.0 - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * tan(fBentEndAngle) / 2.0) * cos(fBentEndAngle) + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0) / tan(fBentEndAngle) - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * sin(fBentEndAngle) + fSuppShift + fSuppressorBackRadius - fSuppressorForwardRadius - fSuppressorShellThickness / 2.0 + distFromOrigin;

   if(fSuppressorPositionSelector == 0 && fIncludeExtensionSuppressors) {
      // If the detectors are forward, put the extensions in the back position
      // the placement of the extensions matches the placement of the sideSuppressor pieces

      x0 += extensionRadialShift;

      z0 += extensionBackShift;

      for(int i = 0; i < 4; i++) {
         rotateExtension[2 * i] = new G4RotationMatrix;
         rotateExtension[2 * i]->rotateZ(M_PI / 2.0);
         rotateExtension[2 * i]->rotateY(fBentEndAngle);
         rotateExtension[2 * i]->rotateX(M_PI / 2.0 - M_PI / 2.0 * i);
         rotateExtension[2 * i]->rotateX(alpha);
         rotateExtension[2 * i]->rotateY(beta);
         rotateExtension[2 * i]->rotateZ(gamma);

         x = x0 * cos((i - 1) * M_PI / 2.0) + y0 * sin((i - 1) * M_PI / 2);
         y = -x0 * sin((i - 1) * M_PI / 2.0) + y0 * cos((i - 1) * M_PI / 2);
         z = z0;

         moveInnerExtension[2 * i] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));
         fRightSuppressorExtensionAssembly->MakeImprint(expHallLog, moveInnerExtension[2 * i], rotateExtension[2 * i], fCopyNumber++, fSurfCheck);

         rotateExtension[2 * i + 1] = new G4RotationMatrix;
         rotateExtension[2 * i + 1]->rotateY(M_PI / 2.0);
         rotateExtension[2 * i + 1]->rotateZ(M_PI / 2.0 + fBentEndAngle);
         rotateExtension[2 * i + 1]->rotateX(-M_PI / 2.0 * i);
         rotateExtension[2 * i + 1]->rotateX(alpha);
         rotateExtension[2 * i + 1]->rotateY(beta);
         rotateExtension[2 * i + 1]->rotateZ(gamma);

         x = x0 * sin(i * M_PI / 2.0) + y0 * cos(i * M_PI / 2);
         y = x0 * cos(i * M_PI / 2.0) - y0 * sin(i * M_PI / 2);
         z = z0;

         moveInnerExtension[2 * i + 1] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));
         fLeftSuppressorExtensionAssembly->MakeImprint(expHallLog, moveInnerExtension[2 * i + 1], rotateExtension[2 * i + 1], fCopyNumberTwo++, fSurfCheck);
      }
   } else if((fSuppressorPositionSelector == 1) && fIncludeExtensionSuppressors) {
      for(int i = 0; i < 4; i++) {
         rotateExtension[2 * i] = new G4RotationMatrix;
         rotateExtension[2 * i]->rotateZ(M_PI / 2.0);
         rotateExtension[2 * i]->rotateY(fBentEndAngle);
         rotateExtension[2 * i]->rotateX(M_PI / 2.0 - M_PI / 2.0 * i);
         rotateExtension[2 * i]->rotateX(alpha);
         rotateExtension[2 * i]->rotateY(beta);
         rotateExtension[2 * i]->rotateZ(gamma);

         x = x0 * cos((i - 1) * M_PI / 2.0) + y0 * sin((i - 1) * M_PI / 2);
         y = -x0 * sin((i - 1) * M_PI / 2.0) + y0 * cos((i - 1) * M_PI / 2);
         z = z0;

         moveInnerExtension[2 * i] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));
         fRightSuppressorExtensionAssembly->MakeImprint(expHallLog, moveInnerExtension[2 * i], rotateExtension[2 * i], fCopyNumber++, fSurfCheck);

         rotateExtension[2 * i + 1] = new G4RotationMatrix;
         rotateExtension[2 * i + 1]->rotateY(M_PI / 2.0);
         rotateExtension[2 * i + 1]->rotateZ(M_PI / 2.0 + fBentEndAngle);
         rotateExtension[2 * i + 1]->rotateX(-M_PI / 2.0 * i);
         rotateExtension[2 * i + 1]->rotateX(alpha);
         rotateExtension[2 * i + 1]->rotateY(beta);
         rotateExtension[2 * i + 1]->rotateZ(gamma);

         x = x0 * sin(i * M_PI / 2.0) + y0 * cos(i * M_PI / 2);
         y = x0 * cos(i * M_PI / 2.0) - y0 * sin(i * M_PI / 2);
         z = z0;

         moveInnerExtension[2 * i + 1] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));
         fLeftSuppressorExtensionAssembly->MakeImprint(expHallLog, moveInnerExtension[2 * i + 1], rotateExtension[2 * i + 1], fCopyNumberTwo++, fSurfCheck);
      }

   }   //end if(detectors back) statement
   return 1;
}   // End PlaceEverythingButCrystals

G4int DetectionSystemGriffin::PlaceDeadLayerSpecificCrystal(G4LogicalVolume* expHallLog, G4int detectorNumber, G4int positionNumber, G4bool posTigress)
{

   G4double theta = fCoords[positionNumber][0] * deg;
   G4double phi   = fCoords[positionNumber][1] * deg;
   G4double alpha = fCoords[positionNumber][2] * deg;   // yaw
   G4double beta  = fCoords[positionNumber][3] * deg;   // pitch
   G4double gamma = fCoords[positionNumber][4] * deg;   // roll

   //In GRIFFIN upstream is now downstream relative to TIGRESS. However, we want to maintain the same numbering scheme as TIGRESS. Net result is that the lampshade angles change by 45 degrees.
   if(posTigress && (positionNumber < 4 || positionNumber > 11)) {
      phi   = fCoords[positionNumber][1] * deg - 45.0 * deg;
      gamma = fCoords[positionNumber][4] * deg - 45.0 * deg;
   }

   // positioning
   G4double distFromOriginDet = fAirBoxBackLengthDet / 2.0 + fAirBoxFrontLengthDet + fNewRhombiRadiusDet;

   G4double x0 = (fGermaniumWidth + fGermaniumSeparation) / 2.0;
   G4double y0 = (fGermaniumWidth + fGermaniumSeparation) / 2.0;
   G4double z0 = fGermaniumLength / 2.0 + fCanFaceThickness / 2.0 + fGermaniumDistFromCanFace + fShift + fAppliedBackShift + distFromOriginDet;

   /////////////////////////////////////////////////////////////////////
   // now we place all 4 of the 1/4 detectors using the LogicalVolume
   // of the 1st, ensuring that they too will have the hole
   /////////////////////////////////////////////////////////////////////
   fCopyNumber = fGermaniumCopyNumber + detectorNumber * 4;

	std::array<G4RotationMatrix*, 4> rotateGermanium;
   std::array<G4ThreeVector, 4>     moveGermanium;

   for(int i = 0; i < 4; i++) {
      rotateGermanium[i] = new G4RotationMatrix;
      rotateGermanium[i]->rotateY(-M_PI / 2.0);
      rotateGermanium[i]->rotateX(M_PI / 2.0 - M_PI / 2.0 * i);
      rotateGermanium[i]->rotateX(alpha);
      rotateGermanium[i]->rotateY(beta);
      rotateGermanium[i]->rotateZ(gamma);

      G4double x = -x0 * pow(-1, floor((i + 1) / 2));
      G4double y = -y0 * pow(-1, floor((i + 2) / 2));
      G4double z = z0;

      moveGermanium[i] = G4ThreeVector(DetectionSystemGriffin::TransX(x, y, z, theta, phi), DetectionSystemGriffin::TransY(x, y, z, theta, phi), DetectionSystemGriffin::TransZ(x, z, theta));

      fGermaniumAssemblyCry[i]->MakeImprint(expHallLog, moveGermanium[i], rotateGermanium[i], fCopyNumber++, fSurfCheck);
   }

   /////////////////////////////////////////////////////////////////////
   // end germaniumBlock1Log
   /////////////////////////////////////////////////////////////////////

   return 1;
}   // end PlaceDeadLayerSpecificCrystal()

///////////////////////////////////////////////////////////////////////
// BuildOneDetector()
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::BuildOneDetector(G4int det)
{
   // Build assembly volumes
   // Holds all pieces that are not a detector (ie. the can, nitrogen tank, cold finger, electrodes, etc.)
   fAssembly = new G4AssemblyVolume();

   // Holds germanium cores
   fGermaniumAssembly = new G4AssemblyVolume();

   // Holds left suppressors
   fLeftSuppressorCasingAssembly = new G4AssemblyVolume();

   // Holds right suppressors
   fRightSuppressorCasingAssembly = new G4AssemblyVolume();

   // holds left suppressor extensions
   fLeftSuppressorExtensionAssembly = new G4AssemblyVolume();

   // holds right suppressor extensions
   fRightSuppressorExtensionAssembly = new G4AssemblyVolume();

   // Holds back suppressors
   fSuppressorBackAssembly = new G4AssemblyVolume();

   // Holds the extension suppressor shells
   fExtensionSuppressorShellAssembly = new G4AssemblyVolume();

   // Holds the back and side suppressor shells
   fBackAndSideSuppressorShellAssembly = new G4AssemblyVolume();

   // Holds the Hevimets
   fHevimetAssembly = new G4AssemblyVolume();

   ConstructComplexDetectorBlockWithDeadLayer();
   BuildelectrodeMatElectrodes();

   ConstructDetector();

   // Include BGOs?
   if(fBGOSelector == 1) {
      ConstructNewSuppressorCasingWithShells(det);
   } else if(fBGOSelector == 0) {
      //G4cout<<"Not building BGO "<<G4endl ;
   } else {
      G4cout << "Error 234235" << G4endl;
      exit(1);
   }

   ConstructColdFinger();

   if(fSuppressorPositionSelector != 0 && fHevimetSelector != 0) {
      ConstructNewHeavyMet();
	}

}   // end BuildOneDetector()

void DetectionSystemGriffin::BuildEverythingButCrystals(G4int det)
{
   // Build assembly volumes
   fAssembly                         = new G4AssemblyVolume();
   fLeftSuppressorCasingAssembly     = new G4AssemblyVolume();
   fRightSuppressorCasingAssembly    = new G4AssemblyVolume();
   fLeftSuppressorExtensionAssembly  = new G4AssemblyVolume();
   fRightSuppressorExtensionAssembly = new G4AssemblyVolume();
   fSuppressorBackAssembly           = new G4AssemblyVolume();

   // Holds the extension suppressor shells
   fExtensionSuppressorShellAssembly = new G4AssemblyVolume();

   // Holds the back and side suppressor shells
   fBackAndSideSuppressorShellAssembly = new G4AssemblyVolume();

   // Holds the hevimets
   fHevimetAssembly = new G4AssemblyVolume();

   BuildelectrodeMatElectrodes();

   ConstructDetector();

   // Include BGOs?
   if(fBGOSelector == 1) {
      ConstructNewSuppressorCasingWithShells(det);
   } else if(fBGOSelector == 0) {
      //G4cout<<"Not building BGO "<<G4endl;
   } else {
      G4cout << "Error 234235" << G4endl;
      exit(1);
   }

   ConstructColdFinger();

   if(fSuppressorPositionSelector != 0 && fHevimetSelector != 0) {
      ConstructNewHeavyMet();
	}

}   // end BuildEverythingButCrystals()

void DetectionSystemGriffin::BuildDeadLayerSpecificCrystal(G4int det)
{
   fGermaniumAssemblyCry[0] = new G4AssemblyVolume();
   fGermaniumAssemblyCry[1] = new G4AssemblyVolume();
   fGermaniumAssemblyCry[2] = new G4AssemblyVolume();
   fGermaniumAssemblyCry[3] = new G4AssemblyVolume();

   for(G4int i = 0; i < 4; i++) {
      ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer(det, i);
   }
}

///////////////////////////////////////////////////////////////////////
// ConstructComplexDetectorBlock builds four quarters of germanium
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructComplexDetectorBlock()
{
   G4Material* materialGe = G4Material::GetMaterial("G4_Ge");
   if(materialGe == nullptr) {
      G4cout << " ----> Material G4_Ge not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   G4Material* materialVacuum = G4Material::GetMaterial("Vacuum");
   if(materialVacuum == nullptr) {
      G4cout << " ----> Material Vacuum not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   auto* germaniumBlock1VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
   germaniumBlock1VisAtt->SetVisibility(true);

   auto* deadLayerVisAtt = new G4VisAttributes(G4Colour(0.80, 0.0, 0.0));
   deadLayerVisAtt->SetVisibility(true);

   auto* airVisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
   airVisAtt->SetVisibility(true);

   auto* electrodeMatVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
   electrodeMatVisAtt->SetVisibility(true);

   G4SubtractionSolid* detector1 = QuarterDetector();

   fGermaniumBlock1Log = new G4LogicalVolume(detector1, materialGe, "germaniumBlock1Log", nullptr, nullptr, nullptr);
   fGermaniumBlock1Log->SetVisAttributes(germaniumBlock1VisAtt);

   fGermaniumAssembly->AddPlacedVolume(fGermaniumBlock1Log, fMoveNull, fRotateNull);

   // now we make the hole that will go in the back of each quarter detector,
   G4double startAngle      = 0.0 * M_PI;
   G4double finalAngle      = 2.0 * M_PI;
   G4double holeRadius      = fGermaniumHoleRadius;
   G4double holeHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace) / 2.0;

   auto* holeTubs = new G4Tubs("holeTubs", 0.0, holeRadius, holeHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveHole(fGermaniumShift, -(fGermaniumShift),
                          -((fGermaniumHoleDistFromFace)));

   fGermaniumHoleLog = new G4LogicalVolume(holeTubs, materialVacuum, "germaniumHoleLog", nullptr, nullptr, nullptr);
   fGermaniumHoleLog->SetVisAttributes(airVisAtt);

   fGermaniumAssembly->AddPlacedVolume(fGermaniumHoleLog, moveHole, fRotateNull);

}   //end ::ConstructComplexDetectorBlock

///////////////////////////////////////////////////////////////////////
// ConstructComplexDetectorBlockWithDeadLayer builds four quarters of
// germanium, with dead layers
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructComplexDetectorBlockWithDeadLayer()
{
   G4Material* materialGe = G4Material::GetMaterial("G4_Ge");
   if(materialGe == nullptr) {
      G4cout << " ----> Material G4_Ge not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   G4Material* materialVacuum = G4Material::GetMaterial("Vacuum");
   if(materialVacuum == nullptr) {
      G4cout << " ----> Material Vacuum not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   auto* germaniumBlock1VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
   germaniumBlock1VisAtt->SetVisibility(true);

   auto* deadLayerVisAtt = new G4VisAttributes(G4Colour(0.80, 0.0, 0.0));
   deadLayerVisAtt->SetVisibility(true);

   auto* airVisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
   airVisAtt->SetVisibility(true);

   auto* electrodeMatVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
   electrodeMatVisAtt->SetVisibility(true);

   G4SubtractionSolid* detector1 = QuarterDetector();

   fGermaniumBlock1Log = new G4LogicalVolume(detector1, materialGe, "germaniumBlock1Log", nullptr, nullptr, nullptr);
   fGermaniumBlock1Log->SetVisAttributes(germaniumBlock1VisAtt);

   fGermaniumAssembly->AddPlacedVolume(fGermaniumBlock1Log, fMoveNull, fRotateNull);

   /////////////////////////////////////////////////////////////////////
   // Since if we are using this method, the inner dead layer must be
   // simulated, make a cylinder of germanium, place the air cylinder
   // inside it, and then place the whole thing in the crystal
   /////////////////////////////////////////////////////////////////////
   G4double startAngle           = 0.0 * M_PI;
   G4double finalAngle           = 2.0 * M_PI;
   G4double deadLayerRadius      = fGermaniumHoleRadius + fInnerDeadLayerThickness;
   G4double deadLayerHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace + fInnerDeadLayerThickness) / 2.0;

   auto* deadLayerTubs = new G4Tubs("deadLayerTubs", fGermaniumHoleRadius,
                                      deadLayerRadius, deadLayerHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveDeadLayer(fGermaniumShift, -(fGermaniumShift),
                               -((fGermaniumHoleDistFromFace - fInnerDeadLayerThickness) / 2.0));

   fInnerDeadLayerLog = new G4LogicalVolume(deadLayerTubs, materialGe, "innerDeadLayerLog", nullptr, nullptr, nullptr);
   fInnerDeadLayerLog->SetVisAttributes(deadLayerVisAtt);

   fGermaniumAssembly->AddPlacedVolume(fInnerDeadLayerLog, moveDeadLayer, fRotateNull);

   // dead layer cap
   deadLayerRadius                  = fGermaniumHoleRadius;
   G4double deadLayerCapHalfLengthZ = (fInnerDeadLayerThickness) / 2.0;

   auto* deadLayerCap = new G4Tubs("deadLayerCap", 0.0,
                                     deadLayerRadius, deadLayerCapHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveDeadLayerCap(fGermaniumShift, -(fGermaniumShift),
                                  -((fGermaniumHoleDistFromFace - fInnerDeadLayerThickness) / 2.0 - deadLayerHalfLengthZ + fInnerDeadLayerThickness / 2.0));

   fInnerDeadLayerCapLog = new G4LogicalVolume(deadLayerCap, materialGe, "innerDeadLayerCapLog", nullptr, nullptr, nullptr);
   fInnerDeadLayerCapLog->SetVisAttributes(deadLayerVisAtt);

   fGermaniumAssembly->AddPlacedVolume(fInnerDeadLayerCapLog, moveDeadLayerCap, fRotateNull);

   // now we make the hole that will go in the back of each quarter detector,
   // and place it inside the dead layer cylinder
   startAngle               = 0.0 * M_PI;
   finalAngle               = 2.0 * M_PI;
   G4double holeRadius      = fGermaniumHoleRadius;
   G4double holeHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace) / 2.0;

   auto* holeTubs = new G4Tubs("holeTubs", 0.0, holeRadius, holeHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveHole(fGermaniumShift, -(fGermaniumShift),
                          -((fGermaniumHoleDistFromFace - fInnerDeadLayerThickness) / 2.0) - (fInnerDeadLayerThickness / 2.0));

   fGermaniumHoleLog = new G4LogicalVolume(holeTubs, materialVacuum, "germaniumHoleLog", nullptr, nullptr, nullptr);
   fGermaniumHoleLog->SetVisAttributes(airVisAtt);

   fGermaniumAssembly->AddPlacedVolume(fGermaniumHoleLog, moveHole, fRotateNull);

}   //end ::ConstructComplexDetectorBlockWithDeadLayer

void DetectionSystemGriffin::ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer(G4int det, G4int cry)
{
   G4String strdet = G4UIcommand::ConvertToString(det);
   G4String strcry = G4UIcommand::ConvertToString(cry);

   G4Material* materialGe = G4Material::GetMaterial("G4_Ge");
   if(materialGe == nullptr) {
      G4cout << " ----> Material G4_Ge not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   G4Material* materialVacuum = G4Material::GetMaterial("Vacuum");
   if(materialVacuum == nullptr) {
      G4cout << " ----> Material Vacuum not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   auto* germaniumBlock1VisAtt = new G4VisAttributes(fGriffinCrystalColours[cry]);
   germaniumBlock1VisAtt->SetVisibility(true);

   auto* deadLayerVisAtt = new G4VisAttributes(fGriffinDeadLayerColours[cry]);
   deadLayerVisAtt->SetVisibility(true);

   auto* airVisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
   airVisAtt->SetVisibility(true);

   auto* electrodeMatVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
   electrodeMatVisAtt->SetVisibility(true);

   G4SubtractionSolid* detector1 = QuarterSpecificDeadLayerDetector(det, cry);

   G4String germaniumBlock1Name = "germaniumDlsBlock1_" + strdet + "_" + strcry + "_log";

   fGermaniumBlock1Log = new G4LogicalVolume(detector1, materialGe, germaniumBlock1Name, nullptr, nullptr, nullptr);
   fGermaniumBlock1Log->SetVisAttributes(germaniumBlock1VisAtt);

   fGermaniumAssemblyCry[cry]->AddPlacedVolume(fGermaniumBlock1Log, fMoveNull, fRotateNull);

   /////////////////////////////////////////////////////////////////////
   // Since if we are using this method, the inner dead layer must be
   // simulated, make a cylinder of germanium, place the air cylinder
   // inside it, and then place the whole thing in the crystal
   /////////////////////////////////////////////////////////////////////
   G4double startAngle           = 0.0 * M_PI;
   G4double finalAngle           = 2.0 * M_PI;
   G4double deadLayerRadius      = fGermaniumHoleRadius + fGriffinDeadLayers[det][cry];
   G4double deadLayerHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace + fGriffinDeadLayers[det][cry]) / 2.0;

   auto* deadLayerTubs = new G4Tubs("deadLayerTubs", fGermaniumHoleRadius,
                                      deadLayerRadius, deadLayerHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveDeadLayer(fGermaniumShift, -(fGermaniumShift),
                               -((fGermaniumHoleDistFromFace - fGriffinDeadLayers[det][cry]) / 2.0));

   fInnerDeadLayerLog = new G4LogicalVolume(deadLayerTubs, materialGe, "innerDeadLayerLog", nullptr, nullptr, nullptr);
   fInnerDeadLayerLog->SetVisAttributes(deadLayerVisAtt);

   fGermaniumAssemblyCry[cry]->AddPlacedVolume(fInnerDeadLayerLog, moveDeadLayer, fRotateNull);

   // dead layer cap
   deadLayerRadius                  = fGermaniumHoleRadius;
   G4double deadLayerCapHalfLengthZ = (fGriffinDeadLayers[det][cry]) / 2.0;

   auto* deadLayerCap = new G4Tubs("deadLayerCap", 0.0,
                                     deadLayerRadius, deadLayerCapHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveDeadLayerCap(fGermaniumShift, -(fGermaniumShift),
                                  -((fGermaniumHoleDistFromFace - fGriffinDeadLayers[det][cry]) / 2.0 - deadLayerHalfLengthZ + fGriffinDeadLayers[det][cry] / 2.0));

   fInnerDeadLayerCapLog = new G4LogicalVolume(deadLayerCap, materialGe, "innerDeadLayerCapLog", nullptr, nullptr, nullptr);
   fInnerDeadLayerCapLog->SetVisAttributes(deadLayerVisAtt);

   fGermaniumAssemblyCry[cry]->AddPlacedVolume(fInnerDeadLayerCapLog, moveDeadLayerCap, fRotateNull);

   // now we make the hole that will go in the back of each quarter detector,
   // and place it inside the dead layer cylinder
   startAngle               = 0.0 * M_PI;
   finalAngle               = 2.0 * M_PI;
   G4double holeRadius      = fGermaniumHoleRadius;
   G4double holeHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace) / 2.0;

   auto* holeTubs = new G4Tubs("holeTubs", 0.0, holeRadius, holeHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveHole(fGermaniumShift, -(fGermaniumShift),
                          -((fGermaniumHoleDistFromFace - fGriffinDeadLayers[det][cry]) / 2.0) - (fGriffinDeadLayers[det][cry] / 2.0));

   fGermaniumHoleLog = new G4LogicalVolume(holeTubs, materialVacuum, "germaniumHoleLog", nullptr, nullptr, nullptr);
   fGermaniumHoleLog->SetVisAttributes(airVisAtt);

   fGermaniumAssemblyCry[cry]->AddPlacedVolume(fGermaniumHoleLog, moveHole, fRotateNull);

}   //end ::ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer

///////////////////////////////////////////////////////////////////////
// Builds a layer of electrodeMat between germanium crystals to
// approximate electrodes, etc. that Eurisys won't tell us
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::BuildelectrodeMatElectrodes()
{

   G4Material* electrodeMat = G4Material::GetMaterial(fElectrodeMaterial);
   if(electrodeMat == nullptr) {
      G4cout << " ----> Electrode material " << fElectrodeMaterial << " not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   auto* electrodeMatVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
   electrodeMatVisAtt->SetVisibility(true);

   // electrodeMat layers between crystals - back part
   auto* rotateElectrodeMatBack = new G4RotationMatrix;
   rotateElectrodeMatBack->rotateZ(M_PI / 2.0);
   G4ThreeVector moveInterCrystalElectrodeMatBack(fGermaniumLength / 2.0 + fGermaniumBentLength / 2.0 + fCanFaceThickness / 2.0 + fGermaniumDistFromCanFace + fShift + fAppliedBackShift, 0, 0);

   G4UnionSolid* interCrystalElectrodeMatBack = InterCrystalelectrodeMatBack();

   fInterCrystalElectrodeMatBackLog = new G4LogicalVolume(interCrystalElectrodeMatBack,
                                                          electrodeMat, "interCrystalElectrodeMatBackLog", nullptr, nullptr, nullptr);
   fInterCrystalElectrodeMatBackLog->SetVisAttributes(electrodeMatVisAtt);

   fAssembly->AddPlacedVolume(fInterCrystalElectrodeMatBackLog, moveInterCrystalElectrodeMatBack, rotateElectrodeMatBack);

   // electrodeMat between crystals - front part
   auto* rotateElectrodeMatFront = new G4RotationMatrix;
   rotateElectrodeMatFront->rotateY(M_PI / 2.0);
   G4UnionSolid* interCrystalElectrodeMatFront = InterCrystalelectrodeMatFront();

   G4ThreeVector moveInterCrystalElectrodeMatFront(fGermaniumBentLength / 2.0 + fElectrodeMatStartingDepth / 2.0 + fCanFaceThickness / 2.0 + fGermaniumDistFromCanFace + fShift + fAppliedBackShift, 0, 0);

   fInterCrystalElectrodeMatFrontLog = new G4LogicalVolume(interCrystalElectrodeMatFront,
                                                           electrodeMat, "interCrystalElectrodeMatFrontLog", nullptr, nullptr, nullptr);
   fInterCrystalElectrodeMatFrontLog->SetVisAttributes(electrodeMatVisAtt);

   fAssembly->AddPlacedVolume(fInterCrystalElectrodeMatFrontLog, moveInterCrystalElectrodeMatFront, rotateElectrodeMatFront);
}   //end ::BuildelectrodeMatElectrodes

///////////////////////////////////////////////////////////////////////
// ConstructDetector builds the structureMat can, cold finger shell,
// and liquid nitrogen tank. Since the can's face is built first,
// everything else is placed relative to it
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructDetector()
{
	std::array<G4RotationMatrix*, 4> rotatePiece;
   std::array<G4ThreeVector, 4>     movePiece;

   G4Material* structureMat = G4Material::GetMaterial(fStructureMaterial);
   if(structureMat == nullptr) {
      G4cout << " ----> Structure material " << fStructureMaterial << " not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   G4Material* liquidN2Material = G4Material::GetMaterial("LiquidN2");
   if(liquidN2Material == nullptr) {
      G4cout << " ----> Material LiquidN2 not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   // first we make the can's front face
   G4Box* frontFace = SquareFrontFace();
   fFrontFaceLog    = new G4LogicalVolume(frontFace, structureMat, "frontFaceLog", nullptr, nullptr, nullptr);

   G4ThreeVector moveFrontFace(fShift + fAppliedBackShift, 0, 0);

   fAssembly->AddPlacedVolume(fFrontFaceLog, moveFrontFace, fRotateNull);

   // now we put on the four angled side pieces

   G4double x0 = ((fBentEndLength) / 2.0) + fShift + fAppliedBackShift;

   G4double z0 = ((fCanFaceThickness - fBentEndLength) * tan(fBentEndAngle) + fDetectorTotalWidth - fCanFaceThickness) / 2.0;

	std::array<G4Para*, 4> sidePiece;

   for(int i = 0; i < 4; i++) {
      // Top, Right, Bottom, Left
      sidePiece[i]   = BentSidePiece();
      rotatePiece[i] = new G4RotationMatrix;

      // The left side is slightly different than the others, so this if statement is required. The order of rotation matters.
      if(i == 3) {
         rotatePiece[i]->rotateX(-M_PI / 2.0 + i * M_PI / 2.0);
         rotatePiece[i]->rotateY(fBentEndAngle);
      } else {
         rotatePiece[i]->rotateY(-fBentEndAngle);
         rotatePiece[i]->rotateX(-M_PI / 2.0 + i * M_PI / 2.0);
      }

      G4double yMov = z0 * cos(i * M_PI / 2.0);
      G4double zMov = z0 * sin(i * M_PI / 2.0);

      movePiece[i] = G4ThreeVector(x0, yMov, zMov);
   }

   fTopBentPieceLog = new G4LogicalVolume(sidePiece[0], structureMat, "topBentPiece", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fTopBentPieceLog, movePiece[0], rotatePiece[0]);

   fRightBentPieceLog = new G4LogicalVolume(sidePiece[1], structureMat, "rightBentPiece", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fRightBentPieceLog, movePiece[1], rotatePiece[1]);

   fBottomBentPieceLog = new G4LogicalVolume(sidePiece[2], structureMat, "bottomBentPiece", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fBottomBentPieceLog, movePiece[2], rotatePiece[2]);

   fLeftBentPieceLog = new G4LogicalVolume(sidePiece[3], structureMat, "leftBentPiece", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fRightBentPieceLog, movePiece[3], rotatePiece[3]);

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // now we add the wedges to the edges of the face. These complete the angled side pieces
	std::array<G4Trap*, 4> sideWedge;

   x0 = fShift + fAppliedBackShift;

   G4double y0 = (fDetectorTotalWidth / 2.0) - (fBentEndLength * tan(fBentEndAngle)) + (fCanFaceThickness / 2.0) * tan(fBentEndAngle) / 2.0;

   for(int i = 0; i < 4; i++) {
      // Top, Right, Bottom, Left
      sideWedge[i]   = CornerWedge();
      rotatePiece[i] = new G4RotationMatrix;

      rotatePiece[i]->rotateX(M_PI / 2.0);
      rotatePiece[i]->rotateY(-M_PI / 2.0);
      rotatePiece[i]->rotateX(-M_PI / 2.0 + i * M_PI / 2.0);

      G4double yMov = y0 * cos(i * M_PI / 2.0);
      G4double zMov = y0 * sin(i * M_PI / 2.0);

      movePiece[i] = G4ThreeVector(x0, yMov, zMov);
   }

   fTopWedgeLog = new G4LogicalVolume(sideWedge[0], structureMat, "topWedgeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fTopWedgeLog, movePiece[0], rotatePiece[0]);

   fRightWedgeLog = new G4LogicalVolume(sideWedge[1], structureMat, "rightWedgeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fRightWedgeLog, movePiece[1], rotatePiece[1]);

   fBottomWedgeLog = new G4LogicalVolume(sideWedge[2], structureMat, "bottomWedgeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fBottomWedgeLog, movePiece[2], rotatePiece[2]);

   fLeftWedgeLog = new G4LogicalVolume(sideWedge[3], structureMat, "leftWedgeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fLeftWedgeLog, movePiece[3], rotatePiece[3]);

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // now we add the rounded corners, made from a quarter of a cone

   // Correction factor used to move the corner cones so as to avoid
   // holes left in the can caused by the bent side pieces

   G4double holeEliminator = fCanFaceThickness * (1.0 - tan(fBentEndAngle));
	std::array<G4Cons*, 4>  coneLocation;

   x0 = fBentEndLength / 2.0 + fShift + fAppliedBackShift;

   y0 = (fDetectorTotalWidth / 2.0) - (fBentEndLength * tan(fBentEndAngle)) - holeEliminator;

   for(int i = 0; i < 4; i++) {
      // Lower Left, Upper Left, Upper Right, Lower Right
      coneLocation[i] = RoundedEndEdge();
      rotatePiece[i]  = new G4RotationMatrix;

      rotatePiece[i]->rotateY(M_PI / 2.0);
      rotatePiece[i]->rotateX(-M_PI / 2.0 + i * M_PI / 2.0);

      G4double yMov = y0 * (-cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      G4double zMov = y0 * (-cos(i * M_PI / 2.0) - sin(i * M_PI / 2.0));

      movePiece[i] = G4ThreeVector(x0, yMov, zMov);
   }

   fLowerLeftConeLog = new G4LogicalVolume(coneLocation[0], structureMat, "lowerLeftConeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fLowerLeftConeLog, movePiece[0], rotatePiece[0]);

   fUpperLeftConeLog = new G4LogicalVolume(coneLocation[1], structureMat, "upperLeftConeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fUpperLeftConeLog, movePiece[1], rotatePiece[1]);

   fUpperRightConeLog = new G4LogicalVolume(coneLocation[2], structureMat, "upperRightConeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fUpperRightConeLog, movePiece[2], rotatePiece[2]);

   fLowerRightConeLog = new G4LogicalVolume(coneLocation[3], structureMat, "lowerRightConeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fLowerRightConeLog, movePiece[3], rotatePiece[3]);

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   /* now we add the corner tubes which extend the rounded corners to the back of the can
	 *
	 * This is extremely similar to the previous loop but they are not exactly the same. While it would be
	 * possible to consolidate them into one, they would only share the rotation variables, so it
	 * might not be a huge benefit.
	 */

	std::array<G4Tubs*, 4> tubeLocation;

   x0 = (fDetectorTotalLength + fBentEndLength - fRearPlateThickness - fCanFaceThickness) / 2.0 + fShift + fAppliedBackShift;

   y0 = (fDetectorTotalWidth / 2.0) - (fBentEndLength * tan(fBentEndAngle));

   for(int i = 0; i < 4; i++) {
      // Lower Left, Upper Left, Upper Right, Lower Right
      tubeLocation[i] = CornerTube();
      rotatePiece[i]  = new G4RotationMatrix;
      rotatePiece[i]->rotateY(M_PI / 2.0);
      rotatePiece[i]->rotateX(-M_PI / 2.0 + i * M_PI / 2.0);

      G4double yMov = y0 * (-cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      G4double zMov = y0 * (-cos(i * M_PI / 2.0) - sin(i * M_PI / 2.0));

      movePiece[i] = G4ThreeVector(x0, yMov, zMov);
   }

   fLowerLeftTubeLog = new G4LogicalVolume(tubeLocation[0], structureMat, "lowerLeftTubeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fLowerLeftTubeLog, movePiece[0], rotatePiece[0]);

   fUpperLeftTubeLog = new G4LogicalVolume(tubeLocation[1], structureMat, "upperLeftTubeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fUpperLeftTubeLog, movePiece[1], rotatePiece[1]);

   fUpperRightTubeLog = new G4LogicalVolume(tubeLocation[2], structureMat, "upperRightTubeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fUpperRightTubeLog, movePiece[2], rotatePiece[2]);

   fLowerRightTubeLog = new G4LogicalVolume(tubeLocation[3], structureMat, "lowerRightTubeLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fLowerRightTubeLog, movePiece[3], rotatePiece[3]);

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // now we add the side panels that extend from the bent pieces to the back of the can

	std::array<G4Box*, 4> panelLocation;

   x0 = (fDetectorTotalLength + fBentEndLength - fRearPlateThickness - fCanFaceThickness) / 2.0 + fShift + fAppliedBackShift;

   y0 = 0;

   z0 = (fDetectorTotalWidth - fCanSideThickness) / 2.0;

   for(int i = 0; i < 4; i++) {
      // Right, Top, Left, Bottom
      panelLocation[i] = SidePanel();
      rotatePiece[i]   = new G4RotationMatrix;
      rotatePiece[i]->rotateX(M_PI / 2.0 - i * M_PI / 2.0);

      G4double yMov = z0 * sin(i * M_PI / 2.0);
      G4double zMov = z0 * cos(i * M_PI / 2.0);

      movePiece[i] = G4ThreeVector(x0, yMov, zMov);
   }

   fRightSidePanelLog = new G4LogicalVolume(panelLocation[0], structureMat, "rightSidePanelLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fRightSidePanelLog, movePiece[0], rotatePiece[0]);

   fTopSidePanelLog = new G4LogicalVolume(panelLocation[1], structureMat, "topSidePanelLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fTopSidePanelLog, movePiece[1], fRotateNull);

   fLeftSidePanelLog = new G4LogicalVolume(panelLocation[2], structureMat, "leftSidePanelLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fLeftSidePanelLog, movePiece[2], rotatePiece[2]);

   fBottomSidePanelLog = new G4LogicalVolume(panelLocation[3], structureMat, "bottomSidePanelLog", nullptr, nullptr, nullptr);
   fAssembly->AddPlacedVolume(fBottomSidePanelLog, movePiece[3], rotatePiece[3]);

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // now we add the rear plate, which has a hole for the cold finger to pass through
   G4SubtractionSolid* rearPlate = RearPlate();

   auto* rotateRearPlate = new G4RotationMatrix;
   rotateRearPlate->rotateY(M_PI / 2.0);

   x0 = fDetectorTotalLength - fCanFaceThickness / 2.0 - fRearPlateThickness / 2.0 + fShift + fAppliedBackShift;

   G4ThreeVector moveRearPlate(fDetectorTotalLength - fCanFaceThickness / 2.0 - fRearPlateThickness / 2.0 + fShift + fAppliedBackShift, 0, 0);

   fRearPlateLog = new G4LogicalVolume(rearPlate, structureMat, "rearPlateLog", nullptr, nullptr, nullptr);

   fAssembly->AddPlacedVolume(fRearPlateLog, moveRearPlate, rotateRearPlate);

   // we know add the cold finger shell, which extends out the back of the can
   G4Tubs* fingerShell = ColdFingerShell();

   auto* rotateFingerShell = new G4RotationMatrix;
   rotateFingerShell->rotateY(M_PI / 2.0);

   G4ThreeVector moveFingerShell(fColdFingerShellLength / 2.0 - fCanFaceThickness / 2.0 + fDetectorTotalLength + fShift + fAppliedBackShift, 0, 0);

   fFingerShellLog = new G4LogicalVolume(fingerShell, structureMat, "fingerShellLog", nullptr, nullptr, nullptr);

   fAssembly->AddPlacedVolume(fFingerShellLog, moveFingerShell, rotateFingerShell);

   // lastly we add the liquid nitrogen tank at the back, past the cold finger shell
   G4Tubs* tank = LiquidNitrogenTank();

   auto* rotateTank = new G4RotationMatrix;
   rotateTank->rotateY(M_PI / 2.0);

   G4ThreeVector moveTank((fCoolantLength - fCanFaceThickness) / 2.0 + fDetectorTotalLength + fColdFingerShellLength + fShift + fAppliedBackShift, 0, 0);

   fTankLog = new G4LogicalVolume(tank, structureMat, "tankLog", nullptr, nullptr, nullptr);

   fAssembly->AddPlacedVolume(fTankLog, moveTank, rotateTank);

   // lids
   G4Tubs* lid = LiquidNitrogenTankLid();

   auto* rotateLid = new G4RotationMatrix;
   rotateLid->rotateY(M_PI / 2.0);

   G4ThreeVector moveLid1(fCoolantLength / 2.0 - fCoolantThickness / 2.0 + (fCoolantLength - fCanFaceThickness) / 2.0 + fDetectorTotalLength + fColdFingerShellLength + fShift + fAppliedBackShift, 0, 0);
   G4ThreeVector moveLid2(-1.0 * fCoolantLength / 2.0 + fCoolantThickness / 2.0 + (fCoolantLength - fCanFaceThickness) / 2.0 + fDetectorTotalLength + fColdFingerShellLength + fShift + fAppliedBackShift, 0, 0);

   fTankLid1Log = new G4LogicalVolume(lid, structureMat, "tankLid1Log", nullptr, nullptr, nullptr);
   fTankLid2Log = new G4LogicalVolume(lid, structureMat, "tankLid2Log", nullptr, nullptr, nullptr);

   fAssembly->AddPlacedVolume(fTankLid1Log, moveLid1, rotateLid);
   fAssembly->AddPlacedVolume(fTankLid2Log, moveLid2, rotateLid);

   // LN2
   G4Tubs* liquidn2 = LiquidNitrogen();

   auto* rotateLiquidn2 = new G4RotationMatrix;
   rotateLiquidn2->rotateY(M_PI / 2.0);

   G4ThreeVector moveLiquidn2((fCoolantLength - fCanFaceThickness) / 2.0 + fDetectorTotalLength + fColdFingerShellLength + fShift + fAppliedBackShift, 0, 0);

   fTankLiquidLog = new G4LogicalVolume(liquidn2, liquidN2Material, "tankLiquidLog", nullptr, nullptr, nullptr);

   fAssembly->AddPlacedVolume(fTankLiquidLog, moveLiquidn2, rotateLiquidn2);
}   // end::ConstructDetector

///////////////////////////////////////////////////////////////////////
// ConstructColdFinger has changed in this version!  It now
// incorporates all the internal details released in Nov 2004 by
// Eurisys. ConstructColdFinger builds the cold finger as well as the
// cold plate. The finger extends to the Liquid Nitrogen tank, while
// the plate is always the same distance from the back of the germanium
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructColdFinger()
{
   G4Material* materialAir = G4Material::GetMaterial("Air");
   if(materialAir == nullptr) {
      G4cout << " ----> Material Air not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   G4Material* structureMat = G4Material::GetMaterial(fStructureMaterial);
   if(structureMat == nullptr) {
      G4cout << " ----> Structure material " << fStructureMaterial << " not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   G4Material* electrodeMat = G4Material::GetMaterial(fElectrodeMaterial);
   if(electrodeMat == nullptr) {
      G4cout << " ----> Electrode material " << fElectrodeMaterial << " not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   auto* coldFingerVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
   coldFingerVisAtt->SetVisibility(true);

   auto* structureMatVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
   structureMatVisAtt->SetVisibility(true);

   auto* airVisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
   airVisAtt->SetVisibility(true);

   G4Box* endPlate = EndPlate();

   // cut out holes
   G4double airHoleDistance = fGermaniumSeparation / 2.0 + (fGermaniumWidth / 2.0 - fGermaniumShift);   //centre of the middle hole

   auto* rotateFetAirHole = new G4RotationMatrix;
   rotateFetAirHole->rotateY(M_PI / 2.0);

   G4Tubs* fetAirHoleCut = AirHoleCut();

	std::array<G4ThreeVector, 8>       movePiece;
   std::array<G4SubtractionSolid*, 8> endPlateCut;

   G4double y0 = airHoleDistance;
   G4double z0 = airHoleDistance;

   for(int i = 0; i < 4; i++) {
      G4double y            = y0 * (-cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      G4double z            = z0 * (cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      movePiece[i] = G4ThreeVector(0, y, z);

      if(i == 0) {
         // Slightly different than the rest. Sorry I couldnt think of a better workaround
         endPlateCut[i] = new G4SubtractionSolid("endPlateCut", endPlate, fetAirHoleCut, rotateFetAirHole, movePiece[i]);
		} else {
         endPlateCut[i] = new G4SubtractionSolid("endPlateCut", endPlateCut[i - 1], fetAirHoleCut, rotateFetAirHole, movePiece[i]);
		}
   }

   G4double x = fColdFingerEndPlateThickness / 2.0 + fCanFaceThickness / 2.0 + fGermaniumDistFromCanFace + fGermaniumLength + fColdFingerSpace + fShift + fAppliedBackShift;

   G4ThreeVector moveEndPlate(x, 0, 0);

   // Fix overlap.

	std::array<G4RotationMatrix*, 4> cutRotatePiece;
   std::array<G4ThreeVector, 4>     cutMovePiece;
   std::array<G4Tubs*, 4>           cutTubeLocation;

   G4double cutTubeY0 = (fDetectorTotalWidth / 2.0) - (fBentEndLength * tan(fBentEndAngle));

   for(int i = 0; i < 4; i++) {
      // Lower Left, Upper Left, Upper Right, Lower Right
      cutTubeLocation[i] = CornerTube();
      cutRotatePiece[i]  = new G4RotationMatrix;
      cutRotatePiece[i]->rotateY(M_PI / 2.0);
      cutRotatePiece[i]->rotateZ(-M_PI / 2.0 + (i - 1) * M_PI / 2.0);

      //cutRotatePiece[i]->rotateX(-M_PI/2.0 + i*M_PI/2.0) ;

      G4double cutYMov = cutTubeY0 * (-cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      G4double cutZMov = cutTubeY0 * (-cos(i * M_PI / 2.0) - sin(i * M_PI / 2.0));

      cutMovePiece[i] = G4ThreeVector(0, cutYMov, cutZMov);

      endPlateCut[i + 4] = new G4SubtractionSolid("endPlateCut", endPlateCut[i + 3], cutTubeLocation[i], cutRotatePiece[i], cutMovePiece[i]);
   }

   fEndPlateLog = new G4LogicalVolume(endPlateCut[7], structureMat, "endPlateLog", nullptr, nullptr, nullptr);

   fEndPlateLog->SetVisAttributes(structureMatVisAtt);

   fAssembly->AddPlacedVolume(fEndPlateLog, moveEndPlate, nullptr);

   // Place holes in the old cold plate
   G4Tubs* fetAirHole = AirHole();

   fFetAirHoleLog = new G4LogicalVolume(fetAirHole, materialAir, "fetAirHoleLog", nullptr, nullptr, nullptr);
   fFetAirHoleLog->SetVisAttributes(airVisAtt);

   for(int i = 0; i < 4; i++) {
      G4double y   = y0 * (-cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      G4double z   = z0 * (cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      movePiece[i] = G4ThreeVector(x, y, z);
      fAssembly->AddPlacedVolume(fFetAirHoleLog, movePiece[i], rotateFetAirHole);
   }

   // Cold Finger
   G4Tubs* fing = Finger();

   auto* rotateFinger = new G4RotationMatrix;
   rotateFinger->rotateY(M_PI / 2.0);

   G4ThreeVector moveFinger((fColdFingerLength + fCanFaceThickness) / 2.0 + fGermaniumDistFromCanFace + fGermaniumLength + fColdFingerSpace + fColdFingerEndPlateThickness + fShift + fAppliedBackShift + fStructureMatColdFingerThickness / 2.0, 0, 0);   //changed Jan 2005

   fFingerLog = new G4LogicalVolume(fing, electrodeMat, "fingerLog", nullptr, nullptr, nullptr);
   fFingerLog->SetVisAttributes(coldFingerVisAtt);

   fAssembly->AddPlacedVolume(fFingerLog, moveFinger, rotateFinger);

   // The extra block of structureMat that's in the diagram "Cut C"
   G4SubtractionSolid* extraColdBlock       = ExtraColdBlock();
   auto*   rotateExtraColdBlock = new G4RotationMatrix;
   rotateExtraColdBlock->rotateY(M_PI / 2.0);

   G4ThreeVector moveExtraColdBlock(fDetectorTotalLength - fCanFaceThickness / 2.0 - fRearPlateThickness + fShift + fAppliedBackShift - fExtraBlockDistanceFromBackPlate - fExtraBlockThickness / 2.0, 0, 0);

   fExtraColdBlockLog = new G4LogicalVolume(extraColdBlock, structureMat, "extraColdBlockLog", nullptr, nullptr, nullptr);

   fAssembly->AddPlacedVolume(fExtraColdBlockLog, moveExtraColdBlock, rotateExtraColdBlock);

   // The structures above the cooling end plate
   G4Tubs* structureMatColdFinger = StructureMatColdFinger();
   fStructureMatColdFingerLog     = new G4LogicalVolume(structureMatColdFinger, structureMat, "structureMatColdFingerLog", nullptr, nullptr, nullptr);
   fStructureMatColdFingerLog->SetVisAttributes(structureMatVisAtt);

   auto* rotateStructureMatColdFinger = new G4RotationMatrix;
   rotateStructureMatColdFinger->rotateY(M_PI / 2.0);
   G4ThreeVector moveStructureMatColdFinger(fStructureMatColdFingerThickness / 2.0 + fColdFingerEndPlateThickness + fCanFaceThickness / 2.0 + fGermaniumDistFromCanFace + fGermaniumLength + fColdFingerSpace + fShift + fAppliedBackShift, 0, 0);

   fAssembly->AddPlacedVolume(fStructureMatColdFingerLog, moveStructureMatColdFinger, rotateStructureMatColdFinger);

   // Cooling Side Block
   G4Box* coolingSideBlock = CoolingSideBlock();
   fCoolingSideBlockLog    = new G4LogicalVolume(coolingSideBlock, structureMat, "coolingSideBlockLog", nullptr, nullptr, nullptr);
   fCoolingSideBlockLog->SetVisAttributes(structureMatVisAtt);

   G4Box* coolingBar = CoolingBar();
   fCoolingBarLog    = new G4LogicalVolume(coolingBar, structureMat, "coolingBarLog", nullptr, nullptr, nullptr);
   fCoolingBarLog->SetVisAttributes(structureMatVisAtt);

	std::array<G4RotationMatrix*, 8> rotatePiece;

   x = fCoolingSideBlockThickness / 2.0 + fColdFingerEndPlateThickness + fCanFaceThickness / 2.0 + fGermaniumDistFromCanFace + fGermaniumLength + fColdFingerSpace + fShift + fAppliedBackShift;

   G4double y01 = fGermaniumWidth - fCoolingSideBlockHorizontalDepth / 2.0;

   G4double y02 = fGermaniumWidth - fCoolingSideBlockHorizontalDepth - (fGermaniumWidth - fCoolingSideBlockHorizontalDepth - fStructureMatColdFingerRadius) / 2.0;

   for(int i = 0; i < 4; i++) {
      // Add cooling side blocks and cooling bars at the same time.
      rotatePiece[2 * i] = new G4RotationMatrix;   // cooling side blocks
      rotatePiece[2 * i]->rotateZ(M_PI / 2.0);
      rotatePiece[2 * i]->rotateX(M_PI / 2.0 - i * M_PI / 2.0);

      G4double y = -y01 * cos(i * M_PI / 2.0);
      G4double z = -y01 * sin(i * M_PI / 2.0);

      movePiece[2 * i] = G4ThreeVector(x, y, z);

      rotatePiece[2 * i + 1] = new G4RotationMatrix;   // cooling bars
      rotatePiece[2 * i + 1]->rotateZ(M_PI / 2.0);
      rotatePiece[2 * i + 1]->rotateX(M_PI / 2.0 - i * M_PI / 2.0);

      y = -y02 * cos(i * M_PI / 2.0);
      z = -y02 * sin(i * M_PI / 2.0);

      movePiece[2 * i + 1] = G4ThreeVector(x, y, z);

      fAssembly->AddPlacedVolume(fCoolingSideBlockLog, movePiece[2 * i], rotatePiece[2 * i]);
      fAssembly->AddPlacedVolume(fCoolingBarLog, movePiece[2 * i + 1], rotatePiece[2 * i + 1]);
   }

   // Triangle posts from "Cut A"
   // First, find how far from the centre to place the tips of the triangles

   G4double distanceOfTheTip = fGermaniumSeparation / 2.0 + (fGermaniumWidth / 2.0 - fGermaniumShift)   //centre of the middle hole
                               + sqrt(pow((fGermaniumOuterRadius + fTrianglePostsDistanceFromCrystals), 2.0) - pow(fGermaniumWidth / 2.0 - fGermaniumShift, 2.0));

   // The distance of the base from the detector centre
   G4double distanceOfTheBase = fGermaniumSeparation / 2.0 + fGermaniumWidth;

   G4double trianglePostLength = fGermaniumLength - fTrianglePostStartingDepth + fColdFingerSpace;

   G4Trd* trianglePost = TrianglePost();
   fTrianglePostLog    = new G4LogicalVolume(trianglePost, structureMat, "trianglePostLog", nullptr, nullptr, nullptr);

   x = fCanFaceThickness / 2.0 + fGermaniumDistFromCanFace + fGermaniumLength + fColdFingerSpace + fShift + fAppliedBackShift - trianglePostLength / 2.0;

   y0 = (distanceOfTheBase + distanceOfTheTip) / 2.0;

   for(int i = 0; i < 4; i++) {
      rotatePiece[i] = new G4RotationMatrix;
      rotatePiece[i]->rotateY(0);                             // Initially included...seems a little redundant.
      rotatePiece[i]->rotateX(M_PI / 2.0 - i * M_PI / 2.0);   // 4, 1, 2, 3

      G4double y = -y0 * cos(i * M_PI / 2.0);
      G4double z = y0 * sin(i * M_PI / 2.0);

      movePiece[i] = G4ThreeVector(x, y, z);

      fAssembly->AddPlacedVolume(fTrianglePostLog, movePiece[i], rotatePiece[i]);
   }
}   //end ::ConstructColdFinger

///////////////////////////////////////////////////////////////////////
// methods used in ConstructColdFinger()
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGriffin::EndPlate() const
{
   G4double halfThicknessX = fColdFingerEndPlateThickness / 2.0;
   G4double halfLengthY    = fGermaniumWidth;   //Since it must cover the back of the detector
   G4double halfLengthZ    = halfLengthY;       //Since it is symmetric

   return new G4Box("endPlate", halfThicknessX, halfLengthY, halfLengthZ);
}   //end ::endPlate

G4Tubs* DetectionSystemGriffin::Finger() const
{
   G4double innerRadius = 0.0 * cm;
   G4double outerRadius = fColdFingerRadius;
   G4double halfLengthZ = (fColdFingerLength - fStructureMatColdFingerThickness) / 2.0;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("finger", innerRadius, outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::finger

G4SubtractionSolid* DetectionSystemGriffin::ExtraColdBlock()
{
   G4double halfWidthX     = fGermaniumWidth;
   G4double halfWidthY     = halfWidthX;
   G4double halfThicknessZ = fExtraBlockThickness / 2.0;

   auto* plate = new G4Box("plate", halfWidthX, halfWidthY, halfThicknessZ);

   G4double innerRadius = 0.0;
   G4double outerRadius = fExtraBlockInnerDiameter / 2.0;

   G4double halfHeightZ = halfThicknessZ + 1.0 * cm;   //+1.0*cm just to make sure the hole goes completely through the plate
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   auto* hole = new G4Tubs("hole", innerRadius, outerRadius, halfHeightZ, startAngle, finalAngle);

   auto* extraColdBlock = new G4SubtractionSolid("extraColdBlock", plate, hole);

   // Fix overlap
	std::array<G4SubtractionSolid*, 4> fixExtraColdBlock;
   std::array<G4RotationMatrix*, 4>   cutRotatePiece;
   std::array<G4ThreeVector, 4>       cutMovePiece;
   std::array<G4Tubs*, 4>             cutTubeLocation;

   G4double cutTubeY0 = (fDetectorTotalWidth / 2.0) - (fBentEndLength * tan(fBentEndAngle));

   for(G4int i = 0; i < 4; i++) {
      // Lower Left, Upper Left, Upper Right, Lower Right
      cutTubeLocation[i] = CornerTube();
      cutRotatePiece[i]  = new G4RotationMatrix;
      cutRotatePiece[i]->rotateZ(-M_PI / 2.0 + (i - 1) * M_PI / 2.0);

      G4double cutYMov = cutTubeY0 * (-cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0));
      G4double cutZMov = cutTubeY0 * (-cos(i * M_PI / 2.0) - sin(i * M_PI / 2.0));

      cutMovePiece[i] = G4ThreeVector(cutZMov, cutYMov, 0);

      if(i == 0) {
         fixExtraColdBlock[i] = new G4SubtractionSolid("extraColdBlock", extraColdBlock, cutTubeLocation[i], cutRotatePiece[i], cutMovePiece[i]);
      } else {
         fixExtraColdBlock[i] = new G4SubtractionSolid("extraColdBlock", fixExtraColdBlock[i - 1], cutTubeLocation[i], cutRotatePiece[i], cutMovePiece[i]);
      }
   }

   return fixExtraColdBlock[3];
}   //end ::extraColdBlock

G4Trd* DetectionSystemGriffin::TrianglePost() const
{
   // Calculations that are also done in ConstructColdFinger for positioning
   // First, find how far from the centre to place the tips of the triangles
   G4double distanceOfTheTip = fGermaniumSeparation / 2.0 + (fGermaniumWidth / 2.0 - fGermaniumShift)   //centre of the middle hole
                               + sqrt(pow((fGermaniumOuterRadius + fTrianglePostsDistanceFromCrystals), 2.0) - pow(fGermaniumWidth / 2.0 - fGermaniumShift, 2.0));

   // The distance of the base from the detector centre
   G4double distanceOfTheBase = fGermaniumSeparation / 2.0 + fGermaniumWidth;

   // The distance away from the boundary between crystals of the side points
   G4double distanceOfTheSidePoints = sqrt(pow((fGermaniumOuterRadius + fTrianglePostsDistanceFromCrystals), 2.0) - pow(fGermaniumWidth / 2.0 + /*notice*/ fGermaniumShift, 2.0));

   // Measurements to make the posts with
   G4double length = fGermaniumLength - fTrianglePostStartingDepth + fColdFingerSpace;

   G4double baseToTipHeight = (distanceOfTheBase - distanceOfTheTip);

   G4double halfWidthOfBase = distanceOfTheSidePoints;

   G4double halfWidthOfTop = fTrianglePostDim;   //the easiest way to make a triangle
   //G4Trd(const G4String& pName,G4double  dx1, G4double dx2, G4double  dy1, G4double dy2,G4double  dz)
   return new G4Trd("trianglePost", length / 2.0, length / 2.0,
                    halfWidthOfTop, halfWidthOfBase, baseToTipHeight / 2.0);
}   //end ::trianglePost

G4Tubs* DetectionSystemGriffin::AirHole() const
{
   G4double innerRadius = 0.0 * cm;
   G4double outerRadius = fFetAirHoleRadius;
   G4double halfLengthZ = fColdFingerEndPlateThickness / 2.0;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("fetAirHole", innerRadius,
                     outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::airHole

G4Tubs* DetectionSystemGriffin::AirHoleCut() const
{
   G4double innerRadius = 0.0 * cm;
   G4double outerRadius = fFetAirHoleRadius;
   G4double halfLengthZ = fColdFingerEndPlateThickness / 2.0 + 1.0 * cm;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("fetAirHole", innerRadius,
                     outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::airHoleCut

G4Box* DetectionSystemGriffin::CoolingBar() const
{
   G4double halfThicknessX = fCoolingBarThickness / 2.0;
   G4double halfLengthY    = fCoolingBarWidth / 2.0;
   G4double halfLengthZ    = (fGermaniumWidth - fCoolingSideBlockHorizontalDepth - fStructureMatColdFingerRadius) / 2.0;

   return new G4Box("coolingBar", halfThicknessX, halfLengthY, halfLengthZ);
}   //end ::coolingBar

G4Box* DetectionSystemGriffin::CoolingSideBlock() const
{
   G4double halfThicknessX = fCoolingSideBlockWidth / 2.0;
   G4double halfLengthY    = fCoolingSideBlockThickness / 2.0;
   G4double halfLengthZ    = fCoolingSideBlockHorizontalDepth / 2.0;

   return new G4Box("coolingSideBlock",
                    halfThicknessX, halfLengthY, halfLengthZ);
}   //end ::coolingSideBlock

G4Tubs* DetectionSystemGriffin::StructureMatColdFinger() const
{
   G4double innerRadius = 0.0 * cm;
   G4double outerRadius = fStructureMatColdFingerRadius;
   G4double halfLengthZ = fStructureMatColdFingerThickness / 2.0;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("structureMatColdFinger",
                     innerRadius, outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::structureMatColdFinger

///////////////////////////////////////////////////////////////////////
// ConstructNewSuppressorCasingWithShells builds the suppressor that
// surrounds the can. It tapers off at the front, and there is a thick
// piece covering the back of the can, which is divided in the middle,
// as per the design in the specifications. Also, there is a layer
// of structureMat that surrounds the physical pieces.
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructNewSuppressorCasingWithShells(G4int det)
{
   G4String strdet = G4UIcommand::ConvertToString(det);

   G4Material* structureMat = G4Material::GetMaterial(fStructureMaterial);
   if(structureMat == nullptr) {
      G4cout << " ----> Structure material " << fStructureMaterial << " not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }
   G4Material* materialBGO = G4Material::GetMaterial(fBGOMaterial);
   if(materialBGO == nullptr) {
      G4cout << " ----> BGO material " << fBGOMaterial << " not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   // Change some values to accomodate the shells
   // Replacement for sideSuppressorLength
   G4double shellSideSuppressorLength = fSideSuppressorLength + (fSuppressorShellThickness * 2.0);

   // Replacement for suppressorExtensionLength
   G4double shellSuppressorExtensionLength = fSuppressorExtensionLength + (fSuppressorShellThickness * 2.0) * (1.0 / tan(fBentEndAngle) - tan(fBentEndAngle));

   // Replacement for suppressorExtensionAngle: must totally recalculate
   G4double shellSuppressorExtensionAngle = atan(((fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0) / tan(fBentEndAngle) - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * sin(fBentEndAngle)) * tan(fBentEndAngle) - (fSuppressorForwardRadius + fHevimetTipThickness) * sin(fBentEndAngle)) / (shellSuppressorExtensionLength));

   auto* SuppressorVisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
   SuppressorVisAtt->SetVisibility(true);
   auto* innardsVisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
   innardsVisAtt->SetVisibility(true);
   auto* backInnardsVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
   backInnardsVisAtt->SetVisibility(true);

   // first we add the four pieces of back suppressor, and their shells

   G4SubtractionSolid* backQuarterSuppressor = BackSuppressorQuarter();

   G4SubtractionSolid* backQuarterSuppressorShell = ShellForBackSuppressorQuarter();

   // Insert the structureMat shell first.  The shells must be given numbers, rather than
   // the suppressor pieces themselves.  As an error checking method, the suppressor
   // pieces are given a copy number value far out of range of any useful copy number.
   if(fIncludeBackSuppressors) {
      G4String backQuarterSuppressorShellName = "backQuarterSuppressor_" + strdet + "_log";

      fBackQuarterSuppressorShellLog = new G4LogicalVolume(
         backQuarterSuppressorShell, structureMat,
         backQuarterSuppressorShellName, nullptr, nullptr, nullptr);

      fBackQuarterSuppressorShellLog->SetVisAttributes(SuppressorVisAtt);

      G4Material* backMaterial = G4Material::GetMaterial(fBackSuppressorMaterial);
      if(backMaterial == nullptr) {
         G4cout << " ----> Back suppressor material " << fBackSuppressorMaterial << " not found, cannot build the detector shell! " << G4endl;
         exit(1);
      }

      G4String backQuarterSuppressorName = "backQuarterSuppressor_" + strdet + "_log";

      fBackQuarterSuppressorLog = new G4LogicalVolume(backQuarterSuppressor, backMaterial,
                                                      backQuarterSuppressorName, nullptr, nullptr, nullptr);
      fBackQuarterSuppressorLog->SetVisAttributes(backInnardsVisAtt);

      fSuppressorBackAssembly->AddPlacedVolume(fBackQuarterSuppressorLog, fMoveNull, fRotateNull);

      G4double x0 = (fBackBGOThickness - fCanFaceThickness) / 2.0 + fSuppressorShellThickness + fDetectorTotalLength + fBGOCanSeperation + fShift + fAppliedBackShift;
      G4double y0 = fDetectorTotalWidth / 4.0;
      G4double z0 = fDetectorTotalWidth / 4.0;

		std::array<G4RotationMatrix*, 4> rotateBackSuppressorShells;
      std::array<G4ThreeVector, 4>     moveBackQuarterSuppressor;

      for(int i = 0; i < 4; i++) {
         rotateBackSuppressorShells[i] = new G4RotationMatrix;
         rotateBackSuppressorShells[i]->rotateX(-M_PI / 2.0 + i * M_PI / 2.0);

         G4double x = x0;
         G4double y = y0 * (sin(i * M_PI / 2.0) - cos(i * M_PI / 2.0));
         G4double z = -z0 * (sin(i * M_PI / 2.0) + cos(i * M_PI / 2.0));

         moveBackQuarterSuppressor[i] = G4ThreeVector(x, y, z);

         fBackAndSideSuppressorShellAssembly->AddPlacedVolume(fBackQuarterSuppressorShellLog, moveBackQuarterSuppressor[i], rotateBackSuppressorShells[i]);
      }
   }

   ////////////////////////////////////////////////////////////////////////////////////////
   // now we add the side pieces of suppressor that taper off towards the front of the can
   ////////////////////////////////////////////////////////////////////////////////////////
   G4String rightSuppressorShellLogName  = "rightSuppressorShell_" + strdet + "_log";
   G4String leftSuppressorShellLogName   = "leftSuppressorShell_" + strdet + "_log";
   G4String rightSuppressorCasingLogName = "rightSuppressorCasing_" + strdet + "_log";
   G4String leftSuppressorCasingLogName  = "leftSuppressorCasing_" + strdet + "_log";

   // Define the structureMat shell logical volume
   G4SubtractionSolid* rightSuppressorShell = ShellForFrontSlantSuppressor("right");

   fRightSuppressorShellLog = new G4LogicalVolume(rightSuppressorShell, structureMat,
                                                  rightSuppressorShellLogName, nullptr, nullptr, nullptr);
   fRightSuppressorShellLog->SetVisAttributes(SuppressorVisAtt);

   G4SubtractionSolid* leftSuppressorShell = ShellForFrontSlantSuppressor("left");

   fLeftSuppressorShellLog = new G4LogicalVolume(leftSuppressorShell, structureMat,
                                                 leftSuppressorShellLogName, nullptr, nullptr, nullptr);
   fLeftSuppressorShellLog->SetVisAttributes(SuppressorVisAtt);

   G4SubtractionSolid* rightSuppressor = FrontSlantSuppressor("right", false);   // Right, non-chopping.

   fRightSuppressorLog = new G4LogicalVolume(rightSuppressor, materialBGO, rightSuppressorCasingLogName, nullptr, nullptr, nullptr);
   fRightSuppressorLog->SetVisAttributes(innardsVisAtt);

   G4SubtractionSolid* leftSuppressor = FrontSlantSuppressor("left", false);   // Left, non-chopping.

   fLeftSuppressorLog = new G4LogicalVolume(leftSuppressor, materialBGO, leftSuppressorCasingLogName, nullptr, nullptr, nullptr);
   fLeftSuppressorLog->SetVisAttributes(innardsVisAtt);

   fRightSuppressorCasingAssembly->AddPlacedVolume(fRightSuppressorLog, fMoveNull, fRotateNull);
   fLeftSuppressorCasingAssembly->AddPlacedVolume(fLeftSuppressorLog, fMoveNull, fRotateNull);

   /////////////////////////////////////////////////////////////////////
   // Note : Left and Right are read from the BACK of the detector
   // Suppressors 1 and 2 cover germanium 1
   // Suppressors 3 and 4 cover germanium 2
   // Suppressors 5 and 6 cover germanium 3
   // Suppressors 7 and 8 cover germanium 4
   /////////////////////////////////////////////////////////////////////

   G4double x0 = shellSideSuppressorLength / 2.0 - fCanFaceThickness / 2.0 + fBentEndLength + (fBGOCanSeperation + fBGOChoppedTip) / tan(fBentEndAngle) + fShift + fAppliedBackShift;
   G4double y0 = fSideBGOThickness / 2.0 + fSuppressorShellThickness + fDetectorTotalWidth / 2.0 + fBGOCanSeperation;
   G4double z0 = fSideBGOThickness / 2.0 + fSuppressorShellThickness + fDetectorTotalWidth / 2.0 + fBGOCanSeperation;

	std::array<G4RotationMatrix*, 8> rotateSuppressorExtensionShell;
   std::array<G4ThreeVector, 8>     moveSuppressorExtensionShell;

   for(int i = 0; i < 4; i++) {
      //********************** RIGHT **********************//
      rotateSuppressorExtensionShell[2 * i] = new G4RotationMatrix;
      rotateSuppressorExtensionShell[2 * i]->rotateZ(M_PI / 2.0);
      rotateSuppressorExtensionShell[2 * i]->rotateY(M_PI / 2.0);
      rotateSuppressorExtensionShell[2 * i]->rotateX(M_PI * (1 - i / 2.0));

      G4double x = x0;
      G4double y = y0 * (-cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0)) / (1 + (i + 1) % 2);   // -y/2, y, y/2, -y
      G4double z = z0 * (cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0)) / (1 + i % 2);   // z, z/2, -z, -z/2

      moveSuppressorExtensionShell[2 * i] = G4ThreeVector(x, y, z);

      fBackAndSideSuppressorShellAssembly->AddPlacedVolume(fRightSuppressorShellLog, moveSuppressorExtensionShell[2 * i], rotateSuppressorExtensionShell[2 * i]);

      //*********************** LEFT ***********************//
      rotateSuppressorExtensionShell[2 * i + 1] = new G4RotationMatrix;
      rotateSuppressorExtensionShell[2 * i + 1]->rotateY(-M_PI / 2.0);
      rotateSuppressorExtensionShell[2 * i + 1]->rotateX(M_PI * (1 - i / 2.0));

      x = x0;
      y = y0 * (cos(i * M_PI / 2.0) - sin(i * M_PI / 2.0)) / (1 + i % 2);   // y, -y/2, -y, y/2
      z = -z0 * (cos(i * M_PI / 2.0) + sin(i * M_PI / 2.0)) / (1 + (i + 1) % 2);   // -z/2, -z, z/2, z

      moveSuppressorExtensionShell[2 * i + 1] = G4ThreeVector(x, y, z);

      fBackAndSideSuppressorShellAssembly->AddPlacedVolume(fLeftSuppressorShellLog, moveSuppressorExtensionShell[2 * i + 1], rotateSuppressorExtensionShell[2 * i + 1]);
   }

   ////////////////////////////////////////////////////////////////////////////////////////
   // now we add the side pieces of suppressor that extend out in front of the can when it's in the back position
   ////////////////////////////////////////////////////////////////////////////////////////
   G4String rightSuppressorExtensionLogName = "rightSuppressorExtension_" + strdet + "_log";
   G4String leftSuppressorExtensionLogName  = "leftSuppressorExtension_" + strdet + "_log";

   G4String rightSuppressorShellExtensionLogName = "rightSuppressorShellExtension_" + strdet + "_log";
   G4String leftSuppressorShellExtensionLogName  = "leftSuppressorShellExtension_" + strdet + "_log";

   // Define the shell right logical volume
   // G4SubtractionSolid* rightSuppressorShellExtension = ShellForRightSuppressorExtension();
   G4SubtractionSolid* rightSuppressorShellExtension = ShellForSuppressorExtension("right");

   fRightSuppressorShellExtensionLog = new G4LogicalVolume(rightSuppressorShellExtension,
                                                           materialBGO, rightSuppressorShellExtensionLogName, nullptr, nullptr, nullptr);
   fRightSuppressorShellExtensionLog->SetVisAttributes(SuppressorVisAtt);

   G4SubtractionSolid* rightSuppressorExtension = SideSuppressorExtension("right", false);   // Right, non-chopping // CALLED

   fRightSuppressorExtensionLog = new G4LogicalVolume(rightSuppressorExtension, materialBGO,
                                                      rightSuppressorExtensionLogName, nullptr, nullptr, nullptr);
   fRightSuppressorExtensionLog->SetVisAttributes(innardsVisAtt);

   // Define the left shell logical volume
   // G4SubtractionSolid* leftSuppressorShellExtension = ShellForLeftSuppressorExtension();
   G4SubtractionSolid* leftSuppressorShellExtension = ShellForSuppressorExtension("left");

   fLeftSuppressorShellExtensionLog = new G4LogicalVolume(leftSuppressorShellExtension,
                                                          materialBGO, leftSuppressorShellExtensionLogName, nullptr, nullptr, nullptr);
   fLeftSuppressorShellExtensionLog->SetVisAttributes(SuppressorVisAtt);

   G4SubtractionSolid* leftSuppressorExtension = SideSuppressorExtension("left", false);   // Left, Non-chopping // CALLED

   fLeftSuppressorExtensionLog = new G4LogicalVolume(leftSuppressorExtension, materialBGO,
                                                     leftSuppressorExtensionLogName, nullptr, nullptr, nullptr);
   fLeftSuppressorExtensionLog->SetVisAttributes(innardsVisAtt);

   fRightSuppressorExtensionAssembly->AddPlacedVolume(fRightSuppressorExtensionLog, fMoveNull, fRotateNull);
   fLeftSuppressorExtensionAssembly->AddPlacedVolume(fLeftSuppressorExtensionLog, fMoveNull, fRotateNull);

   //geometry objects for the following:
	std::array<G4RotationMatrix*, 8> rotateExtension;
   std::array<G4ThreeVector, 8>     moveExtension;

   x0 = -fCanFaceThickness / 2.0 - (shellSuppressorExtensionLength / 2.0 - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * tan(fBentEndAngle) / 2.0) * cos(fBentEndAngle) + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0) / tan(fBentEndAngle) - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * sin(fBentEndAngle) + fSuppShift + fSuppressorBackRadius - fSuppressorForwardRadius;

   y0 = -(shellSuppressorExtensionLength * tan(shellSuppressorExtensionAngle) / 2.0 + (fSuppressorForwardRadius + fHevimetTipThickness) * sin(fBentEndAngle)) / 2.0;

   z0 = -((fSuppressorExtensionThickness / 2.0 + fSuppressorShellThickness) / cos(fBentEndAngle) + (shellSuppressorExtensionLength / 2.0 - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * tan(fBentEndAngle) / 2.0) * sin(fBentEndAngle) - (fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0) / tan(fBentEndAngle) - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * sin(fBentEndAngle)) * tan(fBentEndAngle));

   if(fSuppressorPositionSelector == 0 && fIncludeExtensionSuppressors) {

      // these two parameters are for shifting the extensions back and out when in their BACK position
      G4double extensionBackShift = fAirBoxFrontLength - (fHevimetTipThickness + shellSuppressorExtensionLength + (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * tan(fBentEndAngle)) * cos(fBentEndAngle);

      G4double extensionRadialShift = extensionBackShift * tan(fBentEndAngle);

      // The suppressors are put into the back position

      x0 += extensionBackShift;

      z0 += extensionRadialShift;

      for(int i = 0; i < 4; i++) {
         rotateExtension[i * 2] = new G4RotationMatrix;
         rotateExtension[i * 2]->rotateZ(M_PI / 2.0);
         rotateExtension[i * 2]->rotateY(fBentEndAngle);
         rotateExtension[i * 2]->rotateX(M_PI - M_PI / 2.0 * i);

         G4double x = x0;
         G4double y = y0 * cos(i * M_PI / 2) + z0 * sin(i * M_PI / 2);
         G4double z = z0 * cos(i * M_PI / 2) - y0 * sin(i * M_PI / 2);

         moveExtension[i * 2] = G4ThreeVector(x, y, z);

         fExtensionSuppressorShellAssembly->AddPlacedVolume(fRightSuppressorShellExtensionLog, moveExtension[i * 2], rotateExtension[i * 2]);

         rotateExtension[i * 2 + 1] = new G4RotationMatrix;
         rotateExtension[i * 2 + 1]->rotateY(M_PI / 2.0);
         rotateExtension[i * 2 + 1]->rotateZ(M_PI / 2.0 + fBentEndAngle);
         rotateExtension[i * 2 + 1]->rotateX(M_PI - M_PI / 2.0 * i);

         x = x0;
         y = -z0 * cos(i * M_PI / 2) - y0 * sin(i * M_PI / 2);
         z = -y0 * cos(i * M_PI / 2) + z0 * sin(i * M_PI / 2);

         moveExtension[i * 2 + 1] = G4ThreeVector(x, y, z);

         fExtensionSuppressorShellAssembly->AddPlacedVolume(fLeftSuppressorShellExtensionLog, moveExtension[i * 2 + 1], rotateExtension[i * 2 + 1]);
      }
   } else if(fSuppressorPositionSelector == 1 && fIncludeExtensionSuppressors) {

      for(int i = 0; i < 4; i++) {
         rotateExtension[2 * i] = new G4RotationMatrix;
         rotateExtension[2 * i]->rotateZ(M_PI / 2.0);
         rotateExtension[2 * i]->rotateY(fBentEndAngle);
         rotateExtension[2 * i]->rotateX(M_PI - i * M_PI / 2.0);

         G4double x = x0;
         G4double y = y0 * cos(i * M_PI / 2) + z0 * sin(i * M_PI / 2);
         G4double z = z0 * cos(i * M_PI / 2) - y0 * sin(i * M_PI / 2);

         moveExtension[i * 2] = G4ThreeVector(x, y, z);

         fExtensionSuppressorShellAssembly->AddPlacedVolume(fRightSuppressorShellExtensionLog, moveExtension[2 * i], rotateExtension[2 * i]);

         rotateExtension[2 * i + 1] = new G4RotationMatrix;
         rotateExtension[2 * i + 1]->rotateY(M_PI / 2.0);
         rotateExtension[2 * i + 1]->rotateZ(M_PI / 2.0 + fBentEndAngle);
         rotateExtension[2 * i + 1]->rotateX(M_PI - i * M_PI / 2);

         x = x0;
         y = -z0 * cos(i * M_PI / 2) - y0 * sin(i * M_PI / 2);
         z = -y0 * cos(i * M_PI / 2) + z0 * sin(i * M_PI / 2);

         moveExtension[i * 2 + 1] = G4ThreeVector(x, y, z);

         fExtensionSuppressorShellAssembly->AddPlacedVolume(fLeftSuppressorShellExtensionLog, moveExtension[2 * i + 1], rotateExtension[2 * i + 1]);
      }
   }   //end if(detectors back) statement
}   //end ::ConstructNewSuppressorCasingWithShells

void DetectionSystemGriffin::ConstructNewSuppressorCasingDetectorSpecificDeadLayer(G4int det, G4int cry)
{
   G4String strdet = G4UIcommand::ConvertToString(det);
   G4String strcry = G4UIcommand::ConvertToString(cry);

   G4Material* materialBGO = G4Material::GetMaterial(fBGOMaterial);
   if(materialBGO == nullptr) {
      G4cout << " ----> BGO material " << fBGOMaterial << " not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   auto* SuppressorVisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
   SuppressorVisAtt->SetVisibility(true);
   auto* innardsVisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
   innardsVisAtt->SetVisibility(true);
   auto* backInnardsVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
   backInnardsVisAtt->SetVisibility(true);

   // first we add the four pieces of back suppressor, and their shells

   G4SubtractionSolid* backQuarterSuppressor = BackSuppressorQuarter();

   // Insert the structureMat shell first.  The shells must be given numbers, rather than
   // the suppressor pieces themselves.  As an error checking method, the suppressor
   // pieces are given a copy number value far out of range of any useful copy number.
   if(fIncludeBackSuppressors) {
      G4Material* backMaterial = G4Material::GetMaterial(fBackSuppressorMaterial);

      if(backMaterial == nullptr) {
         G4cout << " ----> Back suppressor material " << fBackSuppressorMaterial << " not found, cannot build the detector shell! " << G4endl;
         exit(1);
      }

      G4String backQuarterSuppressorName = "backDlsQuarterSuppressor_" + strdet + "_" + strcry + "_log";

      fBackQuarterSuppressorLog = new G4LogicalVolume(backQuarterSuppressor, backMaterial,
                                                      backQuarterSuppressorName, nullptr, nullptr, nullptr);

      fBackQuarterSuppressorLog->SetVisAttributes(backInnardsVisAtt);

      fSuppressorBackAssemblyCry[cry]->AddPlacedVolume(fBackQuarterSuppressorLog, fMoveNull, fRotateNull);
   }

   // now we add the side pieces of suppressor that taper off towards the front of the can

   G4String rightSuppressorCasingName = "rightDlsSuppressorCasing_" + strdet + "_" + strcry + "_log";
   G4String leftSuppressorCasingName  = "leftDlsSuppressorCasing_" + strdet + "_" + strcry + "_log";

   G4SubtractionSolid* rightSuppressor = FrontSlantSuppressor("right", false);   // Right, non-chopping. // CALLED

   fRightSuppressorLog = new G4LogicalVolume(rightSuppressor, materialBGO, rightSuppressorCasingName, nullptr, nullptr, nullptr);
   fRightSuppressorLog->SetVisAttributes(innardsVisAtt);

   G4SubtractionSolid* leftSuppressor = FrontSlantSuppressor("left", false);   // Left, non-chopping. // CALLED

   fLeftSuppressorLog = new G4LogicalVolume(leftSuppressor, materialBGO, leftSuppressorCasingName, nullptr, nullptr, nullptr);
   fLeftSuppressorLog->SetVisAttributes(innardsVisAtt);

   fRightSuppressorCasingAssemblyCry[cry]->AddPlacedVolume(fRightSuppressorLog, fMoveNull, fRotateNull);
   fLeftSuppressorCasingAssemblyCry[cry]->AddPlacedVolume(fLeftSuppressorLog, fMoveNull, fRotateNull);

   G4String rightSuppressorExtensionName = "rightDlsSuppressorExtension_" + strdet + "_" + strcry + "_log";
   G4String leftSuppressorExtensionName  = "leftDlsSuppressorExtension_" + strdet + "_" + strcry + "_log";

   // now we add the side pieces of suppressor that extend out in front of the can when it's in the back position

   G4SubtractionSolid* rightSuppressorExtension = SideSuppressorExtension("right", false);   // Right, non-chopping // CALLED

   fRightSuppressorExtensionLog = new G4LogicalVolume(rightSuppressorExtension, materialBGO,
                                                      rightSuppressorExtensionName, nullptr, nullptr, nullptr);
   fRightSuppressorExtensionLog->SetVisAttributes(innardsVisAtt);

   G4SubtractionSolid* leftSuppressorExtension = SideSuppressorExtension("left", false);   // Left, Non-chopping // CALLED

   fLeftSuppressorExtensionLog = new G4LogicalVolume(leftSuppressorExtension, materialBGO,
                                                     leftSuppressorExtensionName, nullptr, nullptr, nullptr);
   fLeftSuppressorExtensionLog->SetVisAttributes(innardsVisAtt);

   fRightSuppressorExtensionAssemblyCry[cry]->AddPlacedVolume(fRightSuppressorExtensionLog, fMoveNull, fRotateNull);
   fLeftSuppressorExtensionAssemblyCry[cry]->AddPlacedVolume(fLeftSuppressorExtensionLog, fMoveNull, fRotateNull);

}   //end ::ConstructNewSuppressorCasingDetectorSpecificDeadLayer

///////////////////////////////////////////////////////////////////////
// This is a messy place to put it, but these are the methods to make
// the electrodeMat layers between the crystals
///////////////////////////////////////////////////////////////////////
G4UnionSolid* DetectionSystemGriffin::InterCrystalelectrodeMatBack() const
{
   G4double distanceOfTheTriangleTips = fGermaniumSeparation / 2.0 + (fGermaniumWidth / 2.0 - fGermaniumShift)   //centre of the middle hole
                                        + sqrt(pow((fGermaniumOuterRadius + fTrianglePostsDistanceFromCrystals), 2.0) - pow(fGermaniumWidth / 2.0 - fGermaniumShift, 2.0));

   G4double extentOfTheElectrodeMatPieces = distanceOfTheTriangleTips - fTrianglePostsDistanceFromCrystals;

   auto* electrodeMatPiece1 = new G4Box("electrodeMatPiece1", extentOfTheElectrodeMatPieces,
                                         (fGermaniumLength - fGermaniumBentLength) / 2.0,
                                         fInterCrystalElectrodeMatThickness / 2.0);

   auto* electrodeMatPiece2 = new G4Box("electrodeMatPiece2", extentOfTheElectrodeMatPieces,
                                         (fGermaniumLength - fGermaniumBentLength) / 2.0,
                                         fInterCrystalElectrodeMatThickness / 2.0);

   auto* rotatePiece2 = new G4RotationMatrix;
   rotatePiece2->rotateY(M_PI / 2.0);

   G4ThreeVector moveZero(0, 0, 0);

   return new G4UnionSolid("interCrystalElectrodeMatBack", electrodeMatPiece1, electrodeMatPiece2, rotatePiece2, moveZero);
}   //end ::interCrystalelectrodeMatBack

G4UnionSolid* DetectionSystemGriffin::InterCrystalelectrodeMatFront() const
{
   G4double distanceOfTheTriangleTips = fGermaniumSeparation / 2.0 + (fGermaniumWidth / 2.0 - fGermaniumShift)   //centre of the middle hole
                                        + sqrt(pow((fGermaniumOuterRadius + fTrianglePostsDistanceFromCrystals), 2.0) - pow(fGermaniumWidth / 2.0 - fGermaniumShift, 2.0));

   auto* electrodeMatPiece1 = new G4Trd("electrodeMatPiece1", distanceOfTheTriangleTips - fTrianglePostsDistanceFromCrystals,
                                         fGermaniumWidth - fGermaniumBentLength * tan(fBentEndAngle),
                                         fGermaniumSeparation / 2.0,
                                         fGermaniumSeparation / 2.0, (fGermaniumBentLength - fElectrodeMatStartingDepth) / 2.0);

   auto* electrodeMatPiece2 = new G4Trd("electrodeMatPiece2", distanceOfTheTriangleTips - fTrianglePostsDistanceFromCrystals,
                                         fGermaniumWidth - fGermaniumBentLength * tan(fBentEndAngle),
                                         fGermaniumSeparation / 2.0,
                                         fGermaniumSeparation / 2.0, (fGermaniumBentLength - fElectrodeMatStartingDepth) / 2.0);

   auto* rotatePiece2 = new G4RotationMatrix;
   rotatePiece2->rotateZ(M_PI / 2.0);

   G4ThreeVector moveZero(0, 0, 0);

   return new G4UnionSolid("interCrystalElectrodeMatFront", electrodeMatPiece1, electrodeMatPiece2, rotatePiece2, moveZero);
}   //end ::interCrystalelectrodeMatFront

///////////////////////////////////////////////////////////////////////
// ConstructNewHeavyMet builds the heavy metal that goes on the front
// of the Suppressor. It only goes on when the extensions are in their
// forward position
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructNewHeavyMet()
{
   G4Material* materialHevimetal = G4Material::GetMaterial("Hevimetal");
   if(materialHevimetal == nullptr) {
      G4cout << " ----> Material Hevimetal not found, cannot build the detector shell! " << G4endl;
      exit(1);
   }

   auto* heviMetVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
   heviMetVisAtt->SetVisibility(true);

   G4SubtractionSolid* hevimet = NewHeavyMet();

   auto* rotateHevimet = new G4RotationMatrix;
   rotateHevimet->rotateY(-1.0 * M_PI / 2.0);

   G4ThreeVector moveHevimet(((fHevimetTipThickness / cos(fHevimetTipAngle)) * cos(fBentEndAngle - fHevimetTipAngle) + fSuppressorForwardRadius * tan(fHevimetTipAngle) * sin(fBentEndAngle)) / 2.0 - fAirBoxBackLength / 2.0 - fAirBoxFrontLength, 0.0, 0.0);

   // NOTE** The hevimet does not require the "shift" parameter, as the airBox has been designed
   // so that the hevimet sits right at the front, and has been moved accordingly. Also, the
   // "appliedBackShift" parameter is missing, as the hevimet only goes on in the "detector back" position

   fHevimetLog = new G4LogicalVolume(hevimet, materialHevimetal, "hevimetLog", nullptr, nullptr, nullptr);

   fHevimetLog->SetVisAttributes(heviMetVisAtt);

   fHevimetAssembly->AddPlacedVolume(fHevimetLog, moveHevimet, rotateHevimet);
}   //end ::ConstructHeavyMet

///////////////////////////////////////////////////////////////////////
// methods used in ConstructDetector()
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGriffin::SquareFrontFace() const
{

   G4double fixOverlap = 250 * nm;

   G4double halfThicknessX = fCanFaceThickness / 2.0;

   G4double halfLengthY = ((fDetectorTotalWidth / 2.0) - (fBentEndLength * tan(fBentEndAngle))) - fixOverlap;

   G4double halfLengthZ = halfLengthY;

   return new G4Box("frontFace", halfThicknessX, halfLengthY, halfLengthZ);
}   //end ::squareFrontFace

G4Trap* DetectionSystemGriffin::CornerWedge() const
{
   G4double heightY          = fCanFaceThickness;
   G4double lengthZ          = fDetectorTotalWidth - 2.0 * (fBentEndLength * tan(fBentEndAngle));
   G4double extensionLengthX = heightY * tan(fBentEndAngle);

   return new G4Trap("wedge", lengthZ, heightY, extensionLengthX, fWedgeDim);
}   //end ::cornerWedge

///////////////////////////////////////////////////////////////////////
// The bent side pieces attach to the front faceface they angle
// outwards, so the face is smaller than the main part of
// the can four are needed for a square-faced detector
///////////////////////////////////////////////////////////////////////
G4Para* DetectionSystemGriffin::BentSidePiece() const
{
   G4double halfLengthX = ((fBentEndLength - fCanFaceThickness) / cos(fBentEndAngle)) / 2.0;

   G4double halfLengthY = (fDetectorTotalWidth / 2.0) - (fBentEndLength * tan(fBentEndAngle)) - fCanFaceThickness * (1.0 - tan(fBentEndAngle));

   // Last two operations in "halfLengthY" are correction
   // factor to keep the bent pieces from overlapping
   G4double halfLengthZ    = fCanFaceThickness * cos(fBentEndAngle) / 2.0;
   G4double alphaAngle     = 0.0 * M_PI;
   G4double polarAngle     = fBentEndAngle;
   G4double azimuthalAngle = 0.0 * M_PI;

   return new G4Para("bentSidePiece", halfLengthX,
                                      halfLengthY, halfLengthZ, alphaAngle, polarAngle, azimuthalAngle);
}   //end ::bentSidePiece

// this is a conic piece that attaches two bent end pieces
G4Cons* DetectionSystemGriffin::RoundedEndEdge() const
{
   G4double holeEliminator = fCanFaceThickness * (1.0 - tan(fBentEndAngle));

   // this is correction factor to avoid holes caused by bent pieces
   G4double halfLengthZ         = (fBentEndLength - fCanFaceThickness) / 2.0;
   G4double insideRadiusAtApex  = 0.0;
   G4double outsideRadiusAtApex = fCanFaceThickness * tan(fBentEndAngle) + holeEliminator;
   G4double outsideRadiusAtBase = fBentEndLength * tan(fBentEndAngle) + holeEliminator;
   G4double insideRadiusAtBase  = outsideRadiusAtBase - fCanFaceThickness;
   G4double startAngle          = 0.0 * M_PI;
   G4double finalAngle          = 0.5 * M_PI;

   return new G4Cons("roundedEndEdge", insideRadiusAtApex,
                                       outsideRadiusAtApex, insideRadiusAtBase, outsideRadiusAtBase,
                                       halfLengthZ, startAngle, finalAngle);
}   //end ::roundedEndEdge

G4Tubs* DetectionSystemGriffin::CornerTube() const
{
   G4double outerRadius = fBentEndLength * tan(fBentEndAngle);
   G4double innerRadius = outerRadius - fCanSideThickness;
   G4double halfHeightZ = (fDetectorTotalLength - fBentEndLength - fRearPlateThickness) / 2.0;

   G4double startAngle = 0.0 * M_PI;
   G4double finalAngle = 0.5 * M_PI;

   return new G4Tubs("cornerTube", innerRadius, outerRadius, halfHeightZ, startAngle, finalAngle);
}   //end ::cornerTube

G4Box* DetectionSystemGriffin::SidePanel() const
{
   G4double halfLengthX    = (fDetectorTotalLength - fBentEndLength - fRearPlateThickness) / 2.0;
   G4double halfThicknessY = (fCanSideThickness) / 2.0;
   G4double halfWidthZ     = (fDetectorTotalWidth) / 2.0 - (fBentEndLength * tan(fBentEndAngle));

   return new G4Box("sidePanel", halfLengthX, halfThicknessY, halfWidthZ);
}   //end ::sidePanel

G4SubtractionSolid* DetectionSystemGriffin::RearPlate() const
{
   G4double halfWidthX     = fDetectorTotalWidth / 2.0;
   G4double halfWidthY     = halfWidthX;
   G4double halfThicknessZ = fRearPlateThickness / 2.0;

   auto* plate = new G4Box("plate", halfWidthX, halfWidthY, halfThicknessZ);

   G4double innerRadius = 0.0;
   G4double outerRadius = fColdFingerOuterShellRadius;
   G4double halfHeightZ = halfThicknessZ + 1.0 * cm;   // +1.0*cm just to make sure the hole goes completely through the plate
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   auto* hole = new G4Tubs("hole", innerRadius, outerRadius, halfHeightZ, startAngle, finalAngle);

   return new G4SubtractionSolid("rearPlate", plate, hole);
}   //end ::rearPlate

G4Tubs* DetectionSystemGriffin::ColdFingerShell() const
{
   G4double outerRadius = fColdFingerOuterShellRadius;
   G4double innerRadius = outerRadius - fColdFingerShellThickness;
   G4double halfLengthZ = fColdFingerShellLength / 2.0;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("coldFingerShell", innerRadius,
                                        outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::coldFingerShell

G4Tubs* DetectionSystemGriffin::LiquidNitrogenTank() const
{
   G4double innerRadius = fCoolantRadius - fCoolantThickness;
   G4double outerRadius = fCoolantRadius;
   G4double halfLengthZ = (fCoolantLength - 2.0 * fCoolantThickness) / 2.0;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("tank", innerRadius, outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::liquidNitrogenTank

G4Tubs* DetectionSystemGriffin::LiquidNitrogenTankLid() const
{
   G4double innerRadius = 0.0 * mm;
   G4double outerRadius = fCoolantRadius;
   G4double halfLengthZ = (fCoolantThickness) / 2.0;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("tank", innerRadius, outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::liquidNitrogenTankLid

G4Tubs* DetectionSystemGriffin::LiquidNitrogen() const
{
   G4double innerRadius = 0.0 * mm;
   G4double outerRadius = fCoolantRadius - fCoolantThickness;
   G4double halfLengthZ = (fCoolantLength - 2.0 * fCoolantThickness) / 2.0;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   return new G4Tubs("tank", innerRadius, outerRadius, halfLengthZ, startAngle, finalAngle);
}   //end ::liquidNitrogen

///////////////////////////////////////////////////////////////////////
// methods used in ConstructBasicDetectorBlock()
// Both the rectangularSegment and the trapezoidalSegment join
// to form a UnionSolid
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGriffin::RectangularSegment() const
{
   G4double halfLengthX = fDetectorBlockHeight / 2.0;
   G4double halfLengthY = fDetectorBlockLength / 2.0;
   G4double halfLengthZ = halfLengthY;   // Since it is symmetric

   return new G4Box("rectangularSegment", halfLengthX, halfLengthY, halfLengthZ);
}   // end ::rectangularSegment

G4Trd* DetectionSystemGriffin::TrapezoidalSegment() const
{
   G4double halfBaseX   = fDetectorBlockLength / 2.0;
   G4double halfTopX    = halfBaseX - (fDetectorBlockTrapezoidalHeight * tan(fBentEndAngle));
   G4double halfBaseY   = halfBaseX;   // Since it is symmetric
   G4double halfTopY    = halfTopX;
   G4double halfHeightZ = fDetectorBlockTrapezoidalHeight / 2.0;

   return new G4Trd("trapezoidalSegment",
                                         halfBaseX, halfTopX, halfBaseY, halfTopY, halfHeightZ);
}   //end ::trapezoidalSegment

///////////////////////////////////////////////////////////////////////
// methods used in ConstructComplexDetectorBlock()
// Starting with a rectangle, it gets chopped by an off-centered
// cylinder, and then the two edges are chopped diagonaly
///////////////////////////////////////////////////////////////////////
G4SubtractionSolid* DetectionSystemGriffin::QuarterDetector() const
{

   G4double halfWidthX  = fGermaniumWidth / 2.0;
   G4double halfWidthY  = halfWidthX;
   G4double halfLengthZ = fGermaniumLength / 2.0;

   auto* rectangularGermanium = new G4Box("rectangularGermanium", halfWidthX, halfWidthY, halfLengthZ);

   G4double outerRadius = ((fGermaniumWidth + fGermaniumShift) * 1.5) / 2.0;

   // 1.5 is almost root 2, so this makes sure the corners are chopped off

   G4double innerRadius = fGermaniumOuterRadius;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   G4ThreeVector moveChoppingCylinder(fGermaniumShift, -(fGermaniumShift), 0);
   auto*       choppingCylinder = new G4Tubs("choppingCylinder", innerRadius,
                                               outerRadius, halfLengthZ + 1.0 * cm, startAngle, finalAngle);

   auto* germaniumRoundedCorners = new G4SubtractionSolid("germaniumRoundedCorners",
                                                                        rectangularGermanium, choppingCylinder, nullptr, moveChoppingCylinder);

   G4double baseInnerRadius = fGermaniumOuterRadius;
   G4double baseOuterRadius = 2.0 * baseInnerRadius;

   // G4double tipInnerRadius = baseInnerRadius - fGermaniumBentLength*tan(fBentEndAngle);
   G4double tipInnerRadius = baseInnerRadius - fGermaniumCornerConeEndLength * tan(fBentEndAngle);
   G4double tipOuterRadius = baseOuterRadius;

   // G4double conHalfLengthZ = fGermaniumBentLength/2.0;
   G4double conHalfLengthZ = fGermaniumCornerConeEndLength / 2.0 + fQuarterDetectorCxn;

   G4double initialAngle = acos((fGermaniumWidth / 2.0 + fGermaniumShift) / fGermaniumOuterRadius);
   G4double totalAngle   = M_PI / 2.0 - 2.0 * initialAngle;

   auto* roundedEdge = new G4Cons("roundedEdge", baseInnerRadius,
                                    baseOuterRadius, tipInnerRadius, tipOuterRadius,
                                    conHalfLengthZ, initialAngle, totalAngle);

   auto* rotateRoundedEdge = new G4RotationMatrix;
   rotateRoundedEdge->rotateZ(-M_PI / 2.0);

   G4ThreeVector moveRoundedEdge(fGermaniumShift, -fGermaniumShift,
                                 fGermaniumLength / 2.0 - fGermaniumCornerConeEndLength / 2.0 + fQuarterDetectorCxnB);

   auto* germaniumRoundedEdge = new G4SubtractionSolid("germaniumRoundedEdge",
                                                                     germaniumRoundedCorners, roundedEdge, rotateRoundedEdge,
                                                                     moveRoundedEdge);

   // now we make the diagonal slices
   G4double halfChopPieceWidthX  = fGermaniumWidth / 2.0;
   G4double halfChopPieceWidthY  = fGermaniumWidth / 2.0;
   G4double halfChopPieceLengthZ = (fGermaniumBentLength / cos(fBentEndAngle)) / 2.0;
   auto*   chopPiece            = new G4Box("chopPiece", halfChopPieceWidthX,
                                             halfChopPieceWidthY, halfChopPieceLengthZ);

   auto* rotateChopPiece1 = new G4RotationMatrix;
   rotateChopPiece1->rotateX(-(fBentEndAngle));

   G4ThreeVector moveChopPiece1(0, fGermaniumWidth / 2.0 + (fGermaniumBentLength / tan(fBentEndAngle) + fGermaniumWidth * cos(fBentEndAngle)) / 2.0 - ((fGermaniumBentLength / cos(fBentEndAngle)) / 2.0) / sin(fBentEndAngle),
                                fGermaniumLength / 2.0 - fGermaniumBentLength + (fGermaniumBentLength + fGermaniumWidth * sin(fBentEndAngle)) / 2.0);

   auto* choppedGermanium1 = new G4SubtractionSolid("choppedGermanium1",
                                                                  germaniumRoundedEdge, chopPiece, rotateChopPiece1, moveChopPiece1);

   auto* rotateChopPiece2 = new G4RotationMatrix;
   rotateChopPiece2->rotateY(-(fBentEndAngle));

   G4ThreeVector moveChopPiece2(-(fGermaniumWidth / 2.0 + (fGermaniumBentLength / tan(fBentEndAngle) + fGermaniumWidth * cos(fBentEndAngle)) / 2.0 - ((fGermaniumBentLength / cos(fBentEndAngle)) / 2.0) / sin(fBentEndAngle)), 0,
                                fGermaniumLength / 2.0 - fGermaniumBentLength + (fGermaniumBentLength + fGermaniumWidth * sin(fBentEndAngle)) / 2.0);

   auto* choppedGermanium2 = new G4SubtractionSolid("choppedGermanium2",
                                                                  choppedGermanium1, chopPiece, rotateChopPiece2, moveChopPiece2);

   // now we make the hole that will go in the back of each quarter detector, /////////////////////////////////////////////////////////
   startAngle               = 0.0 * M_PI;
   finalAngle               = 2.0 * M_PI;
   G4double holeRadius      = fGermaniumHoleRadius;
   G4double holeHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace) / 2.0;

   auto*       holeTubs = new G4Tubs("holeTubs", 0.0, holeRadius, holeHalfLengthZ, startAngle, finalAngle);
   G4ThreeVector moveHole(fGermaniumShift, -(fGermaniumShift), -((fGermaniumHoleDistFromFace)));

   auto* choppedGermanium3 = new G4SubtractionSolid("choppedGermanium3",
                                                                  choppedGermanium2, holeTubs, nullptr, moveHole);

   startAngle                    = 0.0 * M_PI;
   finalAngle                    = 2.0 * M_PI;
   G4double deadLayerRadius      = fGermaniumHoleRadius + fInnerDeadLayerThickness;
   G4double deadLayerHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace + fInnerDeadLayerThickness) / 2.0;

   // now dead layer
   auto* deadLayerTubs = new G4Tubs("deadLayerTubs", 0.0,
                                      deadLayerRadius, deadLayerHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveDeadLayer(fGermaniumShift, -(fGermaniumShift),
                               -((fGermaniumHoleDistFromFace - fInnerDeadLayerThickness) / 2.0));

   return new G4SubtractionSolid("choppedGermanium4",
                                                                  choppedGermanium3, deadLayerTubs, nullptr, moveDeadLayer);
}   //end ::quarterDetector

G4SubtractionSolid* DetectionSystemGriffin::QuarterSpecificDeadLayerDetector(G4int det, G4int cry)
{

   G4double halfWidthX  = fGermaniumWidth / 2.0;
   G4double halfWidthY  = halfWidthX;
   G4double halfLengthZ = fGermaniumLength / 2.0;

   auto* rectangularGermanium = new G4Box("rectangularGermanium", halfWidthX, halfWidthY, halfLengthZ);

   G4double outerRadius = ((fGermaniumWidth + fGermaniumShift) * 1.5) / 2.0;

   // 1.5 is almost root 2, so this makes sure the corners are chopped off

   G4double innerRadius = fGermaniumOuterRadius;
   G4double startAngle  = 0.0 * M_PI;
   G4double finalAngle  = 2.0 * M_PI;

   G4ThreeVector moveChoppingCylinder(fGermaniumShift, -(fGermaniumShift), 0);
   auto*       choppingCylinder = new G4Tubs("choppingCylinder", innerRadius,
                                               outerRadius, halfLengthZ + 1.0 * cm, startAngle, finalAngle);

   auto* germaniumRoundedCorners = new G4SubtractionSolid("germaniumRoundedCorners",
                                                                        rectangularGermanium, choppingCylinder, nullptr, moveChoppingCylinder);

   G4double baseInnerRadius = fGermaniumOuterRadius;
   G4double baseOuterRadius = 2.0 * baseInnerRadius;

   G4double tipInnerRadius = baseInnerRadius - fGermaniumCornerConeEndLength * tan(fBentEndAngle);
   G4double tipOuterRadius = baseOuterRadius;

   G4double conHalfLengthZ = fGermaniumCornerConeEndLength / 2.0 + fQuarterDetectorCxn;

   G4double initialAngle = acos((fGermaniumWidth / 2.0 + fGermaniumShift) / fGermaniumOuterRadius);
   G4double totalAngle   = M_PI / 2.0 - 2.0 * initialAngle;

   auto* roundedEdge = new G4Cons("roundedEdge", baseInnerRadius,
                                    baseOuterRadius, tipInnerRadius, tipOuterRadius,
                                    conHalfLengthZ, initialAngle, totalAngle);

   auto* rotateRoundedEdge = new G4RotationMatrix;
   rotateRoundedEdge->rotateZ(-M_PI / 2.0);

   G4ThreeVector moveRoundedEdge(fGermaniumShift, -fGermaniumShift,
                                 fGermaniumLength / 2.0 - fGermaniumCornerConeEndLength / 2.0 + fQuarterDetectorCxnB);

   auto* germaniumRoundedEdge = new G4SubtractionSolid("germaniumRoundedEdge",
                                                                     germaniumRoundedCorners, roundedEdge, rotateRoundedEdge,
                                                                     moveRoundedEdge);

   // now we make the diagonal slices
   G4double halfChopPieceWidthX  = fGermaniumWidth / 2.0;
   G4double halfChopPieceWidthY  = fGermaniumWidth / 2.0;
   G4double halfChopPieceLengthZ = (fGermaniumBentLength / cos(fBentEndAngle)) / 2.0;
   auto*   chopPiece            = new G4Box("chopPiece", halfChopPieceWidthX,
                                             halfChopPieceWidthY, halfChopPieceLengthZ);

   auto* rotateChopPiece1 = new G4RotationMatrix;
   rotateChopPiece1->rotateX(-(fBentEndAngle));

   G4ThreeVector moveChopPiece1(0, fGermaniumWidth / 2.0 + (fGermaniumBentLength / tan(fBentEndAngle) + fGermaniumWidth * cos(fBentEndAngle)) / 2.0 - ((fGermaniumBentLength / cos(fBentEndAngle)) / 2.0) / sin(fBentEndAngle),
                                fGermaniumLength / 2.0 - fGermaniumBentLength + (fGermaniumBentLength + fGermaniumWidth * sin(fBentEndAngle)) / 2.0);

   auto* choppedGermanium1 = new G4SubtractionSolid("choppedGermanium1",
                                                                  germaniumRoundedEdge, chopPiece, rotateChopPiece1, moveChopPiece1);

   auto* rotateChopPiece2 = new G4RotationMatrix;
   rotateChopPiece2->rotateY(-(fBentEndAngle));

   G4ThreeVector moveChopPiece2(-(fGermaniumWidth / 2.0 + (fGermaniumBentLength / tan(fBentEndAngle) + fGermaniumWidth * cos(fBentEndAngle)) / 2.0 - ((fGermaniumBentLength / cos(fBentEndAngle)) / 2.0) / sin(fBentEndAngle)), 0,
                                fGermaniumLength / 2.0 - fGermaniumBentLength + (fGermaniumBentLength + fGermaniumWidth * sin(fBentEndAngle)) / 2.0);

   auto* choppedGermanium2 = new G4SubtractionSolid("choppedGermanium2",
                                                                  choppedGermanium1, chopPiece, rotateChopPiece2, moveChopPiece2);

   // now we make the hole that will go in the back of each quarter detector, /////////////////////////////////////////////////////////
   startAngle               = 0.0 * M_PI;
   finalAngle               = 2.0 * M_PI;
   G4double holeRadius      = fGermaniumHoleRadius;
   G4double holeHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace) / 2.0;

   auto*       holeTubs = new G4Tubs("holeTubs", 0.0, holeRadius, holeHalfLengthZ, startAngle, finalAngle);
   G4ThreeVector moveHole(fGermaniumShift, -(fGermaniumShift), -((fGermaniumHoleDistFromFace)));

   auto* choppedGermanium3 = new G4SubtractionSolid("choppedGermanium3",
                                                                  choppedGermanium2, holeTubs, nullptr, moveHole);

   startAngle                    = 0.0 * M_PI;
   finalAngle                    = 2.0 * M_PI;
   G4double deadLayerRadius      = fGermaniumHoleRadius + fGriffinDeadLayers[det][cry];
   G4double deadLayerHalfLengthZ = (fGermaniumLength - fGermaniumHoleDistFromFace + fGriffinDeadLayers[det][cry]) / 2.0;

   // now dead layer
   auto* deadLayerTubs = new G4Tubs("deadLayerTubs", 0.0,
                                      deadLayerRadius, deadLayerHalfLengthZ, startAngle, finalAngle);

   G4ThreeVector moveDeadLayer(fGermaniumShift, -(fGermaniumShift),
                               -((fGermaniumHoleDistFromFace - fGriffinDeadLayers[det][cry]) / 2.0));

   auto* choppedGermanium4 = new G4SubtractionSolid("choppedGermanium4",
                                                                  choppedGermanium3, deadLayerTubs, nullptr, moveDeadLayer);

   return choppedGermanium4;

}   //end ::quarterSpecificDeadLayerDetector

///////////////////////////////////////////////////////////////////////
// methods used in ConstructNewSuppressorCasingWithShells
///////////////////////////////////////////////////////////////////////

G4SubtractionSolid* DetectionSystemGriffin::ShellForBackSuppressorQuarter() const
{
   G4double halfThicknessX         = fBackBGOThickness / 2.0 + fSuppressorShellThickness;
   G4double halfLengthY            = fDetectorTotalWidth / 4.0;
   G4double halfLengthZ            = halfLengthY;
   auto*   quarterSuppressorShell = new G4Box("quarterSuppressorShell", halfThicknessX, halfLengthY, halfLengthZ);

   G4double shellHoleRadius = fColdFingerOuterShellRadius;
   auto*  shellHole       = new G4Tubs("shellHole", 0, shellHoleRadius, halfThicknessX + 1.0 * cm, 0.0 * M_PI, 2.0 * M_PI);

   auto* rotateShellHole = new G4RotationMatrix;
   rotateShellHole->rotateY(M_PI / 2.0);
   G4ThreeVector moveShellHole(0, -halfLengthY, halfLengthZ);

   auto* quarterSuppressorShellWithHole = new G4SubtractionSolid("quarterSuppressorShellWithHole",
                                                                               quarterSuppressorShell, shellHole, rotateShellHole, moveShellHole);

   //now we need to cut out inner cavity, first we define the cavity
   halfThicknessX           = fBackBGOThickness / 2.0;
   halfLengthY              = fDetectorTotalWidth / 4.0 - fSuppressorShellThickness;
   halfLengthZ              = halfLengthY;
   auto* quarterSuppressor = new G4Box("quarterSuppressor", halfThicknessX, halfLengthY, halfLengthZ);

   G4ThreeVector moveCut(0, 0, 0);

   // cut
   auto* quarterSuppressorShellWithHoleAndCavity = new G4SubtractionSolid("quarterSuppressorShellWithHoleAndCavity",
                                                                                        quarterSuppressorShellWithHole, quarterSuppressor, nullptr, moveCut);

   return quarterSuppressorShellWithHoleAndCavity;

}   //end ::shellForBackSuppressorQuarter

G4SubtractionSolid* DetectionSystemGriffin::ShellForFrontSlantSuppressor(const G4String& sidePosition)
{
   // Change some values to accomodate the shells
   // Replacement for sideSuppressorLength
   G4double shellSideSuppressorShellLength = fSideSuppressorLength + (fSuppressorShellThickness * 2.0);

   G4double lengthZ        = shellSideSuppressorShellLength;
   G4double lengthY        = fSideBGOThickness + fSuppressorShellThickness * 2.0;
   G4double lengthLongerX  = fDetectorTotalWidth / 2.0 + fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0;
   G4double lengthShorterX = lengthLongerX - fSideBGOThickness - fSuppressorShellThickness * 2.0;

   auto* suppressorShell = new G4Trap("suppressorShell", lengthZ, lengthY, lengthLongerX, lengthShorterX);

   G4double halfLengthX    = lengthLongerX / 2.0 + 1.0 * cm;
   G4double halfThicknessY = fSideBGOThickness / 2.0 + fSuppressorShellThickness;
   G4double halfLengthZ    = ((fSideBGOThickness + fSuppressorShellThickness * 2.0 - fBGOChoppedTip) / sin(fBentEndAngle)) / 2.0;

   auto* choppingShellBox = new G4Box("choppingShellBox", halfLengthX, halfThicknessY, halfLengthZ);

   auto* rotateChoppingShellBox = new G4RotationMatrix;

   G4double y0 = fSideBGOThickness / 2.0 - fSuppressorShellThickness * sin(fBentEndAngle) - fBGOChoppedTip - 0.5 * sqrt(pow((2.0 * halfLengthZ), 2.0) + pow((2.0 * halfThicknessY), 2.0)) * cos(M_PI / 2.0 - fBentEndAngle - atan(halfThicknessY / halfLengthZ));

   G4double      z0 = lengthZ / 2.0 - 0.5 * sqrt(pow((2.0 * halfLengthZ), 2.0) + pow((2.0 * halfThicknessY), 2.0)) * sin(M_PI / 2.0 - fBentEndAngle - atan(halfThicknessY / halfLengthZ));
   G4ThreeVector moveChoppingShellBox;

   if(sidePosition == "left") {
      rotateChoppingShellBox->rotateX(fBentEndAngle);
      moveChoppingShellBox = G4ThreeVector(0, y0, z0);
   } else if(sidePosition == "right") {
      rotateChoppingShellBox->rotateX(-fBentEndAngle);
      moveChoppingShellBox = G4ThreeVector(0, y0, -z0);
   }

   auto* sideSuppressorShell = new G4SubtractionSolid("sideSuppressorShell",
                                                                    suppressorShell, choppingShellBox, rotateChoppingShellBox, moveChoppingShellBox);

   G4SubtractionSolid* sideSuppressorShellWithCavity = nullptr;

   // cut out cavity
   if(sidePosition == "left") {
      G4ThreeVector moveCut(-(fSuppressorShellThickness + (fExtraCutLength / 2.0 - fSuppressorShellThickness) / 2.0), 0, (fSuppressorShellThickness / 2.0));

      sideSuppressorShellWithCavity = new G4SubtractionSolid("sideSuppressorShellWithCavity", sideSuppressorShell,
                                                             FrontSlantSuppressor("left", true), nullptr, moveCut);   // chopping, left from frontSlantSuppressor.
   } else if(sidePosition == "right") {
      G4ThreeVector moveCut(-(fSuppressorShellThickness + (fExtraCutLength / 2.0 - fSuppressorShellThickness) / 2.0), 0, -(fSuppressorShellThickness / 2.0));

      sideSuppressorShellWithCavity = new G4SubtractionSolid("sideSuppressorShellWithCavity", sideSuppressorShell,
                                                             FrontSlantSuppressor("right", true), nullptr, moveCut);   // chopping, right from frontSlantSuppressor.
   }

   return sideSuppressorShellWithCavity;
}   // end ::shellForFrontSlantSuppressor

G4SubtractionSolid* DetectionSystemGriffin::ShellForSuppressorExtension(const G4String& sidePosition)
{

   // Replacement for suppressorExtensionLength
   G4double shellSuppressorShellExtensionLength = fSuppressorExtensionLength + (fSuppressorShellThickness * 2.0) * (1.0 / tan(fBentEndAngle) - tan(fBentEndAngle));

   G4double thicknessZ = fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0;
   G4double lengthY    = shellSuppressorShellExtensionLength;

   G4double longerLengthX = (fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness + fSuppressorShellThickness * 2.0) / tan(fBentEndAngle) - (fSuppressorExtensionThickness + fSuppressorShellThickness * 2.0) * sin(fBentEndAngle)) * tan(fBentEndAngle);

   G4double shorterLengthX = (fSuppressorForwardRadius + fHevimetTipThickness) * sin(fBentEndAngle);

   auto* uncutExtensionShell = new G4Trap("uncutExtensionShell", thicknessZ, lengthY, longerLengthX, shorterLengthX);

   // because these pieces are rotated in two planes, there are multiple angles that need to be calculated to make sure
   // all of the extensions join up

   G4double beta = atan((longerLengthX - shorterLengthX) / (lengthY));
   G4double phi  = atan(1 / cos(fBentEndAngle));

   G4double choppingHalfLengthX = (thicknessZ / sin(phi)) / 2.0;
   G4double choppingHalfLengthY = lengthY / (2.0 * cos(beta));
   G4double choppingHalfLengthZ = choppingHalfLengthX;
   G4double yAngle              = -beta;
   G4double xAngle              = 0.0;
   G4double zAngle              = 0.;

   if(sidePosition == "left") {
      zAngle = M_PI / 2.0 - phi;
	} else if(sidePosition == "right") {
      zAngle = phi - M_PI / 2.0;
	}

   auto* choppingParaShell = new G4Para("choppingParaShell", choppingHalfLengthX,
                                          choppingHalfLengthY, choppingHalfLengthZ, yAngle, zAngle, xAngle);

   G4ThreeVector moveParaShell(((longerLengthX - shorterLengthX) / 2.0 + shorterLengthX) / 2.0 + choppingHalfLengthX - choppingHalfLengthX * cos(phi), 0.0, 0.0);

   auto* extensionShell = new G4SubtractionSolid("extensionShell",
                                                               uncutExtensionShell, choppingParaShell, nullptr, moveParaShell);

   G4ThreeVector       moveCut;
   G4SubtractionSolid* extensionSuppressorShellWithCavity    = nullptr;
   G4SubtractionSolid* extensionSuppressorShellWithCavityPre = nullptr;

   // cut out cavity ----
   // We make two cuts, one "crude" one from a G4Trap, and a more refined cut using a G4SubtractionSolid. There are two reasons why we make two cuts.
   // Firstly, the shape of these suppressors are strange. To make the inside cut large enough to have an opening such that the left and right suppressors touch would
   // require some clever math, and that's a pain. The second (less lazy) reason is Geant4 doesn't like to make a G4SubtractionSolid out of two G4SubtractionSolids.
   // My experience is that the visulization will fail if two sides of the G4SubtractionSolids are close, that is do not have a LARGE (>~10mm) overlap. To get around this,
   // I simply make a "crude" first cut using a G4Trap.
   // The cuts are slightly larger than suppressor to allow for small gaps between the shell and the suppressor.
   if(sidePosition == "left") {
      moveCut                               = G4ThreeVector(-30.0 * mm, 0.40 * mm, 0.25 * mm);
      extensionSuppressorShellWithCavityPre = new G4SubtractionSolid("extensionSuppressorShellWithCavityPre", extensionShell, SideSuppressorExtensionUncut(), nullptr, (moveCut) / 2.0);   // Left, Chopping sideSuppressorExtension

      moveCut                            = G4ThreeVector(-1.15 * mm, 0.40 * mm, 0.25 * mm);
      extensionSuppressorShellWithCavity = new G4SubtractionSolid("extensionSuppressorShellWithCavity", extensionSuppressorShellWithCavityPre, SideSuppressorExtension("left", true), nullptr, (moveCut) / 2.0);
   } else if(sidePosition == "right") {
      moveCut                               = G4ThreeVector(-30.0 * mm, 0.40 * mm, -0.25 * mm);
      extensionSuppressorShellWithCavityPre = new G4SubtractionSolid("extensionSuppressorShellWithCavityPre", extensionShell, SideSuppressorExtensionUncut(), nullptr, (moveCut) / 2.0);   // Left, Chopping sideSuppressorExtension

      moveCut                            = G4ThreeVector(-1.15 * mm, 0.40 * mm, -0.25 * mm);
      extensionSuppressorShellWithCavity = new G4SubtractionSolid("extensionSuppressorShellWithCavity", extensionSuppressorShellWithCavityPre, SideSuppressorExtension("right", true), nullptr, (moveCut) / 2.0);
   }

   return extensionSuppressorShellWithCavity;
}   // end ::shellForSuppressorExtension

// Back to making the suppressor volumes themselves
G4SubtractionSolid* DetectionSystemGriffin::BackSuppressorQuarter() const
{
   G4double halfThicknessX    = fBackBGOThickness / 2.0;
   G4double halfLengthY       = fDetectorTotalWidth / 4.0 - fSuppressorShellThickness;
   G4double halfLengthZ       = halfLengthY;
   auto*   quarterSuppressor = new G4Box("quarterSuppressor", halfThicknessX, halfLengthY, halfLengthZ);

   G4double holeRadius = fColdFingerOuterShellRadius;
   auto*  hole       = new G4Tubs("hole", 0, holeRadius, halfThicknessX + 1.0 * cm, 0.0 * M_PI, 2.0 * M_PI);

   auto* rotateHole = new G4RotationMatrix;
   rotateHole->rotateY(M_PI / 2.0);
   G4ThreeVector moveHole(0, -halfLengthY, halfLengthZ);

   auto* quarterSuppressorWithHole = new G4SubtractionSolid("quarterSuppressorWithHole",
                                                                          quarterSuppressor, hole, rotateHole, moveHole);

   return quarterSuppressorWithHole;

}   //end ::backSuppressorQuarter

G4SubtractionSolid* DetectionSystemGriffin::FrontSlantSuppressor(const G4String& sidePosition, G4bool choppingSuppressor) const
{
   // If chop is true, it will be a chopping suppressor.

   G4double lengthZ       = fSideSuppressorLength;
   G4double lengthY       = fSideBGOThickness;
   G4double lengthLongerX = 0;

   if(choppingSuppressor) {
      lengthLongerX = fDetectorTotalWidth / 2.0 + fBGOCanSeperation + fSideBGOThickness + fExtraCutLength / 2.0;
	} else {
      lengthLongerX = fDetectorTotalWidth / 2.0 + fBGOCanSeperation + fSideBGOThickness;
	}

   G4double lengthShorterX = lengthLongerX - fSideBGOThickness;

   auto* suppressor = new G4Trap("suppressor", lengthZ, lengthY, lengthLongerX, lengthShorterX);

   G4double halfLengthX    = lengthLongerX / 2.0 + 1.0 * cm;
   G4double halfThicknessY = fSideBGOThickness / 2.0;
   G4double halfLengthZ    = ((fSideBGOThickness - fBGOChoppedTip) / sin(fBentEndAngle)) / 2.0;

   auto* choppingBox = new G4Box("choppingBox", halfLengthX, halfThicknessY, halfLengthZ);

   auto* rotateChoppingBox = new G4RotationMatrix;
   G4ThreeVector     moveChoppingBox;

   G4double y0 = fSideBGOThickness / 2.0 - fBGOChoppedTip - 0.5 * sqrt(pow((2.0 * halfLengthZ), 2.0) + pow((2.0 * halfThicknessY), 2.0)) * cos(M_PI / 2.0 - fBentEndAngle - atan(halfThicknessY / halfLengthZ));

   G4double z0 = lengthZ / 2.0 - 0.5 * sqrt(pow((2.0 * halfLengthZ), 2.0) + pow((2.0 * halfThicknessY), 2.0)) * sin(M_PI / 2.0 - fBentEndAngle - atan(halfThicknessY / halfLengthZ));

   if(sidePosition == "left") {
      rotateChoppingBox->rotateX(fBentEndAngle);
      moveChoppingBox = G4ThreeVector(0, y0, z0);
   } else if(sidePosition == "right") {
      rotateChoppingBox->rotateX(-fBentEndAngle);
      moveChoppingBox = G4ThreeVector(0, y0, -z0);
   }

   auto* sideSuppressor = new G4SubtractionSolid("sideSuppressor",
                                                               suppressor, choppingBox, rotateChoppingBox, moveChoppingBox);

   return sideSuppressor;

}   // end ::frontSlantSuppressor

G4SubtractionSolid* DetectionSystemGriffin::SideSuppressorExtension(const G4String& sidePosition, G4bool choppingSuppressor) const
{
   G4double thicknessZ = fSuppressorExtensionThickness;
   G4double lengthY    = fSuppressorExtensionLength;

   G4double longerLengthX = ((fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness) / tan(fBentEndAngle) - fSuppressorExtensionThickness * sin(fBentEndAngle)) * tan(fBentEndAngle) - fBGOCanSeperation * 2.0);   // - fBGOCanSeperation

   G4double shorterLengthX = ((fSuppressorForwardRadius + fHevimetTipThickness) * sin(fBentEndAngle) - fBGOCanSeperation * 2.0);   // - fBGOCanSeperation

   if(choppingSuppressor) {
      thicknessZ += 0.1 * mm;
      lengthY += 0.1 * mm;
   }

   auto* uncutExtension = new G4Trap("uncutExtension", thicknessZ, lengthY, longerLengthX, shorterLengthX);

   // because these pieces are rotated in two planes, there are multiple angles that need to be calculated to make sure
   // all of the extensions join up
   G4double beta = atan((longerLengthX - shorterLengthX) / (lengthY));
   G4double phi  = atan(1 / cos(fBentEndAngle));

   G4double choppingHalfLengthX = thicknessZ / (2.0 * sin(phi));
   G4double choppingHalfLengthY = lengthY / (2.0 * cos(beta));
   G4double choppingHalfLengthZ = choppingHalfLengthX;

   G4double xAngle = 0.0;
   G4double yAngle = -beta;
   G4double zAngle = phi - M_PI / 2.0;

   if(sidePosition == "left") {
      zAngle *= -1;
	}

   auto* choppingPara = new G4Para("choppingPara", choppingHalfLengthX,
                                     choppingHalfLengthY, choppingHalfLengthZ,
                                     yAngle, zAngle, xAngle);

   G4ThreeVector movePara(((longerLengthX - shorterLengthX) / 2.0 + shorterLengthX) / 2.0 + choppingHalfLengthX - choppingHalfLengthX * cos(phi), 0.0, 0.0);

   auto* rightExtension = new G4SubtractionSolid("rightExtension", uncutExtension, choppingPara, nullptr, movePara);

   return rightExtension;
}   // end ::sideSuppressorExtension()

G4Trap* DetectionSystemGriffin::SideSuppressorExtensionUncut() const
{

   G4double thicknessZ = fSuppressorExtensionThickness;
   G4double lengthY    = fSuppressorExtensionLength;

   G4double longerLengthX = ((fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness) / tan(fBentEndAngle) - fSuppressorExtensionThickness * sin(fBentEndAngle)) * tan(fBentEndAngle) - fBGOCanSeperation * 2.0);   // - fBGOCanSeperation

   G4double shorterLengthX = ((fSuppressorForwardRadius + fHevimetTipThickness) * sin(fBentEndAngle) - fBGOCanSeperation * 2.0);   // - fBGOCanSeperation

   // Similar to the "true" choppingSuppressor, add a bit on the z and y lengths to give a bit of a gap.
   thicknessZ += 0.1 * mm;
   lengthY += 0.1 * mm;

   auto* uncutExtension = new G4Trap("uncutExtension", thicknessZ, lengthY, longerLengthX, shorterLengthX);

   return uncutExtension;
}   // end ::sideSuppressorExtensionUncut()

///////////////////////////////////////////////////////////////////////
// methods used in ConstructNewHeavyMet
///////////////////////////////////////////////////////////////////////
G4SubtractionSolid* DetectionSystemGriffin::NewHeavyMet() const
{

   G4double halfThicknessZ = ((fHevimetTipThickness / cos(fHevimetTipAngle)) * cos(fBentEndAngle - fHevimetTipAngle) + fSuppressorForwardRadius * tan(fHevimetTipAngle) * sin(fBentEndAngle)) / 2.0;

   G4double halfLengthShorterX = fSuppressorForwardRadius * sin(fBentEndAngle);
   G4double halfLengthLongerX  = halfLengthShorterX + 2.0 * halfThicknessZ * tan(fBentEndAngle);
   G4double halfLengthShorterY = halfLengthShorterX;
   G4double halfLengthLongerY  = halfLengthLongerX;

   auto* uncutHevimet = new G4Trd("uncutHevimet", halfLengthLongerX,
                                   halfLengthShorterX, halfLengthLongerY, halfLengthShorterY, halfThicknessZ);

   G4double halfShorterLengthX = halfLengthShorterX - fSuppressorForwardRadius * tan(fHevimetTipAngle) * cos(fBentEndAngle);

   G4double halfHeightZ = (2.0 * halfThicknessZ + (halfLengthLongerX - ((fSuppressorForwardRadius + fHevimetTipThickness) * tan(fHevimetTipAngle)) / cos(fBentEndAngle) - halfShorterLengthX) * tan(fBentEndAngle)) / 2.0;

   G4double halfLongerLengthX  = halfShorterLengthX + 2.0 * halfHeightZ / tan(fBentEndAngle);
   G4double halfShorterLengthY = halfShorterLengthX;
   G4double halfLongerLengthY  = halfLongerLengthX;

   auto* intersector = new G4Trd("intersector", halfShorterLengthX,
                                  halfLongerLengthX, halfShorterLengthY, halfLongerLengthY, halfHeightZ);

   auto* chopper = new G4Trd("chopper", halfShorterLengthX, halfLongerLengthX,
                              halfShorterLengthY, halfLongerLengthY, halfHeightZ);

   G4ThreeVector moveChopper(0.0, 0.0, halfHeightZ + halfThicknessZ - fSuppressorForwardRadius * tan(fHevimetTipAngle) * sin(fBentEndAngle));

   auto* choppedHevimet = new G4SubtractionSolid("choppedHevimet",
                                                               uncutHevimet, chopper, nullptr, moveChopper);

   G4ThreeVector moveIntersector(0.0, 0.0, -halfHeightZ + halfThicknessZ);

   auto* intersectedHevimet = new G4IntersectionSolid("intersectedHevimet",
                                                                     choppedHevimet, intersector, nullptr, moveIntersector);

   G4double longerHalfLength = halfLengthLongerX - ((fSuppressorForwardRadius + fHevimetTipThickness) * tan(fHevimetTipAngle)) / cos(fBentEndAngle);

   G4double shorterHalfLength = longerHalfLength - 2.0 * (halfThicknessZ + 0.1 * mm) * tan(fBentEndAngle - fHevimetTipAngle);

   auto* middleChopper = new G4Trd("middleChopper", longerHalfLength, shorterHalfLength,
                                    longerHalfLength, shorterHalfLength, halfThicknessZ + 0.1 * mm);

   auto* hevimet = new G4SubtractionSolid("hevimet", intersectedHevimet, middleChopper);

   return hevimet;
}   //end ::newHeavyMet
