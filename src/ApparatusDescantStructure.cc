#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemDescant.hh"
#include "ApparatusDescantStructure.hh"

#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"

#include "G4VSolid.hh"
#include "G4Polyhedra.hh"

#include <vector>

ApparatusDescantStructure::ApparatusDescantStructure()
   // DESCANT structure properties
	// using freeCAD to calculate thickness of the shell from STEP file
	: fStructureMaterial("G4_Al"), fStructureInnerRadius(915.82 * mm), fStructureOuterRadius(1007.9566 * mm),
   fStartPhi(0. * deg), fDeltaPhi((360. / 5) * deg), fStartTheta(13. * deg), fDeltaTheta(67.8743593 * deg - fStartTheta),
   fStructureRadialDistance(85.75 * cm),
   fNumSides(6), fNumZplanes(2),
   fInnerRadiusRed({{0. * mm, 0. * mm}}),
   fZPlanes({{-500. * mm, 500. * mm}}),
   fInnerRadiusWhite({{0. * mm, 0. * mm}}),
   fInnerRadiusGreenYellow({{0. * mm, 0. * mm}}),
   fInnerRadiusBack({{0. * mm, 0. * mm}}),
   fOuterRadiusBack({{150. * mm, 150. * mm}}),
   fZPlanesBack({{-9. * mm, 9. * mm}}),
   // for cylinder subtraction from back plates
   fMinRadius(0. * mm), fMaxRadius(50. * mm), fHalfLength(2000. * mm), fSubtractionStartPhi(0. * deg), fSubtractionDeltaPhi(360. * deg),
   // for subtraction box
   fSubBoxX(1000. * mm), fSubBoxY(1000. * mm), fSubBoxZ(1000. * mm), fMoveCut(100. * mm)

{
   // for red polyhedra subtraction
	G4double cutExtra = 4. * mm;
	fPolyRadiusRed = 171.57 * mm / 2. + cutExtra;   // taken from furthest distance between opposide sides on a red detector
   fOuterRadiusRed = {{fPolyRadiusRed, fPolyRadiusRed}};

   // for green/yellow polyhedra subtraction
   fPolyRadiusGreenYellow = 150.96 * mm / 2. + cutExtra;   // taken from furthest distance between opposide sides on a red detector
   fOuterRadiusGreenYellow = {{fPolyRadiusGreenYellow, fPolyRadiusGreenYellow}};

   // for white polyhedra subtraction
   cutExtra         = 10. * mm;
   fPolyRadiusWhite = 150.86 * mm / 2. + cutExtra;
   fOuterRadiusWhite = {{fPolyRadiusWhite, fPolyRadiusWhite}};

   // The Euler angles from James' MSc thesis which gives us the detector positions
   // Some of the angles for the green and yellow detectors are wrong in James' thesis,
   // note the +180 on a few angles.

   fBlueAlphaBetaGamma = {{
		{{-22.39 * deg, -33.76 * deg, -71.19 * deg}},
      {{-49.61 * deg, -33.76 * deg, 71.19 * deg}},
      {{-36.00 * deg, -46.06 * deg, 180.00 * deg}},
      {{85.61 * deg, 33.76 * deg, 108.81 * deg}},
      {{58.39 * deg, 33.76 * deg, 251.19 * deg}},
      {{72.00 * deg, 46.06 * deg, 360.00 * deg}},
      {{13.61 * deg, 33.76 * deg, 108.81 * deg}},
      {{-13.61 * deg, 33.76 * deg, 251.19 * deg}},
      {{0.00 * deg, 46.06 * deg, 360.00 * deg}},
      {{-58.39 * deg, 33.76 * deg, 108.81 * deg}},
      {{-85.61 * deg, 33.76 * deg, 251.19 * deg}},
      {{-72.00 * deg, 46.06 * deg, 360.00 * deg}},
      {{49.61 * deg, -33.76 * deg, -71.19 * deg}},
      {{22.39 * deg, -33.76 * deg, 71.19 * deg}},
      {{36.00 * deg, -46.06 * deg, 180.00 * deg}}
	}};

   fGreenAlphaBetaGamma = {{
		{{-12.45 * deg, -60.47 * deg, -12.45 * deg}},
      {{-27.73 * deg, -58.55 * deg, -4.54 * deg}},
      {{-84.45 * deg, -60.47 * deg, -12.45 * deg}},
      {{80.27 * deg, 58.55 * deg, (-4.54 + 180) * deg}},
      {{23.55 * deg, 60.47 * deg, 167.55 * deg}},
      {{8.27 * deg, 58.55 * deg, 175.46 * deg}},
      {{-48.45 * deg, 60.47 * deg, 167.55 * deg}},
      {{-63.73 * deg, 58.55 * deg, 175.46 * deg}},
      {{59.55 * deg, -60.47 * deg, -12.45 * deg}},
      {{44.27 * deg, -58.55 * deg, -4.54 * deg}}
	}};

   fRedAlphaBetaGamma = {{
		{{-36.00 * deg, -20.39 * deg, 180.00 * deg}},
      {{-16.04 * deg, -47.84 * deg, 45.10 * deg}},
      {{-55.96 * deg, -47.84 * deg, 314.90 * deg}},
      {{72.00 * deg, 20.39 * deg, 0.00 * deg}},
      {{-88.04 * deg, -47.84 * deg, 45.10 * deg}},
      {{52.04 * deg, 47.84 * deg, 134.90 * deg}},
      {{0.00 * deg, 20.39 * deg, 0.00 * deg}},
      {{19.96 * deg, 47.84 * deg, 225.10 * deg}},
      {{-19.96 * deg, 47.84 * deg, 134.90 * deg}},
      {{-72.00 * deg, 20.39 * deg, 0.00 * deg}},
      {{-52.04 * deg, 47.84 * deg, 225.10 * deg}},
      {{88.04 * deg, -47.84 * deg, 314.90 * deg}},
      {{36.00 * deg, -20.39 * deg, 180.00 * deg}},
      {{55.96 * deg, -47.84 * deg, 45.10 * deg}},
      {{16.04 * deg, -47.84 * deg, 314.90 * deg}}
	}};

   fWhiteAlphaBetaGamma = {{
		{{0.00 * deg, -11.37 * deg, 90.00 * deg}},
      {{0.00 * deg, -24.67 * deg, 90.00 * deg}},
      {{0.00 * deg, -38.76 * deg, -90.00 * deg}},
      {{0.00 * deg, -52.06 * deg, -90.00 * deg}},
      {{-72.00 * deg, -11.37 * deg, 90.00 * deg}},
      {{-72.00 * deg, -24.67 * deg, 90.00 * deg}},
      {{-72.00 * deg, -38.76 * deg, -90.00 * deg}},
      {{-72.00 * deg, -52.06 * deg, -90.00 * deg}},
      {{36.00 * deg, 11.37 * deg, 270.00 * deg}},
      {{36.00 * deg, 24.67 * deg, 270.00 * deg}},
      {{36.00 * deg, 38.76 * deg, 90.00 * deg}},
      {{36.00 * deg, 52.06 * deg, 90.00 * deg}},
      {{-36.00 * deg, 11.37 * deg, 270.00 * deg}},
      {{-36.00 * deg, 24.67 * deg, 270.00 * deg}},
      {{-36.00 * deg, 38.76 * deg, 90.00 * deg}},
      {{-36.00 * deg, 52.06 * deg, 90.00 * deg}},
      {{72.00 * deg, -11.37 * deg, 90.00 * deg}},
      {{72.00 * deg, -24.67 * deg, 90.00 * deg}},
      {{72.00 * deg, -38.76 * deg, -90.00 * deg}},
      {{72.00 * deg, -52.06 * deg, -90.00 * deg}}   //59
	}};

   fYellowAlphaBetaGamma = {{
		{{-44.27 * deg, -58.55 * deg, (4.54 + 180) * deg}},
      {{-59.55 * deg, -60.47 * deg, 192.45 * deg}},
      {{63.73 * deg, 58.55 * deg, 4.54 * deg}},
      {{48.45 * deg, 60.47 * deg, (192.45 + 180) * deg}},
      {{-8.27 * deg, 58.55 * deg, (184.54 + 180) * deg}},
      {{-23.55 * deg, 60.47 * deg, 372.45 * deg}},
      {{-80.27 * deg, 58.55 * deg, (184.54 + 180) * deg}},
      {{84.45 * deg, -60.47 * deg, (372.45 + 180) * deg}},
      {{27.73 * deg, -58.55 * deg, (4.54 + 180) * deg}},
      {{12.45 * deg, -60.47 * deg, 192.45 * deg}}
	}};

   fGreyColour = G4Colour(0.9, 0.9, 0.9);
}

ApparatusDescantStructure::~ApparatusDescantStructure()
{
   // LogicalVolumes
   delete fDescantStructureLog;
}

G4int ApparatusDescantStructure::Build()
{
   fAssemblyDescantStructure = new G4AssemblyVolume();

   BuildDescantShell();

   return 1;
}

G4int ApparatusDescantStructure::PlaceDescantStructure(G4LogicalVolume* expHallLog)
{
   G4ThreeVector move = G4ThreeVector(0., 0., 0.);

   for(G4int i = 0; i < 5; i++) {
      auto* rotate = new G4RotationMatrix;
      rotate->rotateZ((i * 72.) * deg);
      fAssemblyDescantStructure->MakeImprint(expHallLog, move, rotate);
   }

   return 0;
}

G4int ApparatusDescantStructure::BuildDescantShell()
{
   G4ThreeVector move           = G4ThreeVector(0., 0., 0.);
   auto*         rotate         = new G4RotationMatrix;
   G4int         detectorNumber = 70;
   G4int         idx            = 0;
   fStructureRadialDistance     = fStructureOuterRadius + 8. * mm;

   // for cutting the detector shapes
   auto* subtractionBox = new G4Box("subtractionBox", fSubBoxX, fSubBoxY, fSubBoxZ);

   // for construction of main DESCANT shell structure
   auto* alumSphere = new G4Sphere("descant sphere",
                                       fStructureInnerRadius, fStructureOuterRadius, fStartPhi, fDeltaPhi, fStartTheta, fDeltaTheta);

   // for the subtraction of each detector shape from the solid shell
   //red - cut
   auto* polyRed = new G4Polyhedra("polyRed", fStartPhi, fDeltaPhi, fNumSides, fNumZplanes, fZPlanes.data(), fInnerRadiusRed.data(), fOuterRadiusRed.data());
   rotate               = new G4RotationMatrix;
   rotate->rotateZ(20. * deg);
   fMoveCut = 120. * mm;
   move     = G4ThreeVector(fSubBoxX + fMoveCut, 0., 0.);
   rotate   = new G4RotationMatrix;
   rotate->rotateZ(-20. * deg);
   fMoveCut = 120. * mm;
   move     = G4ThreeVector(fSubBoxX + fMoveCut, 0., 0.);
   //white
   auto* polyWhite = new G4Polyhedra("polyWhite", fStartPhi, fDeltaPhi, fNumSides, fNumZplanes, fZPlanes.data(), fInnerRadiusWhite.data(), fOuterRadiusWhite.data());
   //green or yellow before cut
   auto* polyGreenYellow = new G4Polyhedra("polyGreenYellow", fStartPhi, fDeltaPhi, fNumSides, fNumZplanes, fZPlanes.data(), fInnerRadiusGreenYellow.data(), fOuterRadiusGreenYellow.data());
   //green after cut
   rotate                           = new G4RotationMatrix;
   fMoveCut                         = 65. * mm;
   move                             = G4ThreeVector(fSubBoxX + fMoveCut, 0., 0.);
   auto* polyGreenCut = new G4SubtractionSolid("polyGreenCut", polyGreenYellow, subtractionBox, rotate, move);
   //yellow after cut
   rotate                            = new G4RotationMatrix;
   fMoveCut                          = 65. * mm;
   move                              = G4ThreeVector(-fSubBoxX - fMoveCut, 0., 0.);
   auto* polyYellowCut = new G4SubtractionSolid("polyYellowCut", polyGreenYellow, subtractionBox, rotate, move);

   // for construction of the back plates of each DESCANT detector
   // general
   rotate                            = new G4RotationMatrix;
   move                              = G4ThreeVector(0., 0., 0.);
   auto*        shellBack     = new G4Polyhedra("shellBack", fStartPhi, fDeltaPhi, fNumSides, fNumZplanes, fZPlanesBack.data(), fInnerRadiusBack.data(), fOuterRadiusBack.data());
   auto*             cylinder      = new G4Tubs("cylinder", fMinRadius, fMaxRadius, fHalfLength, fSubtractionStartPhi, fSubtractionDeltaPhi);
   auto* shellBackHole = new G4SubtractionSolid("shellBackWithHoles", shellBack, cylinder, rotate, move);
   // yellow
   fMoveCut                            = 90. * mm;
   rotate                              = new G4RotationMatrix;
   move                                = G4ThreeVector(-fSubBoxX - fMoveCut, 0., 0.);
   auto* shellBackYellow = new G4SubtractionSolid("shellBackWithHolesYellow", shellBackHole, subtractionBox, rotate, move);
   // green
   fMoveCut                           = 90. * mm;
   rotate                             = new G4RotationMatrix;
   move                               = G4ThreeVector(fSubBoxX + fMoveCut, 0., 0.);
   auto* shellBackGreen = new G4SubtractionSolid("shellBackWithHolesGreen", shellBackHole, subtractionBox, rotate, move);

   // Colour
   auto* greyVisAtt = new G4VisAttributes(fGreyColour);
   greyVisAtt->SetVisibility(true);

   // Material
   G4Material* structureG4material = G4Material::GetMaterial(fStructureMaterial);
   if(structureG4material == nullptr) {
      G4cout << " ----> Material " << fStructureMaterial << " not found, cannot build! " << G4endl;
      return 0;
   }

   auto* logicBack = new G4LogicalVolume(shellBackHole, structureG4material, "logicBack", nullptr, nullptr, nullptr);
   logicBack->SetVisAttributes(fGreyColour);

   auto* logicBackYellow = new G4LogicalVolume(shellBackYellow, structureG4material, "logicBackYellow", nullptr, nullptr, nullptr);
   logicBackYellow->SetVisAttributes(fGreyColour);

   auto* logicBackGreen = new G4LogicalVolume(shellBackGreen, structureG4material, "logicBackGreen", nullptr, nullptr, nullptr);
   logicBackYellow->SetVisAttributes(fGreyColour);

   // for making the DESCANT shell
   G4SubtractionSolid*              preSubtraction = nullptr;
   G4SubtractionSolid*              postSubtraction = nullptr;
   std::vector<G4SubtractionSolid*> descantShell;

   // subtract Blue Detector
   for(G4int i = 0; i < (detectorNumber - 55); i++) {
      idx = i;
      //
      move = G4ThreeVector(0., 0., fStructureRadialDistance);
      move.rotateZ(fBlueAlphaBetaGamma[idx][2]);
      move.rotateY(fBlueAlphaBetaGamma[idx][1]);
      move.rotateZ(fBlueAlphaBetaGamma[idx][0]);
      rotate = new G4RotationMatrix;
      rotate->rotateY(M_PI);   // flip the detector so that the face is pointing upstream.
      rotate->rotateZ(fBlueAlphaBetaGamma[idx][2]);
      rotate->rotateY(fBlueAlphaBetaGamma[idx][1]);
      rotate->rotateZ(fBlueAlphaBetaGamma[idx][0]);
      if(i == 4 || i == 6 || i == 8) { fAssemblyDescantStructure->AddPlacedVolume(logicBack, move, rotate); }
      rotate->invert();
      std::stringstream temp;
      temp << "blueSubtraction_" << i;
      G4String name = temp.str();
      if(i == 4 || i == 5 || i == 6 || i == 8) {
         if(i == 4) {
            postSubtraction = new G4SubtractionSolid(name, alumSphere, polyWhite, rotate, move);
            descantShell.push_back(postSubtraction);
         } else {
            preSubtraction  = descantShell.back();
            postSubtraction = new G4SubtractionSolid(name, preSubtraction, polyWhite, rotate, move);
            descantShell.push_back(postSubtraction);
         }
      }
   }

   // subtract Green Detector
   for(G4int i = 15; i < (detectorNumber - 45); i++) {
      idx  = i - 15;
      move = G4ThreeVector(0., 0., fStructureRadialDistance);
      move.rotateZ(fGreenAlphaBetaGamma[idx][2]);
      move.rotateY(fGreenAlphaBetaGamma[idx][1]);
      move.rotateZ(fGreenAlphaBetaGamma[idx][0]);
      rotate = new G4RotationMatrix;
      rotate->rotateY(M_PI);   // flip the detector so that the face is pointing upstream.
      rotate->rotateZ(fGreenAlphaBetaGamma[idx][2]);
      rotate->rotateY(fGreenAlphaBetaGamma[idx][1]);
      rotate->rotateZ(fGreenAlphaBetaGamma[idx][0]);
      if(i == 19 || i == 20) { fAssemblyDescantStructure->AddPlacedVolume(logicBackGreen, move, rotate); }
      rotate->invert();
      std::stringstream temp;
      temp << "greenSubtraction_" << i;
      G4String name = temp.str();
      if(i == 19 || i == 20) {
         preSubtraction  = descantShell.back();
         postSubtraction = new G4SubtractionSolid(name, preSubtraction, polyGreenCut, rotate, move);
         descantShell.push_back(postSubtraction);
      }
   }

   // subtract Red Detector
   for(G4int i = 25; i < (detectorNumber - 30); i++) {
      idx  = i - 25;
      move = G4ThreeVector(0., 0., fStructureRadialDistance);
      move.rotateZ(fRedAlphaBetaGamma[idx][2]);
      move.rotateY(fRedAlphaBetaGamma[idx][1]);
      move.rotateZ(fRedAlphaBetaGamma[idx][0]);
      rotate = new G4RotationMatrix;
      rotate->rotateY(M_PI);   // flip the detector so that the face is pointing upstream.
      rotate->rotateZ(fRedAlphaBetaGamma[idx][2]);
      rotate->rotateY(fRedAlphaBetaGamma[idx][1]);
      rotate->rotateZ(fRedAlphaBetaGamma[idx][0]);
      if(i == 30 || i == 31 || i == 32) { fAssemblyDescantStructure->AddPlacedVolume(logicBack, move, rotate); }
      rotate->invert();
      std::stringstream temp;
      temp << "redSubtraction_" << i;
      G4String name = temp.str();
      if(i == 28 || i == 30 || i == 31 || i == 32) {
         preSubtraction  = descantShell.back();
         postSubtraction = new G4SubtractionSolid(name, preSubtraction, polyRed, rotate, move);
         //postSubtraction = new G4SubtractionSolid(name, preSubtraction, polyRedCut2, rotate, move);
         descantShell.push_back(postSubtraction);
      }
   }

   // subtract White Detector
   for(G4int i = 40; i < (detectorNumber - 10); i++) {
      idx  = i - 40;
      move = G4ThreeVector(0., 0., fStructureRadialDistance);
      move.rotateZ(fWhiteAlphaBetaGamma[idx][2]);
      move.rotateY(fWhiteAlphaBetaGamma[idx][1]);
      move.rotateZ(fWhiteAlphaBetaGamma[idx][0]);
      rotate = new G4RotationMatrix;
      rotate->rotateY(M_PI);   // flip the detector so that the face is pointing upstream.
      rotate->rotateZ(fWhiteAlphaBetaGamma[idx][2]);
      rotate->rotateY(fWhiteAlphaBetaGamma[idx][1]);
      rotate->rotateZ(fWhiteAlphaBetaGamma[idx][0]);
      if(i == 48 || i == 49 || i == 50 || i == 51) { fAssemblyDescantStructure->AddPlacedVolume(logicBack, move, rotate); }
      rotate->invert();
      std::stringstream temp;
      temp << "whiteSubtraction_" << i;
      G4String name = temp.str();
      if(i == 48 || i == 49 || i == 50 || i == 51) {
         preSubtraction  = descantShell.back();
         postSubtraction = new G4SubtractionSolid(name, preSubtraction, polyWhite, rotate, move);
         descantShell.push_back(postSubtraction);
      }
   }
   // subtract Yellow Detector
   for(G4int i = 60; i < (detectorNumber); i++) {
      idx  = i - 60;
      move = G4ThreeVector(0., 0., fStructureRadialDistance);
      move.rotateZ(fYellowAlphaBetaGamma[idx][2]);
      move.rotateY(fYellowAlphaBetaGamma[idx][1]);
      move.rotateZ(fYellowAlphaBetaGamma[idx][0]);
      rotate = new G4RotationMatrix;
      rotate->rotateY(M_PI);   // flip the detector so that the face is pointing upstream.
      rotate->rotateZ(fYellowAlphaBetaGamma[idx][2]);
      rotate->rotateY(fYellowAlphaBetaGamma[idx][1]);
      rotate->rotateZ(fYellowAlphaBetaGamma[idx][0]);
      if(i == 62 || i == 63) { fAssemblyDescantStructure->AddPlacedVolume(logicBackYellow, move, rotate); }
      rotate->invert();
      std::stringstream temp;
      temp << "yellowSubtraction_" << i;
      G4String name = temp.str();
      if(i == 62 || i == 63) {
         preSubtraction  = descantShell.back();
         postSubtraction = new G4SubtractionSolid(name, preSubtraction, polyYellowCut, rotate, move);
         descantShell.push_back(postSubtraction);
      }
   }

   // Logical Volume
   if(fDescantStructureLog == nullptr) {
      fDescantStructureLog = new G4LogicalVolume(descantShell.back(), structureG4material, "descantStructureLog", nullptr, nullptr, nullptr);
      fDescantStructureLog->SetVisAttributes(greyVisAtt);
   }

   rotate = new G4RotationMatrix;
   move   = G4ThreeVector(0., 0., 0.);
   fAssemblyDescantStructure->AddPlacedVolume(fDescantStructureLog, move, rotate);

   return 1;
}
