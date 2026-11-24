#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

//#include "FieldSetup.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4Material.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemGrid.hh"

#include "G4SystemOfUnits.hh"   // new version geant4.10 requires units

DetectionSystemGrid::DetectionSystemGrid(G4double xLengthIn, G4double yLengthIn, G4double zLengthIn, G4double cubeSizeIn, G4String gridMat, G4ThreeVector gridColour, const G4ThreeVector& posOffset)
   : fXLength(xLengthIn), fYLength(yLengthIn), fZLength(zLengthIn),
     fXPosOffset(posOffset.x()), fYPosOffset(posOffset.y()), fZPosOffset(posOffset.z()),
     fCubeSize(cubeSizeIn),
     fCellsXRow(static_cast<G4int>(fXLength / fCubeSize)), fCellsYRow(static_cast<G4int>(fYLength / fCubeSize)), fCellsZRow(static_cast<G4int>(fZLength / fCubeSize)),
     fNumberOfCells(fCellsXRow * fCellsYRow * fCellsZRow), fCellsPerRow(10), fCellWidth(fCubeSize),
     fCellMaterial(std::move(gridMat)), fCellColour(std::move(gridColour))
{
   //fEmFieldSetup = theFieldSetup;

   G4cout << "xLength = " << fXLength / cm << G4endl;
   G4cout << "yLength = " << fYLength / cm << G4endl;
   G4cout << "zLength = " << fZLength / cm << G4endl;

   G4cout << "cellsXRow = " << fCellsXRow << G4endl;
}

DetectionSystemGrid::~DetectionSystemGrid()
{
   // LogicalVolumes
   delete fGridcellLog;
}

G4int DetectionSystemGrid::Build()
{

   // Build assembly volume
   fAssembly = new G4AssemblyVolume();

   G4cout << "BuildCellVolume" << G4endl;
   BuildCellVolume();

   return 1;
}

G4int DetectionSystemGrid::PlaceDetector(G4LogicalVolume* expHallLog)
{
   G4double startX = (-1.0 * fXLength / 2.0);
   G4double startY = (-1.0 * fYLength / 2.0);
   G4double startZ = (-1.0 * fZLength / 2.0);

   G4int cellNumber = 0;

   G4ThreeVector     position = G4ThreeVector(0.0 * mm, 0.0 * mm, 0.0 * mm);
   auto* rotate   = new G4RotationMatrix;

   for(G4int rowX = 0; rowX < fCellsXRow; rowX++) {
      G4double posX = startX + (fCellWidth * rowX) + (fCellWidth / 2.0);
      position.setX(posX + fXPosOffset);

      for(G4int rowY = 0; rowY < fCellsYRow; rowY++) {
         G4double posY = startY + (fCellWidth * rowY) + (fCellWidth / 2.0);
         position.setY(posY + fYPosOffset);

         for(G4int rowZ = 0; rowZ < fCellsZRow; rowZ++) {
            G4double posZ = startZ + (fCellWidth * rowZ) + (fCellWidth / 2.0);
            position.setZ(posZ + fZPosOffset);

            cellNumber++;
            //G4cout<<"--------------------------- Grid Detector Number = "<<cellNumber<<G4endl;

            fAssembly->MakeImprint(expHallLog, position, rotate, cellNumber);
         }
      }
   }

   return 1;
}

G4int DetectionSystemGrid::BuildCellVolume()
{
   G4Material* material = G4Material::GetMaterial(fCellMaterial);
   if(material == nullptr) {
      G4cout << " ----> Material " << fCellMaterial << " not found, cannot build the detector shell! " << G4endl;
      return 0;
   }

   // Set visualization attributes
   auto* visAtt = new G4VisAttributes(G4Colour(fCellColour));
   visAtt->SetVisibility(true);

   G4Box* gridcell = BuildCell();

   // Define rotation and movement objects
   G4ThreeVector     direction = G4ThreeVector(0, 0, 0);
   G4double          zPosition = 0.0 * mm;
   G4ThreeVector     move      = zPosition * direction;
   auto* rotate    = new G4RotationMatrix;

   //logical volume
   if(fGridcellLog == nullptr) {
      fGridcellLog = new G4LogicalVolume(gridcell, material, "gridcellLog", nullptr, nullptr, nullptr);
      fGridcellLog->SetVisAttributes(visAtt);

      //// Set local field manager and local field in radiator and its daughters:
      //G4bool allLocal = true ;
      //gridcellLog->SetFieldManager(fEmFieldSetup->GetLocalFieldManager(), allLocal) ;
   }

   fAssembly->AddPlacedVolume(fGridcellLog, move, rotate);

   return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGrid::BuildCell() const
{
   G4double halfLengthX = fCellWidth / 2.0;
   G4double halfLengthY = fCellWidth / 2.0;
   G4double halfLengthZ = fCellWidth / 2.0;

   return new G4Box("cell", halfLengthX, halfLengthY, halfLengthZ);
}   //end ::BuildCell
