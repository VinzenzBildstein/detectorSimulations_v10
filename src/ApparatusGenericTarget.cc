#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ApparatusGenericTarget.hh"

ApparatusGenericTarget::ApparatusGenericTarget()
   : fTargetMaterial("Gold")
{}

ApparatusGenericTarget::~ApparatusGenericTarget()
{
   delete fTargetLog;
	delete fAssembly;
}

void ApparatusGenericTarget::Build(const G4String& targetMaterialIn, G4double targetLengthXIn, G4double targetLengthYIn, G4double targetLengthZIn)
{
   fTargetMaterial = targetMaterialIn;
   fTargetLengthX  = targetLengthXIn * mm;
   fTargetLengthY  = targetLengthYIn * mm;
   fTargetLengthZ  = targetLengthZIn * mm;

   // Build assembly volume
   fAssembly = new G4AssemblyVolume();

   G4cout << "BuildTargetVolume" << G4endl;
   BuildTargetVolume();
}

void ApparatusGenericTarget::PlaceApparatus(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate)
{
   G4int copy_ID = 0;

   fAssembly->MakeImprint(exp_hall_log, move, rotate, copy_ID);
}

void ApparatusGenericTarget::BuildTargetVolume()
{
   G4Material* material = G4Material::GetMaterial(fTargetMaterial);
   if(material == nullptr) {
      G4cout << " ----> Material " << fTargetMaterial << " not found, cannot build!" << G4endl;
      throw std::invalid_argument("Material not found, cannot build!");
   }

   // Set visualization attributes
   auto* vis_att = new G4VisAttributes(G4Colour(0.8, 0.0, 0.2));
   vis_att->SetVisibility(true);

   G4ThreeVector     move;
   auto* rotate = new G4RotationMatrix;

   G4Box* target = BuildTarget();

   //logical volume
   if(fTargetLog == nullptr) {
      fTargetLog = new G4LogicalVolume(target, material, "target_log", nullptr, nullptr, nullptr);
      fTargetLog->SetVisAttributes(vis_att);
   }

   fAssembly->AddPlacedVolume(fTargetLog, move, rotate);
}

G4Box* ApparatusGenericTarget::BuildTarget() const
{
   return new G4Box("target", fTargetLengthX / 2., fTargetLengthY / 2., fTargetLengthZ / 2.);
}

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector ApparatusGenericTarget::GetDirectionXYZ(G4double theta, G4double phi)
{
   G4double x = sin(theta) * cos(phi);
   G4double y = sin(theta) * sin(phi);
   G4double z = cos(theta);

   return {x, y, z};
}
