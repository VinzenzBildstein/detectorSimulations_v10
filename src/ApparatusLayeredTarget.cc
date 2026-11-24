#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"   //individual step-size per logical volume/region

#include "ApparatusLayeredTarget.hh"

ApparatusLayeredTarget::ApparatusLayeredTarget(G4double beaminput)
   : fBeamPos(beaminput * mm), fTargetRadius(6.0 * mm), fStepLimit(new G4UserLimits(10. * CLHEP::um))
{
}

ApparatusLayeredTarget::~ApparatusLayeredTarget()
{
   // LogicalVolumes
   for(auto& layer : fTargetLayerLog) {
      delete layer;
   }
   delete fStepLimit;
}

G4int ApparatusLayeredTarget::BuildTargetLayer(const G4String& LayerMaterial, G4double LayerArealDensity)
{
   //Get Materials
   G4Material* material = G4Material::GetMaterial(LayerMaterial);
   if(material == nullptr) {
      G4cout << " ----> Target Material " << LayerMaterial << " not found, cannot build!" << G4endl;
      return 0;
   }

   G4double layerMaterialDensity = material->GetDensity();
   G4cout << "Target LayerMaterialDensity: " << layerMaterialDensity / (g / cm3) << " g/cm3" << G4endl;
   G4double layerThickness = ((LayerArealDensity * (mg / cm2)) / layerMaterialDensity);
   G4cout << "Target Thickness (um): " << layerThickness / (um) << G4endl;

   // Build assembly volume
   auto* vis_att = new G4VisAttributes(G4Colour(0.8 + static_cast<double>(fTargetLayerLog.size()) * 0.1, 0.0, 0.2));
   vis_att->SetVisibility(true);

   std::stringstream ss;
   ss << "TargetLayer" << fTargetLayerLog.size();
   auto* target = new G4Tubs(ss.str(), 0.0, fTargetRadius, std::abs(layerThickness) / 2., 0.0, 360 * deg);

   //logical volume
   ss.str("");
	ss.clear();
   ss << "TargetLayerLog" << fTargetLayerLog.size();
   auto* layerLog = new G4LogicalVolume(target, material, ss.str(), nullptr, nullptr, nullptr);
   layerLog->SetVisAttributes(vis_att);
   fTargetLayerLog.push_back(layerLog);
   fTargetLayerThick.push_back(layerThickness);

   layerLog->SetUserLimits(fStepLimit);
   return 1;
}

void ApparatusLayeredTarget::PlaceTarget(G4LogicalVolume* oexpHallLog)
{
   if(fTargetLayerLog.size() <= fTargetLayerPhys.size()) { return; }
   G4cout << "Adding Target Layer " << G4endl;

   G4double placement = LayerStart(fTargetLayerPhys.size()) - fTargetLayerThick[fTargetLayerPhys.size()] / 2.;
   G4cout << "Placement target : " << placement << " " << placement / mm << G4endl;

   G4ThreeVector     move    = G4ThreeVector(0., 0., placement);   //0 for now until read in
   auto* sRotate = new G4RotationMatrix;

   std::stringstream ss;
   ss << "TargetLayer" << fTargetLayerLog.size();
   auto* TargetPhys = new G4PVPlacement(sRotate,                                    //no rotation
                                        move,                                       //at position
                                        fTargetLayerLog[fTargetLayerPhys.size()],   //its logical volume
                                        ss.str(),                                   //its name
                                        oexpHallLog,                                //its mother  volume
                                        false,                                      //no boolean operation
                                        0,                                          //copy number
                                        false);                                     //overlaps checking

   fTargetLayerPhys.push_back(TargetPhys);
}

G4double ApparatusLayeredTarget::LayerStart(int layeri)
{
	if(fTargetLayerThick.empty()) { return 0; }

	G4double ret = fBeamPos;
	for(int i = 0; i < layeri; i++) {
		if(i < static_cast<int>(fTargetLayerThick.size())) { ret -= fTargetLayerThick[i]; }
	}

	return ret;
}
