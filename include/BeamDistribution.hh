#ifndef BEAMDISTRIBUTION_HH
#define BEAMDISTRIBUTION_HH

#include "globals.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>
#include <string>

class BeamDistribution {
public:
   BeamDistribution() = default;
   explicit BeamDistribution(const G4String& filename);
   ~BeamDistribution() = default;

   void     LoadDistribution(const G4String& filename);
   G4double GetRandom();

   bool Good() const { return fLoadGood; }

private:
   bool                  fLoadGood{false};
   std::vector<G4double> fDistVec;      //these vectors will hold all of the data that is read in
   std::vector<G4double> fRelProbVec;   //these vectors will hold all of the data that is read in
   G4double              fDistStep{0.};
};
#endif
