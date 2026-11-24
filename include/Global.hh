#ifndef GLOBAL_HH
#define GLOBAL_HH

#include "G4Types.hh"

inline G4double TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi)
{
   return (std::pow(x * x + y * y + z * z, 0.5) * std::sin(theta) * std::cos(phi));
}

inline G4double TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi)
{
   return (std::pow(x * x + y * y + z * z, 0.5) * std::sin(theta) * std::sin(phi));
}

inline G4double TransZ(G4double x, G4double y, G4double z, G4double theta)
{
   return (std::pow(x * x + y * y + z * z, 0.5) * std::cos(theta));
}

#endif
