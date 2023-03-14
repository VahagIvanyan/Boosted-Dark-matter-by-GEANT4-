
#ifndef __RAT_VertexGen_BDS__
#define __RAT_VertexGen_BDS__

#include <RAT/GLG4VertexGen.hh>
#include <CLHEP/Random/RandGeneral.h>
#include <G4TransportationManager.hh>
#include <G4NistManager.hh>
#include <G4ParticleTable.hh>
#include <RAT/DB.hh>


namespace RAT {

class VertexGen_BDS : public GLG4VertexGen {
public:
  VertexGen_BDS(const char *arg_dbname="bds");
  virtual ~VertexGen_BDS();
  virtual void GeneratePrimaryVertex( G4Event *argEvent,
                                      G4ThreeVector &dx,
                                      G4double dt);
  /** State format wil be added later for the kinematics */
        virtual void SetState( G4String newValues );
        virtual G4String GetState();


protected:
  void Setup();
  std::string fNucleusName;

    G4Navigator* theNavigator;
    G4NistManager* man;
    G4ParticleDefinition *fGeantino;
    G4double fmassGeantino; //GeV (the mass)

//G4ParticleDefinition *fElectron;// fro electron
//G4ParticleDEfinition *fPositron;// for electron and positron

 //CLHEP::RandGeneral *fEnergyDist;
 //double fEnergyDistLo, fEnergyDistHi;
};


} // namespace RAT

#endif


