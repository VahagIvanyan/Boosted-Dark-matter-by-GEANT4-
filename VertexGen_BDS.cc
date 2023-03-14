////
///     In the code dx and dt should be checked and understood
////
#include <RAT/VertexGen_BDS.hh>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <Randomize.hh>
#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4ThreeVector.hh>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <G4VPhysicalVolume.hh>
#include <G4TransportationManager.hh>
#include <G4Material.hh>

#include <G4Element.hh>
#include <G4NistManager.hh>

#include <ctime>


namespace RAT {

 VertexGen_BDS::VertexGen_BDS(const char *arg_dbname)
    : GLG4VertexGen(arg_dbname)
  {

  }

 VertexGen_BDS::~VertexGen_BDS()
  {

  }
  void VertexGen_BDS::GeneratePrimaryVertex(G4Event *event, G4ThreeVector &dx, G4double dt)
  {
    /**
     * Current production version of code used to implement the generator
    **/
   G4ThreeVector Position, decayPos;
   //G4cout << "1" << G4endl;
   G4ThreeVector Momentum,Momentum1,Momentum2;
   //G4cout << "1" << G4endl;

   G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");//->FindParticle("e-");// Result of the Recoiled Argon ionization with electron emission
   //G4cout << "1" << G4endl;
 G4ParticleDefinition* particle1 = G4ParticleTable::GetParticleTable()->FindParticle("geantino");//->FindParticle("e-");
   //G4cout << "1" << G4endl;
   G4ParticleDefinition* particle2 = G4ParticleTable::GetParticleTable()->FindParticle("geantino");//->FindParticle("e+");
   //G4cout << "1" << G4endl;
   //G4ParticleDefinition* fGeantino = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
 G4cout<<"Initial Position:\t"<<Position<<G4endl;
   //G4VPhysicalVolume* phisVol = theNavigator->LocateGlobalPointAndSetup(Position, (const G4ThreeVector*) 0, false, true); // g$
//DBLinkPtr dbMaterial = db->GetLink("MATERIAL", pv->GetLogicalVolume()->GetMaterial()->GetName()); // gets the properties of the material from MATERIALS.ratdb
   G4cout << "Checking: 1" << G4endl;

 //if (phisVol->GetName() == "cryoliquid")// if we are in LAr volume
 //{
   G4cout << "Checking: 1" << G4endl;

     //   G4cout<<phisVol->GetName()<<G4endl;
        ///energy of the elasticlly scattered electron        WHILE for now for geantino
        G4double energy = 1000*keV;
        G4double energy1 = 511*keV;
        G4double energy2 = 511*keV;

        // position should be randomly selected from the inside of the LAr (random selection of Radius value(1.7 meter))

        G4double PosX,PosY,PosZ;
        G4double bettaPosX, bettaPosY, bettaPosZ;
        G4double R_Sphere = 850*mm;// Would be better to find the R_Sphere from the "geometry file (DetectorConstruction or DB or wherever it is defined)" of the LAr detector
        PosX = G4UniformRand() * R_Sphere;
        PosY = G4UniformRand() * R_Sphere;
        PosZ = G4UniformRand() * R_Sphere;//should be selected randomly and no need to be equal each other
        bettaPosX = G4UniformRand() * R_Sphere;// X -> (e-,e+) e- and e+ randomly selected positions (but same for both partcile)
        bettaPosY =  G4UniformRand() * R_Sphere;
        bettaPosZ = G4UniformRand() * R_Sphere;

//&dx argument of the vertex generator must be somehow connected/assigned to PosX,PosY,PosZ and as its G4ThreeVector, then similar to:
   G4cout << "Checking: 1" << G4endl;
        Position = G4ThreeVector(PosX,PosY,PosZ);
        //Pos.setX(PosX); Pos.setY(PosY); Pos.setZ(PosZ);

        //dt = maybe would be better to connect with timer
        G4double firstScatStart = G4UniformRand()*1000000*ns;

        // Pick direction isotropically

        double theta = acos(2.0 * G4UniformRand() - 1.0);
        double phi = 2.0 * G4UniformRand() * pi;
        Momentum.setRThetaPhi(energy, theta, phi); // Momentum == energy units in GEANT4
        G4ThreeVector dir = Momentum.unit();
        dx = Position;
        dt = firstScatStart;
        G4PrimaryVertex* vertex= new G4PrimaryVertex(dx, dt);
        G4PrimaryParticle* geantino = new G4PrimaryParticle(particle,Momentum.x(),Momentum.y(),Momentum.z());
        // Generate random polarization
        phi = (G4UniformRand()*2.0-1.0)*M_PI;
       G4ThreeVector e1 = dir.orthogonal().unit();
        G4ThreeVector e2 = dir.cross(e1);
        G4ThreeVector rpol = e1*cos(phi)+e2*sin(phi);
        geantino->SetPolarization(rpol.x(), rpol.y(), rpol.z());

        // Something for e-,e+ decay
        //dx
        decayPos = G4ThreeVector(bettaPosX,bettaPosY,bettaPosZ);//setX(bettaPosX);decayPos.setY(bettaPosY);decayPos.setY(bettaPosZ);
        //dt
        G4double decayStart = G4UniformRand()*1000000*ns;
        // for electron
        G4double theta1 = acos(2.0 * G4UniformRand() - 1.0);
        G4double phi1 = 2.0 * G4UniformRand() * pi;

        Momentum1.setRThetaPhi(energy1, theta1, phi1); // Momentum1 == energy units in GEANT4
        G4ThreeVector dir1 = Momentum1.unit();

        // for positron
        G4double theta2 = acos(2.0 * G4UniformRand() - 1.0);
        G4double phi2 = 2.0 * G4UniformRand() * pi;

        Momentum2.setRThetaPhi(energy2, theta2, phi2); // Momentum1 == energy units in GEANT4
        G4ThreeVector dir2 = Momentum2.unit();

        //for decay electron
                dx = decayPos;// stays same for positron, as well. This is decay position of X virtual photon.
               dt = firstScatStart + decayStart;//dt can be higher than 1*ms
                G4PrimaryVertex* vertex1= new G4PrimaryVertex(dx, dt); ////be carefull with this (dx, dt) must be different from the initial value
                G4PrimaryParticle* geantino1 = new G4PrimaryParticle(particle1, Momentum1.x(),Momentum1.y(),Momentum1.z());// Must be electron

        // Generate random polarization
        phi1 = (G4UniformRand()*2.0-1.0)*M_PI;
        G4ThreeVector el1 = dir.orthogonal().unit();
        G4ThreeVector el2 = dir.cross(el1);
        G4ThreeVector rpol1 = el1*cos(phi1)+el2*sin(phi1);
        geantino1->SetPolarization(rpol1.x(), rpol1.y(), rpol1.z());

        // for positron

                G4PrimaryVertex* vertex2 = new G4PrimaryVertex(dx, dt); //be carefull with this (dx, dt) must be different from the initial value
                G4PrimaryParticle* geantino2 = new G4PrimaryParticle(particle2, Momentum2.x(),Momentum2.y(),Momentum2.z());// Must be positron

         // Generate random polarization
         phi2 = (G4UniformRand()*2.0-1.0)*M_PI;
         G4ThreeVector elp1 = dir.orthogonal().unit();
         G4ThreeVector elp2 = dir.cross(elp1);
        G4ThreeVector rpol2 = elp1*cos(phi2)+elp2*sin(phi2);
        geantino2->SetPolarization(rpol2.x(), rpol2.y(), rpol2.z());



        //    vertex->SetPrimary(geantino);
       //    event->AddPrimaryVertex(vertex);

        // Here should be something else such as an algorithm to define hit or miss method precisely for (e-,e+ pair production) decay of X (chi) DM intermediator resulting X1 and electron positron pair$
        // for now we have the below randomization code
        G4double probDecay = G4UniformRand(); //generates double values in (0-1) range

        if(probDecay < 0.5){

                vertex->SetPrimary(geantino);
                event->AddPrimaryVertex(vertex);
                G4cout<<"BDM: Elastic Scattering of X1 at: "<<Position<<G4endl;
        } else {
                vertex->SetPrimary(geantino);  //1) must be electron (from Ar ionization)
                vertex1->SetPrimary(geantino1); //2a)must be electron (from X intermediator decay)
                vertex2->SetPrimary(geantino2); //2b)must be positron (from X intermediator decay)

                event->AddPrimaryVertex(vertex);  // for 1)
                event->AddPrimaryVertex(vertex1); //for 2a)
                event->AddPrimaryVertex(vertex2); //for 2b)
                G4cout<<"iBDM: X1 scattering at: "<<Position<<G4endl;
                G4cout<<"iBDM: X1->(X2*->X(e-e+)X1. Virtual Photon Decay Position at: "<<dx<<G4endl;
        // try the code below with changes of
        /*
        G4PrimaryVertex* vertex =  new G4PrimaryVertex(G4ThreeVector(0.,0.,0.),GetTime());

        if(GetVerboseLevel() > 0)
               G4cout << "Creating primaries and assigning to vertex" << G4endl;

         for( G4int i=0; i<GetNumberOfParticles(); i++ )
                 {
                         G4PrimaryParticle* particle = new G4PrimaryParticle(m_pd,0. ,0. , 0.);
                         vertex->SetPrimary( particle );
                 }

                aEvent->AddPrimaryVertex( vertex );

        */
        }

//} else {
//      G4cout<<"You are out of the LAr volume. Find a way to enter inside of the LAr volume..."<<G4endl;
//}
}
  void VertexGen_BDS::SetState(G4String newValues)
  {
    if (newValues.length() == 0) {
      // print help and current state
      G4cout << "Current state of this VertexGen_BDS:\n"
             << " \"" << GetState() << "\"\n" << G4endl;
      G4cout << "Format of argument to VertexGen_BDS::SetState: \n"
        " \"ejectile mass in_MeV\"\n" << G4endl;
      return;
  }
    Setup(); //Should Configure energy distribution
  }

  G4String VertexGen_BDS::GetState()
  {
    return ("DEBUGED the GET Function.");
    //return dformat("%s\t%f", fNucleusName.c_str(), fmassGeantino);
  }
  void VertexGen_BDS::Setup()
  {
        G4cout<<"Sucessful Setup!"<<G4endl;
   // DBLinkPtr lwimp = DB::Get()->GetLink("WIMP");

    //const int nsamples = 10000;
    //double dRdQ[nsamples];

    //fEnergyDistLo = lwimp->GetD("energy_lo") * keV;
//  //fEnergyDistHi = lwimp->GetD("energy_hi") * keV;


    // Target properties
    //const double m_nucleus = fGeantino->GetPDGMass();
    //const double m_nucleus = fNucleus->GetPDGMass();

  }

}
