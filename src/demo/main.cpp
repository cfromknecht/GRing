#include <ring/ibe.hpp>
#include <ring/abe.hpp>

#include <cstdlib>
#include <sys/time.h>

const size_t N = 256;          // main security parameter
const size_t K = 24;           // lg( Q )
const size_t L = 5;            // Length of ABE id
const long SS = 4;             // smoothing parameter

std::default_random_engine generator;

std::string msg = "Lorem ipsum dolor sit amet, verterem delicata qui an, "
                  "iudico utinam eum ne. Ut vide zril corrumpit qui, ius "
                  "dolorem pertinax ex. Quo modus conceptam cu. Vivendo "
                  "conclusionemque cu has. Sanctus dolorem dissentiet ut vim, "
                  "quidam eloquentiam ut eos.  Eam everti scripta dissentiet "
                  "ad, ignota libris accusata his at. Legere populo elaboraret "
                  "ut sea, in nibh vivendo splendide pro. Pri te wisi ferri "
                  "ullamcorper, novum vitae feugait cu his. Per ut dolor "
                  "graece, at nec eripuit interesset. Malis velit quo ei, "
                  "utinam eripuit ne nec.  Id justo volutpat mei, eos ei veri "
                  "dolores invidunt, per at graece putent causae. Te sea "
                  "aperiam eleifend sententiae. Ne laudem bonorum volutpat "
                  "sit. Qui eu erat aperiam. Ius malis graece animal an, "
                  "numquam scribentur eam et. Imperdiet efficiantur "
                  "definitionem mea ei, vix enim vide evertitur in.  No sed "
                  "movet iracundia, vel veritus lucilius honestatis id, enim "
                  "nusquam an eos. Graeci contentiones at eos, id aliquid "
                  "noluisse iracundia nec. Te pertinax elaboraret per. Duo eu "
                  "alia insolens repudiandae.  Fabulas quaerendum sit eu, ne "
                  "mentitum postulant consequuntur sit, illum mazim vix an. Ea "
                  "prima accusamus vim. Diam mollis volutpat eu usu, sed ut "
                  "illum verear prodesset, usu no tation diceret ancillae. Mel "
                  "aeque dolore ei. Ius interesset complectitur an, pro vitae "
                  "antiopam ut. Quem impedit vim eu.";


void 
printTimes( std::string action, const struct timeval t1, const struct 
    timeval t2 ) {
  long millis = t2.tv_usec - t1.tv_usec + 1000000.0*(t2.tv_sec - t1.tv_sec);
  std::cout << action << " time: " << millis << std::endl;
  std::cout << action << "s/sec: " << 1000000.0/millis << std::endl;
}

void IBEDemo() {
  std::cout << "IBE" << std::endl;

  const size_t Q = size_t(1) << K;

  struct timeval t1, t2;
  std::string idString = "conner";

  // Make MSK and PK
  gettimeofday( &t1, 0 );
  gring::IBEMasterSecret<N,Q,K>* MSK = new gring::IBEMasterSecret<N,Q,K>{};
  gettimeofday( &t2, 0 );
  printTimes( "setup", t1, t2 );

  gring::IBEPublic<N,Q,K>* PK = MSK->publicKey();

  // Encrypt to ID
  gettimeofday( &t1, 0 );
  gring::IBECiphertext<N,Q,K>* ctxt = PK->encrypt( idString, "something much "
      "longer than 32 bytes so that maybe it gets cut off?" );
  gettimeofday( &t2, 0 );
  printTimes( "encryption", t1, t2 );

  // Invert Secret Key
  gettimeofday( &t1, 0 );
  gring::IBESecret<N,Q,K>* userSecret = MSK->secretKeyForID( idString );
  gettimeofday( &t2, 0 );
  printTimes( "inversion", t1, t2 );

  // Decrypt
  gettimeofday( &t1, 0 );
  std::string messagePrime = userSecret->decrypt( ctxt );
  gettimeofday( &t2, 0 );
  printTimes( "decryption", t1, t2 );

  std::cout << "messagePrime: " << messagePrime << std::endl;
}

void ABEDemo() {
  std::cout << "ABE" << std::endl;

  const size_t Q = size_t(1) << K;

  std::cout << "generating samples ..." << std::endl;
  long* samples = gring::readOrBuildSamples<N,Q,K>( "samples", SS );

  std::cout << "creating master secret ..." << std::endl;
  auto msk = new gring::ABEMasterSecret<N,Q,K>{L, samples};
  std::cout << "creating public..." << std::endl;
  auto pk = msk->publicKey();

  std::vector< size_t > addWeights{1, 2, 3, 4, 5};
  auto gAdd = new gring::ABEAddGate<N,Q,K>{addWeights};
  gAdd->print();

  std::vector< size_t > mulWeights{1}; 
  auto gMul = new gring::ABEMulGate<N,Q,K>{mulWeights};
  gMul->print();

  std::vector< size_t > badWeights{5, 4, 3, 2, 1};
  auto gBad = new gring::ABEAddGate<N,Q,K>{badWeights};
  gBad->print();

  std::cout << "creating secret keys ..." << std::endl;
  auto skAdd = msk->keyGen( gAdd );
  auto skMul = msk->keyGen( gMul );
  auto skBad = msk->keyGen( gBad );

  std::vector< size_t > ident{0, Q-1, Q-1, 0, 1};
  auto abeId = gring::ABEId{ident};
  abeId.print();

  std::cout << "encrypting ..." << std::endl;
  auto ctxt = pk->encrypt( abeId, msg );

  std::cout << "decrypting with skAdd ..." << std::endl;
  std::string msgPrimeAdd = skAdd->decrypt( ctxt );
  std::cout << "msgPrimeAdd: " << msgPrimeAdd << std::endl;
  std::cout << std::endl;

  std::cout << "decrypting with skMul ..." << std::endl;
  std::string msgPrimeMul = skMul->decrypt( ctxt );
  std::cout << "msgPrimeMul: " << msgPrimeMul << std::endl;
  std::cout << std::endl;

  std::cout << "decrypting with skBad ..." << std::endl;
  std::string msgPrimeBad = skBad->decrypt( ctxt );
  std::cout << "msgPrimeBad: " << msgPrimeBad << std::endl;

  delete [] samples;
}

int main( void ) {
  srand( time( NULL ) );

//  IBEDemo();
  ABEDemo();
  
  return 0;
}

