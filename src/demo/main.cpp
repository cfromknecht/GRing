#include <ring/ibe.hpp>
#include <ring/abe.hpp>

#include <cstdlib>
#include <sys/time.h>

const size_t N = 256;          // main security parameter
const size_t K = 30;           // lg( Q )
const size_t L = 5;            // Length of ABE id
const long SS = 4;             // smoothing parameter

std::default_random_engine generator;

void
printVector( std::vector< size_t > vec ) {
  for ( size_t i = 0; i < vec.size(); ++i )
    std::cout << vec[i] << " ";
  std::cout << std::endl;
}

void
buildCircuit() {
  /*
  std::vector< gring::ABEGate<N,Q,K>* > gates;
  gates.push_back( f );

  std::vector< size_t > fWire;
  fWire.push_back( 0 );
  fWire.push_back( 1 );
  fWire.push_back( 2 );
  fWire.push_back( 3 );
  fWire.push_back( 4 );
  
  std::vector< std::vector< size_t > > wires;
  wires.push_back( fWire );

  gring::ABECircuit<N,Q,K> circ{ gates, wires, L };
  std::cout << "evaluating circuit ..." << std::endl;
  std::cout << "eval: " << circ.evalId( abeId ) << std::endl;

  */
}

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
  auto abeMSK = new gring::ABEMasterSecret<N,Q,K>{L, samples};
  std::cout << "creating public..." << std::endl;
  auto abePK = abeMSK->publicKey();


  std::vector< size_t > weights{1}; // mul or 1 add
//  std::vector< size_t > weights{1, 2, 3, 4, 5}; //5 add

  auto f = new gring::ABEMulGate<N,Q,K>{weights};
  f->print();


  std::cout << "creating skF ..." << std::endl;
  auto abeSK = abeMSK->keyGen( f );

//  std::vector< size_t > ident{0};              // 1 mul or 1 add
  std::vector< size_t > ident{2, 3, 5, 7, Q/2};  // 5 mul
//  std::vector< size_t > ident{0, Q-1, Q-1, 0, 1};  // 5 add

  auto abeId = gring::ABEId{ident};
  std::cout << "ABEId: ";
  printVector( ident );

  std::string message = "Lorem ipsum dolor sit amet, verterem delicata qui an, "
                        "iudico utinam eum ne. Ut vide zril corrumpit qui, ius "
                        "dolorem pertinax ex. Quo modus conceptam cu. Vivendo "
                        "conclusionemque cu has. Sanctus dolorem dissentiet ut "
                        "vim, quidam eloquentiam ut eos.  Eam everti scripta "
                        "dissentiet ad, ignota libris accusata his at. Legere "
                        "populo elaboraret ut sea, in nibh vivendo splendide "
                        "pro. Pri te wisi ferri ullamcorper, novum vitae "
                        "feugait cu his. Per ut dolor graece, at nec eripuit "
                        "interesset. Malis velit quo ei, utinam eripuit ne "
                        "nec.  Id justo volutpat mei, eos ei veri dolores "
                        "invidunt, per at graece putent causae. Te sea aperiam "
                        "eleifend sententiae. Ne laudem bonorum volutpat sit. "
                        "Qui eu erat aperiam. Ius malis graece animal an, "
                        "numquam scribentur eam et. Imperdiet efficiantur "
                        "definitionem mea ei, vix enim vide evertitur in.  No "
                        "sed movet iracundia, vel veritus lucilius honestatis "
                        "id, enim nusquam an eos. Graeci contentiones at eos, "
                        "id aliquid noluisse iracundia nec. Te pertinax "
                        "elaboraret per. Duo eu alia insolens repudiandae.  "
                        "Fabulas quaerendum sit eu, ne mentitum postulant "
                        "consequuntur sit, illum mazim vix an. Ea prima "
                        "accusamus vim. Diam mollis volutpat eu usu, sed ut "
                        "illum verear prodesset, usu no tation diceret "
                        "ancillae. Mel aeque dolore ei. Ius interesset "
                        "complectitur an, pro vitae antiopam ut. Quem impedit "
                        "vim eu.";

  std::cout << "encrypting ..." << std::endl;
  auto abeCtxt = abePK->encrypt( abePK, abeId, message, abeSK->skF() );

  std::cout << "decrypting ..." << std::endl;
  std::string abeMsgPrime = abeSK->decrypt( abeCtxt );
  std::cout << "abeMsgPrime: " << abeMsgPrime << std::endl;

  delete [] samples;
}

int main( void ) {
  srand( time( NULL ) );

//  IBEDemo();
  ABEDemo();
  
  return 0;
}

