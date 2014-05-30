#ifndef _TRAPDOOR_H_
#define _TRAPDOOR_H_

#include <ring/GRingArray.hpp>
#include <string>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>

namespace gring {
   
  template< size_t N, size_t Q, size_t K >
  void 
  trapGen( GRingArray<N,Q,K>* pub, GRingArray<N,Q,K>* secret, const size_t 
      start ) {
    assert( start + 1 < pub->len() );
    assert( start + 1 < secret->len() );

    secret->ternaryInit( start, start + 1 );
    secret->identityInit( start + 1 );

    pub->uniformInit( start, start + 1 );
    pub->mulmodGRingSection( start, start + 1, *secret, start, *pub, start + 1 );

    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );

    ring_t coeffPoly;
    fmpz_mod_poly_init2( coeffPoly, q_z, 1 );
    
    ring_t* ARPolys = pub->data();
    for ( size_t j = K*(start + 1); j < K*(start + 2); ++j ) {
      fmpz_mod_poly_set_coeff_ui( coeffPoly, 0, (1 << j % K) );
      fmpz_mod_poly_sub( ARPolys[j], coeffPoly, ARPolys[j] );
    }

    fmpz_mod_poly_clear( coeffPoly );
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K >
  ring_t*
  dotWithG( GRingArray<N,Q,K>* X, const size_t index ) {
    ring_t* xData = X->data() + index*K;
    fmpz_t q_z; fmpz_init_set_ui( q_z, Q );
    fmpz_t coeff; fmpz_init( coeff );

    ring_t* u = new ring_t[1];
    ring_t temp;

    fmpz_mod_poly_init2( *u, q_z, N );
    fmpz_mod_poly_init2( temp, q_z, N );
    for ( size_t i = 0; i < K; ++i ) {
      fmpz_set_ui( coeff, 1 << i );
      fmpz_mod_poly_scalar_mul_fmpz( temp, xData[i], coeff );
      fmpz_mod_poly_add( *u, *u, temp );
    }
    fmpz_mod_poly_clear( temp );
    fmpz_clear( coeff );
    fmpz_clear( q_z );
    
    return u;
  }

  inline float
  uniformRand() {
    return rand() / float(RAND_MAX);
  }

  size_t
  sampleDiscreteZeroNormWithGaussParam( const long Q, const float param ) {
    const float c = 1.32;
    float denom = 2 * param * param;
    std::random_device rd;
    std::mt19937 generator{ rd() };
    std::uniform_int_distribution<long> distribution{0, Q};
    long Y;
    float U, threshold;
    do {
      Y = distribution( generator ) - Q/2;
      U = uniformRand();
      threshold = exp( -Y*Y/denom ) / c;
    } while ( U > threshold );
    return Y;
  }

  template< size_t N, size_t Q, size_t K >
  long*
  buildSamplesWithGuassianParameter( const size_t num, const float param ) {
    long* samples = new long[2*num];
    size_t evenIndex = 0, oddIndex = num;

    size_t t;
    while ( evenIndex < num || oddIndex < 2*num ) {
      t = sampleDiscreteZeroNormWithGaussParam( long(Q), param );
      if ( t & 1 ) { 
        if ( oddIndex < 2*num ) samples[oddIndex++] = t;
      } else if ( evenIndex < N*K ) samples[evenIndex++] = t;
    }
    return samples;
  }

  template< size_t N, size_t Q, size_t K >
  void
  sampleZZ( GRingArray<N,Q,K>* X, long* samples ) {
    std::random_device rd;
    std::mt19937 gen{ rd() };
    std::uniform_int_distribution<int> dist{0, 2*N*K};

    ring_t* xData = X->data();

    long coeff;
    for ( size_t i = 0; i < X->size(); ++i )
      for ( size_t j = 0; j < N; ++j ) {
        coeff = samples[dist( gen )];
        fmpz_mod_poly_set_coeff_ui( xData[i], j, coeff >= 0 ? coeff : coeff + long(Q) );
      }
  }

  template< size_t N, size_t Q, size_t K >
  void
  sampleZZWithSupport( GRingArray<N,Q,K>* X, const GRingArray<N,Q,K>* T, long* samples ) {
    assert( X->len() == T->len() );

    sampleZZ( X, samples );
    *X *= *T;
  }

  template< size_t N, size_t Q, size_t K >
  void
  solveRWeighted( GRingArray<N,Q,K>* R, std::vector< size_t > weights ) {
    assert( R->len()/(4*K) == weights.size() );

    ring_t* rData = R->data();
    long coeff;
    for ( size_t i = 0; i < weights.size(); ++i )
      for ( size_t k = 0; k < 2*K; ++k ) {
        coeff = weights[i] * (size_t(1) << k);
        for ( size_t j = 0; j < 2*K; ++j ) {
          fmpz_mod_poly_set_coeff_ui( rData[4*K*K*i + 2*K*k + j], 0, coeff & 1 );
          coeff /= 2;
        }
      }
  }

  template< size_t N, size_t Q, size_t K >
  void
  solveRWeightedSection( GRingArray<N,Q,K>* R, size_t weight ) {
    ring_t* rData = R->data();
    long coeff;
    for ( size_t i = 0; i < 2*K; ++i ) {
      coeff = weight * (size_t(1) << i );;
      for ( size_t j = 0; j < 2*K; ++j ) {
        fmpz_mod_poly_set_coeff_ui( rData[2*K*i + j], 0, coeff & 1 );
        coeff /= 2;
      }
    }
  }
  
  template< size_t N, size_t Q, size_t K >
  void
  solveRTarget( GRingArray<N,Q,K>* R, GRingArray<N,Q,K>* D, const size_t section ) {
    assert( 4*K*section + 2*K*D->len() <= R->len() );


    fmpz_t coeff_z; fmpz_init( coeff_z );
    long coeff;
    
    ring_t* rData = R->data();
    ring_t* dData = D->data();
    for ( size_t i = 0; i < D->size(); ++i )
      for ( size_t j = 0; j < N; ++j ) {
        fmpz_mod_poly_get_coeff_fmpz( coeff_z, dData[i], j );
        coeff = fmpz_get_ui( coeff_z );
        for ( size_t k = 0; k < 2*K; ++k ) {
          fmpz_mod_poly_set_coeff_ui( rData[4*K*K*section + 2*K*i + k], j, coeff & 1 );
          coeff /= 2;
        }
      }
    fmpz_clear( coeff_z );
  }

  template< size_t N, size_t Q, size_t K >
  void
  sampleUFromG( ring_t* xData, ring_t& u, long* samples ) {
    long* evenSamples = samples;
    long* oddSamples = samples + N*K;

    fmpz_t coeff_z; fmpz_init( coeff_z );
    long coeff, xi;
    
    size_t evenIndex = 0, oddIndex = 0;
    for ( size_t i = 0; i < N; ++i ) {
      fmpz_mod_poly_get_coeff_fmpz( coeff_z, u, i );
      coeff = fmpz_get_ui( coeff_z );
      for ( size_t j = 0; j < K; ++j ) {
        if ( coeff & 1 ) { xi = oddSamples[oddIndex++]; oddIndex %= N*K; } 
        else { xi = evenSamples[evenIndex++]; evenIndex %= N*K; }
        fmpz_mod_poly_set_coeff_ui( xData[j], i, xi >= 0 ? xi : xi + long(Q) );
        coeff = (coeff - xi)/2;
      }
    }
    fmpz_clear( coeff_z );
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  sampleD( GRingArray<N,Q,K>* R, GRingArray<N,Q,K>* A, GRingArray<N,Q,K>* U, long* samples ) {
    assert( A->len() == R->len() );

    size_t rLen = R->len();
    size_t rSize = R->size();
    size_t uLen = U->len();
    size_t uSize = U->size();

    auto p = GRingArray<N,Q,K>{rLen};
    auto Ap = GRingArray<N,Q,K>{rLen};
    auto v = GRingArray<N,Q,K>{uLen};
    auto z = GRingArray<N,Q,K>{uSize};
    auto X = new GRingArray<N,Q,K>{rLen*uSize};
    
    // draw fresh p with support for R 
    sampleZZWithSupport<N,Q,K>( &p, R, samples );
    // calculate Ap
    A->mulmodGRingSection( 0, rLen, p, 0, Ap, 0 );
    for ( size_t i = 1; i < rSize; ++i )
      fmpz_mod_poly_add( Ap.data()[0], Ap.data()[0], Ap.data()[i] );
    // v = U - Ap
    U->minusGRingSection( 0, rLen, Ap.data()[0], v, 0 );
    // sample z <- Dv,r*sqrt(Eg) 
//#pragma omp parallel for
    for ( size_t i = 0; i < uSize; ++i )
      sampleUFromG<N,Q,K>( z.data() + K*i, v.data()[i], samples );

    /*
    std::cout << "z: " << std::endl;
    fmpz_mod_poly_print( z.data()[0] );
    std::cout << std::endl;
    */
    
    // X = p + R * z
//#pragma omp parallel for
    for ( size_t i = 0; i < uSize; ++i )
      for ( size_t j = 0; j < rLen; ++j )
        R->mulmodGRingSection( j, j+1, z, i, *X, i*rLen + j );
//#pragma omp parallel for
    for ( size_t i = 0; i < rLen*uSize; i += rLen )
      p.plusGRingSection( 0, rLen, *X, i, *X, i );

    return X;
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  extendRightSampleD( GRingArray<N,Q,K>* R, GRingArray<N,Q,K>* A, 
      GRingArray<N,Q,K>* Abar, GRingArray<N,Q,K>* U, long* samples ) {
    assert( R->len() == U->len() );
    assert( R->len() == Abar->len() );

    auto vBar = GRingArray<N,Q,K>{4*K};
    auto AbarVBarTemp = GRingArray<N,Q,K>{4*K};
    auto yBar = GRingArray<N,Q,K>{2};
    auto v = new GRingArray<N,Q,K>{8*K};

//    std::cout << "sampling ZZ ..." << std::endl;
    // sample vBar from integers 
    sampleZZ<N,Q,K>( &vBar, samples );

//    std::cout << "mulmod ..." << std::endl;
    // calculate yBar = Abar * vBar 
    for ( size_t i = 0; i < 4*K; i += 2 )
      Abar->mulmodGRingSection( 0, 2, vBar, i, AbarVBarTemp, i );
//    std::cout << "sum ..." << std::endl;
#pragma omp parallel for
    for ( size_t i = 0; i < 2*K; ++i )
      for ( size_t j = 0; j < 2*K; ++j )
        fmpz_mod_poly_add( yBar.data()[i], yBar.data()[i], 
            AbarVBarTemp.data()[2*K*i+j] );
//    std::cout << "subtract ..." << std::endl;
    // calculate y - yBar 
    U->minusGRingSection( 0, 2, yBar, 0, yBar, 0 );
    // sample v from y - yBar
//    std::cout << "sampleD ..." << std::endl;
    auto vTemp = sampleD( R, A, &yBar, samples );

    /*
    std::cout << "vTemp: " << std::endl;
    fmpz_mod_poly_print( vTemp->data()[0] );
    std::cout << std::endl;
    */

//    std::cout << "copy ..." << std::endl;
    // copy into v
    v->identityInit( 0, v->len() );
    v->mulmodGRingSection( 0, v->len()/2, *vTemp, 0, *v, 0 );
    v->mulmodGRingSection( v->len()/2, v->len(), vBar, 0, *v, v->len()/2 );

    return v;
  }

  template< size_t N, size_t Q, size_t K >
  size_t
  writeSamples( std::string filename, long* samples, const long num, const long param) {
    std::stringstream ss;
    ss << "_" << param << ".dat";
    std::ofstream fout{ filename + ss.str() };
    if ( !fout.is_open() ) return 0;

    auto start = fout.tellp();

    fout << num << " " << param << "\n";
    for ( size_t i = 0; i < size_t(num) - 1; ++i )
      fout << samples[i] << " ";
    fout << samples[num-1] << "\n";
    for ( size_t i = num; i < 2*size_t(num)-1; ++i )
      fout << samples[i] << " ";
    fout << samples[2*num-1] << "\n";
    fout << "\0";

    auto end = fout.tellp();
    fout.close();
    return end - start;
  }

  template< size_t N, size_t Q, size_t K >
  long*
  readSamples( std::string filename, const long param ) {
    std::stringstream extn;
    extn << "_" << param << ".dat";
    std::ifstream fin{ filename + extn.str() };
    if ( !fin.is_open() ) return nullptr;

    long n, p;
    fin >> n >> p;
    if ( n < long(N*K) ||  p != param ) return nullptr;

    long* samples = new long[2*N*K];
    for ( size_t i = 0; i < N*K; ++i )
      fin >> samples[i];
    for ( size_t i = N*K; i < size_t(n); ++i )
      fin >> samples[N*K];
    for ( size_t i = N*K; i < 2*N*K; ++i )
      fin >> samples[i];

//    std::cout << "samples " << samples[0] << " " << samples[N*K] << std::endl;

    fin.close();
    return samples;
  }

  template< size_t N, size_t Q, size_t K >
  long*
  readOrBuildSamples( std::string filename, const long param ) {
    long* samples = readSamples<N,Q,K>( filename, param );
    if ( samples ) return samples;
    samples = buildSamplesWithGuassianParameter<N,Q,K>( N*K, param );
    writeSamples<N,Q,K>( filename, samples, N*K, param );
    return samples;
  }


  template< size_t N, size_t Q, size_t K >
  void
  invert( ring_t s, ring_t* e, const ring_t* b ) {
    fmpz_t bj_z;
    fmpz_init( bj_z );

    long bj, r;
    size_t sj, notClose;
    for ( size_t j  = 0; j < N; ++j ) {
      sj = 0;
      for ( size_t i = 0; i < K; ++i ) {
        fmpz_mod_poly_get_coeff_fmpz( bj_z, b[i], j );
        bj = fmpz_get_ui( bj_z );
        r = (bj - (1 << i) * sj) % Q;
        notClose = (r >> (K - 1)) ^ ((r >> (K - 2)) & 1 );
        sj += (1 << (K - i - 1)) * notClose;
        fmpz_mod_poly_set_coeff_ui( e[i], j, bj - (1 << i) * sj );
      }
      fmpz_mod_poly_set_coeff_ui( s, j, sj );
    }
    fmpz_clear( bj_z );
  }

  template< size_t N, size_t Q, size_t K >
  void
  invertY( ring_t y, GRingArray<N,Q,K>* target, const size_t start ) {
    assert( (start + 1)*K <= K*target->len() );
    size_t yj;
    fmpz_t y_z;
    fmpz_init( y_z );

    for ( size_t j = 0; j < N; ++j ) {
      fmpz_mod_poly_get_coeff_fmpz( y_z, y, j );
      yj = fmpz_get_ui( y_z );
      for ( size_t i = start*K; i < (start + 1)*K; ++i ) {
        fmpz_mod_poly_set_coeff_ui( target->data()[i], j, yj & 1 );
        yj >>= 1;
      }
    }
    fmpz_clear( y_z );
  }

   // Computes the SHA512 hash of str and returns it as hexidecimal characters         
                                                                                      
  char*                                                                              
  sha512Chain( const char* str, const size_t num_hashes ) {                                       
    size_t hash_len = SHA512_DIGEST_LENGTH * num_hashes + 1;                        
    char* hashed = new char[hash_len];
    unsigned char digest[SHA512_DIGEST_LENGTH];                                     
    for ( size_t j = 0; j < num_hashes; j++ ) {                                     
        if ( j == 0 )                                                               
            SHA512( (unsigned char*) str, strlen( str ), digest );                  
        else                                                                        
            SHA512( (unsigned char*) hashed, j * SHA512_DIGEST_LENGTH, digest );  
        for ( size_t i = 0; i < SHA512_DIGEST_LENGTH; i++ )                         
            sprintf( &hashed[j * SHA512_DIGEST_LENGTH + i], "%0x", (unsigned char) digest[i] );
    }                                                                               
    hashed[hash_len] = '\0';                                                        
    return hashed;                                                                  
  }                                                                                   
                                                                                      
   // Converts str into a polynomial of length n by first hashing and then setting  
   // the coefficients of y in ZZ_q. Accepts k in [8, 16]                              

  void stringToModQ( fmpz_mod_poly_t y, const char* str, const size_t N, 
      const size_t K ) {                                        
    size_t i = 0, offset = 0;                                                    
    size_t mask, val, tempint;                                                   
                                                                                 
    fmpz_t temp;                                                                 
    fmpz_init( temp );                                                           
                                                                                 
    for ( size_t polyi = 0; polyi < N; ) {                                       
        if ( offset == 0 ) {                                                     
          if ( K == 16 ) {
            fmpz_mod_poly_set_coeff_ui( y, polyi, str[i] + (str[i+1] << 8) );
            i += 2;                                                          
          }                                                                    
          else                                                                 
            fmpz_mod_poly_set_coeff_ui( y, polyi, str[i++] );                
          if ( K % 8 == 0 )                                                    
            polyi++;                                                         
          offset = (offset + K) % 8;                                           
        } else {                                                                 
          mask = ((1 << offset) - 1);                                          
          val = (mask & str[i]) << (K - offset);                               
          fmpz_mod_poly_get_coeff_fmpz( temp, y, polyi );                      
          tempint = fmpz_get_ui( temp );                                       
          fmpz_mod_poly_set_coeff_ui( y, polyi, tempint + val );               
          polyi++;                                                             
          if ( polyi < N ) {                                                   
              mask = (1 << (8 - offset)) << offset;                            
              val = (mask & str[i]) >> offset;                                 
              fmpz_mod_poly_set_coeff_ui( y, polyi, val );                     
              i++;                                                             
              offset = (offset + K) % 8;                                       
          }                                                                    
        }                                                                        
    }                                                                            
                                                                                 
    fmpz_clear( temp );                                                          
  } 

   // Creates y in R_q = H(id)                                                         
                                                                                      
  void idToY( ring_t y, const std::string id, const size_t N, const size_t K ) {                                          
    size_t num_hashes = std::ceil( N * K / 512.0 );                                      
    char* hashed = sha512Chain( id.c_str(), num_hashes );                              
    stringToModQ( y, hashed, N, K );                                          
    delete [] hashed;
  }                                                                                   
                                                                                      
  template< size_t N, size_t Q, size_t K >
  void invertID( ring_t y, const std::string id, GRingArray<N,Q,K>* target, 
      const size_t start ) {                                        
    idToY( y, id, N, K );                                                      
    invertY( y, target, start );                                                
  }    

  template< size_t N, size_t Q, size_t K >
  void
  encode( ring_t& mu, std::string message ) {
    size_t numbytes = N/8;
    size_t slen = message.size();
    size_t len = slen < numbytes ? slen : numbytes;
    size_t si;
    size_t halfq = Q/2;
    size_t randbuff = 0;
    for ( size_t i = 0; i < len; ++i ) {
      if ( randbuff % 4 == 0 ) randbuff = rand();
      si = (size_t) message.c_str()[i];
      for ( size_t j = 0; j < 8; ++j ) {
        fmpz_mod_poly_set_coeff_ui( mu, 8*i + j, halfq * (si & 1) + (randbuff&1) + ((randbuff>>1)&1)*(Q-1));
        si >>= 1;
        randbuff >>= 2;
      }
    }
  }

  template< size_t N, size_t Q, size_t K >
  std::string
  decode( ring_t d ) {
    size_t quarter = Q/4;
    size_t numbytes = N/8;
    size_t msgi, di;

    std::string retString;

    fmpz_t di_z;
    fmpz_init( di_z );
    for ( size_t i = 0; i < numbytes; ++i ) {
      msgi = 0;
      for ( size_t j = 0; j < 8; ++j ) {
        fmpz_mod_poly_get_coeff_fmpz( di_z, d, 8*i + j );
        di = fmpz_get_ui( di_z );
        if ( di > quarter && di <= Q - quarter )
          msgi += (1 << j);
      }
      retString += (char) msgi;
    }
    fmpz_clear( di_z );
    
    return retString;
  }

} // namespace gring

#endif // _TRAPDOOR_H_

