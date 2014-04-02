#ifndef _G_RING_ARRAY_
#define _G_RING_ARRAY_

#include <flint/fmpz_mod_poly.h>
#include <openssl/sha.h>
#include <omp.h>
#include <sys/time.h>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <cstring>

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

namespace gring {

  typedef fmpz_mod_poly_t ring_t;

  template< size_t N, size_t Q, size_t K, size_t LEN >
  class GRingArray {
  private:
    ring_t* const _polys;
    ring_t& _FF;
    ring_t& makeF() const; 
  public:
    GRingArray();
    GRingArray( const GRingArray& other );
    GRingArray& operator=( const GRingArray& other );
    ~GRingArray();
    // Accessors
    inline size_t NN() const { return N; }
    inline size_t QQ() const { return Q; }
    inline size_t KK() const { return K; }    
    inline size_t len() const { return LEN; }
    inline size_t size() const { return K*LEN; }
    inline ring_t* data() const { return _polys; }
    inline ring_t& FF() const { return _FF; }
    // Initializers
    void uniformInit( const size_t start, const size_t end );
    void ternaryInit( const size_t start, const size_t end );
    // Arithmetic Operations
    GRingArray operator*( const ring_t& s ) const;
    void operator*=( const ring_t& s );
    GRingArray operator*( const GRingArray& other ) const;
    void operator*=( const GRingArray& other );
    GRingArray operator+( const GRingArray& other ) const;
    void operator+=( const GRingArray& other );
    GRingArray operator-( const ring_t& other ) const;
    void operator-=( const ring_t& other );
    void plusGRingSection( const size_t start, const size_t end, const 
        GRingArray& other );
    void minusGRingSection( const size_t start, const size_t end, const 
        GRingArray& other );
    void mulmodGRingSection( const size_t start, const size_t end, const 
        GRingArray& other, GRingArray& target );
    // Cryptographic Operations
    void dualRegevEncode();
    void dualRegevEncode( ring_t& s );
    void dualRegevEncrypt( const unsigned char* const data );
    void dualRegevEncrypt( ring_t& s, const unsigned char* const data );
    void ternaryPerturb( ring_t& r );
    // Factory Methods
    void makeUniformPoly( ring_t& u );
    void makeParityCheckForSecret( GRingArray& AR, const size_t start, 
        const size_t end, const GRingArray& secret );
    void sampleD( const size_t s, const ring_t& u, const size_t start, const 
        size_t end, const ring_t p1, const ring_t p2, const GRingArray& Ta, 
        const GRingArray& sigma );
    void invert( ring_t s, ring_t* e, const ring_t* b );
  };

  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>::GRingArray() : _polys{new ring_t[K*LEN]}, _FF{makeF()} {
    assert( (1 << K) >= Q );
    assert( Q > (1 << (K-1)) );

    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
#pragma omp parallel for
    for ( size_t i = 0; i < K*LEN; ++i )
      fmpz_mod_poly_init2( _polys[i], q_z, N );
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>::GRingArray( const GRingArray& other ) : 
      _polys{new ring_t[K*LEN]}, _FF{makeF()} {
    assert( (1 << K) >= Q );
    assert( Q > (1 << (K-1)) );

    ring_t* otherPolys = other.data();
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
#pragma omp parallel for
    for ( size_t i = 0; i < K*LEN; ++i ) {
      fmpz_mod_poly_init2( _polys[i], q_z, N );
      fmpz_mod_poly_set( _polys[i], otherPolys[i] );
    }
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>::~GRingArray() {
#pragma omp parallel for
    for ( size_t i = 0; i < K*LEN; ++i )
      fmpz_mod_poly_clear( _polys[i] );
    delete [] _polys;
    fmpz_mod_poly_clear( _FF );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::uniformInit( const size_t start, const size_t end ) {
#pragma omp parallel for 
    for ( size_t i = K*start; i < K*end; ++i )
      for ( size_t j = 0; j < N; ++j )
        fmpz_mod_poly_set_coeff_ui( _polys[i], j, rand() );
    return;
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::ternaryInit( const size_t start, const size_t end ) {
#pragma omp parallel for
    for ( size_t i = K*start; i < K*end; ++i )
      for ( size_t j = 0; j < N; ++j )
        fmpz_mod_poly_set_coeff_ui( _polys[i], j, (rand()&1) + (rand()&1)*(Q-1) );
    return;
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>& 
  GRingArray<N,Q,K,LEN>::operator=( const GRingArray& other ) {
    ring_t* otherPolys = other.data();
#pragma omp parallel for
    for ( size_t i = 0; i < K*LEN; ++i )
      fmpz_mod_poly_set( _polys[i], otherPolys[i] );
    return *this;
  }


  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>
  GRingArray<N,Q,K,LEN>::operator*( const ring_t& s ) const {
    GRingArray rv{*this};
    rv *= s;
    return rv;
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void 
  GRingArray<N,Q,K,LEN>::operator*=( const ring_t& s ) {
    ring_t& F = FF();
#pragma omp parallel for
    for ( size_t i = 0; i < K*LEN; ++i )
      fmpz_mod_poly_mulmod( _polys[i], _polys[i], s, F );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>
  GRingArray<N,Q,K,LEN>::operator*( const GRingArray& other ) const {
    GRingArray rv{*this};
    rv *= other;
    return rv;
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void 
  GRingArray<N,Q,K,LEN>::operator*=( const GRingArray& other ) {
    mulmodGRingSection( 0, LEN, other, *this);
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>
  GRingArray<N,Q,K,LEN>::operator+( const GRingArray& other ) const {
    GRingArray rv{*this};
    rv += other;
    return rv;
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void 
  GRingArray<N,Q,K,LEN>::operator+=( const GRingArray& other ) {
    plusGRingSection( 0, LEN, other );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  GRingArray<N,Q,K,LEN>
  GRingArray<N,Q,K,LEN>::operator-( const ring_t& other ) const {
    GRingArray rv{*this};
    rv -= other;
    return rv;
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void 
  GRingArray<N,Q,K,LEN>::operator-=( const ring_t& other ) {
    minusGRingSection( 0, LEN, other );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::plusGRingSection( const size_t start, const size_t end, 
      const GRingArray& other ) {
    if ( end > other.len() ) throw std::runtime_error{ "end > other.len()" };
    ring_t* otherPolys = other.data();
    for ( size_t i = K*start; i <= K*end; i++ )
      fmpz_mod_poly_add( _polys[i], _polys[i], otherPolys[i] );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::minusGRingSection( const size_t start, const size_t end, 
      const GRingArray& other ) {
    if ( end > other.len() ) throw std::runtime_error{ "end > other.len()" };
    ring_t* otherPolys = other.data();
    for ( size_t i = K*start; i <= K*end; i++ )
      fmpz_mod_poly_sub( _polys[i], _polys[i], otherPolys[i] );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::mulmodGRingSection( size_t start, size_t end, const 
      GRingArray& other, GRingArray& target) {
    if ( end > other.len() ) throw std::runtime_error{ "end > other.len()" };
    ring_t& F = FF();
    ring_t* otherPolys = other.data();
    ring_t* targetPolys = target.data();
    for ( size_t i = K*start; i < K*end; i++ )
      fmpz_mod_poly_mulmod( targetPolys[i], _polys[i], otherPolys[i], F );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::dualRegevEncode() {
    ring_t s;
    makeUniformPoly( s );
    dualRegevEncode( s );
    fmpz_mod_poly_clear( s );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::dualRegevEncode( ring_t& s ) {
    ring_t& F = FF();
#pragma omp parallel for
    for ( size_t i = 0; i < K*LEN; ++i ) {
      fmpz_mod_poly_mulmod( _polys[i], _polys[i], s, F );
      ternaryPerturb( _polys[i] );
    }
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::dualRegevEncrypt( const unsigned char* const mu ) {
    size_t blocksize = N/8;
    ring_t& F = FF();
    ring_t s;
    makeUniformPoly( s );
#pragma omp parallel for
    for ( size_t i = 0; i < K*LEN; ++i ) {
      fmpz_mod_poly_mulmod( _polys[i], _polys[i], s, F );
      fmpz_t coeff, cj;
      fmpz_init( coeff ); fmpz_init( cj );
      size_t byte, bit;
      for ( size_t j = 0; j < blocksize; ++j ) {
        byte = mu[blocksize*i + j];
        for ( size_t k = 0; k < 8; ++k ) {
          bit = byte & 1;
          fmpz_set_ui(coeff, bit *(Q/2) );
          fmpz_mod_poly_get_coeff_fmpz( cj, _polys[i], j );
          fmpz_add( cj, cj, coeff );
          fmpz_set_ui( cj, (rand()&1) + (rand()&1)*(Q-1) );
          fmpz_add( cj, cj, coeff );
          fmpz_mod_poly_set_coeff_fmpz( _polys[i], coeff, j );
          byte >>= 1;
        }
      }
    }
  }

  /* Instead of creating a random ternary error polynomial and adding, we can 
   * avoid creating ephemeral polynomials by perturbing each coefficient 
   * individually.
   */

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void 
  GRingArray<N,Q,K,LEN>::ternaryPerturb( ring_t& ri ) {
    fmpz_t coeff, eij;
    fmpz_init( eij ); fmpz_init( coeff );

    for ( size_t j = 0; j < N; ++j ) {
      fmpz_set_ui( eij, (rand()&1) + (rand()&1)*(Q-1) );
      fmpz_mod_poly_get_coeff_fmpz( coeff, ri, j );
      fmpz_add( coeff, coeff, eij );
      fmpz_mod_poly_set_coeff_fmpz( ri, j, coeff );
    }

    fmpz_clear( coeff ); fmpz_clear( eij );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::makeUniformPoly( ring_t& r ) {
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
    fmpz_mod_poly_init2( r, q_z, N );
    for ( size_t i = 0; i < N; ++i )
      fmpz_mod_poly_set_coeff_ui( r, i, rand() );
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::makeParityCheckForSecret( GRingArray& AR, const size_t 
      start, const size_t end, const GRingArray& secret ) {
    if ( end > secret.len() ) throw std::runtime_error{ "end > secret.len()" };
    if ( end > AR.len() ) throw std::runtime_error{ "end > AR.len()" };

    uniformInit( start, end );
    mulmodGRingSection( start, end, secret, AR);

    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );

    ring_t coeffPoly;
    fmpz_mod_poly_init( coeffPoly, q_z );
    ring_t* ARPolys = AR.data();
    for ( size_t i = start; i < end; ++i ) {
      for ( size_t j = 0; j < K; ++j ) {
        fmpz_mod_poly_set_coeff_ui( coeffPoly, 0, 1 << j );
        fmpz_mod_poly_sub( ARPolys[K*i+j], coeffPoly, ARPolys[K*i+j] );
      }
    }

    fmpz_mod_poly_clear( coeffPoly );
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  GRingArray<N,Q,K,LEN>::sampleD( const size_t s, const ring_t& u, const size_t 
      start, const size_t end, const ring_t p1, const ring_t p2, const 
      GRingArray& Ta, const GRingArray& sigma ) {
    //GRingArray wBar = (*this) * (p1 - (Ta * p2));
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
  ring_t&
  GRingArray<N,Q,K,LEN>::makeF() const {
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
    ring_t* F = new ring_t[1];
    fmpz_mod_poly_init2( F[0], q_z, N );
    fmpz_mod_poly_set_coeff_ui( F[0], 0, 1 );
    fmpz_mod_poly_set_coeff_ui( F[0], N, 1 );
    return F[0];
  }

  template< size_t N, size_t Q, size_t K, size_t LEN >
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

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void
  invertY( GRingArray<N,Q,K,LEN>& rPrime, const ring_t y ) {
    size_t yj;
    fmpz_t y_z;
    fmpz_init( y_z );

    ring_t* rPrimePolys = rPrime.data();

    for ( size_t j = 0; j < N; ++j ) {
      fmpz_mod_poly_get_coeff_fmpz( y_z, y, j );
      yj = fmpz_get_ui( y_z );
      for ( size_t i = 0; i < K; ++i ) {
        fmpz_mod_poly_set_coeff_ui( rPrimePolys[i], j, yj & 1 );
        yj >>= 1;
      }
    }
    fmpz_clear( y_z );
  }

  /**                                                                                 
   * Computes the SHA512 hash of str and returns it as hexidecimal characters         
   */                                                                                 
                                                                                      
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
                                                                                      
  /**                                                                                 
   * Converts str into a polynomial of length n by first hashing and then setting  
   * the coefficients of y in ZZ_q. Accepts k in [8, 16]                              
   */                                                                                 

  template< size_t N, size_t Q, size_t K, size_t LEN >
  void                                                                             
  stringToModQ( fmpz_mod_poly_t y, const char* str ) {                                        
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

  /**                                                                                 
   * Creates y in R_q = H(id)                                                         
   */                                                                                 
                                                                                      
  template< size_t N, size_t Q, size_t K, size_t LEN >
  void                                                                                
  idToY( fmpz_mod_poly_t y, const char* id ) {                                          
      size_t num_hashes = std::ceil( N * K / 512.0 );                                      
      char* hashed = sha512Chain( id, num_hashes );                              
      stringToModQ<N,Q,K,LEN>( y, hashed );                                          
      delete [] hashed;
  }                                                                                   
                                                                                      
  template< size_t N, size_t Q, size_t K, size_t LEN >
  void                                                                                
  invertId( ring_t y, GRingArray<N,Q,K,LEN>& rPrime, const char* id ) {                                        
    idToY<N,Q,K,LEN>( y, id );                                                      
    invertY<N,Q,K,LEN>( rPrime, y );                                                
  }    
     

} // namespace gring

#endif

