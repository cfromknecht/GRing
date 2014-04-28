#ifndef _TRAPDOOR_H_
#define _TRAPDOOR_H_

#include <string>

namespace gring {

    template< size_t N, size_t Q, size_t K >
    void 
    trapGen( GRingArray<N,Q,K>* pub, GRingArray<N,Q,K>* secret, const size_t start ) {
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
      fmpz_mod_poly_set_coeff_ui( coeffPoly, 0, 1 << (j % K) );
      fmpz_mod_poly_sub( ARPolys[j], coeffPoly, ARPolys[j] );
    }

    fmpz_mod_poly_clear( coeffPoly );
    fmpz_clear( q_z );
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

  /**                                                                                 
   * Creates y in R_q = H(id)                                                         
   */                                                                                 
                                                                                      
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
    size_t len = (numbytes + ((slen - numbytes) & -(slen < numbytes)));
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

