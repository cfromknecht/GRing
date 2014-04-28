#ifndef _IBE_H_
#define _IBE_H_

#include <ring/GRingArray.hpp>

namespace gring {

  template< size_t N, size_t Q, size_t K >
  class IBEMasterPublic {
  protected:
    GRingArray<N,Q,K,2> _mpk;
  public:
    IBEMasterPublic() : _mpk{GRingArray<N,Q,K,2>{}} {}
    IBEMasterPublic( GRingArray<N,Q,K,2>& existingMPK ) : _mpk{existingMPK} {}
    IBEMasterPublic( const IBEMasterPublic& other ) : IBEMasterPublic{other.mpk()} {}
    IBEMasterPublic& operator=( const IBEMasterPublic& other );
    virtual ~IBEMasterPublic() {}
    GRingArray<N,Q,K,2> mpk() const { return _mpk; }
  };

  template< size_t N, size_t Q, size_t K >
  class IBEMasterSecret : public IBEMasterPublic<N,Q,K> {
  private:
    GRingArray<N,Q,K,2> _msk;
  public:
    IBEMasterSecret();
    IBEMasterSecret( const IBEMasterSecret& other );
    IBEMasterSecret& operator=( const IBEMasterSecret& other );
    ~IBEMasterSecret() {}
    GRingArray<N,Q,K,2> msk() const { return _msk; }
    IBEMasterPublic<N,Q,K> master() const { return IBEMasterPublic<N,Q,K>{this->mpk()}; }
  };

  template< size_t N, size_t Q, size_t K >
  IBEMasterSecret<N,Q,K>::IBEMasterSecret() : IBEMasterPublic<N,Q,K>{}, 
      _msk{GRingArray<N,Q,K,2>{}} {
    this->mpk().trapGen( 0, _msk );
  }

  template< size_t N, size_t Q, size_t K >
  IBEMasterSecret<N,Q,K>::IBEMasterSecret( const IBEMasterSecret& other ) : 
      IBEMasterPublic<N,Q,K>{other.mpk()} {
    this->_msk = other.msk();
  }

  template< size_t N, size_t Q, size_t K >
  IBEMasterSecret<N,Q,K>&
  IBEMasterSecret<N,Q,K>::operator=( const IBEMasterSecret& other ) {
    this->_mpk = other.mpk();
    this->_msk = other.msk();
  }

} // namespace gring

#endif // _IBE_H_

