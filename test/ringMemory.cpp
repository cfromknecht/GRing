#include <ring/GRingArray.hpp>
#include <gtest/gtest.h>
#include <iostream>

TEST( memory_test, allocation_test ) {
  const size_t K = 14;
  const size_t Q = 1 << K;
  gring::GRingArray<256,Q,K,100>* rv1 = new gring::GRingArray<256,Q,K,100>;
  delete rv1;

  gring::GRingArray<256,Q,K,100>* rv2 = new gring::GRingArray<256,Q,K,100>[100];
  delete [] rv2;
}



