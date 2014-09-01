// $Id: ensem_random.cc,v 2.0 2008/12/05 04:43:34 edwards Exp $
//
// Random number generator support


#include "ensem/ensem.h"

namespace ENSEM {

// Random number generator namespace
/* 
 * A collection of routines and data for supporting random numbers
 * 
 * It is a linear congruential with modulus m = 2**47, increment c = 0,
 * and multiplier a = (2**36)*m3 + (2**24)*m2 + (2**12)*m1 + m0.  
 */

namespace RNG
{
  //! Global (current) seed
  Seed ran_seed;

  //! RNG multiplier
  Seed ran_mult;

  //! Hook to get code to link
  const bool registered = initDefaultRNG();

  
  //! Find the number of bits required to represent x.
  int numbits(int x)
  {
    int num = 1;
    int iceiling = 2;
    while (iceiling <= x)
    {
      num++;
      iceiling *= 2;
    }

    return num;
  }


  //! Initialize the random number generator with a default seed
  bool initDefaultRNG()
  {
    RNG::initRNG();

    Seed seed = 11;
    RNG::setrn(seed);
    return true;
  }


  //! Initialize the internals of the random number generator
  void initRNG()
  {
    /* Multiplier used. Use big integer arithmetic */
    Seed seed_tmp3;
    Seed seed_tmp2;
    Seed seed_tmp1;
    Seed seed_tmp0;

    seed_tmp3 = 1222;
    seed_tmp2 = (seed_tmp3 << Int(12)) | Int(1498);
    seed_tmp1 = (seed_tmp2 << Int(12)) | Int(712);
    seed_tmp0 = (seed_tmp1 << Int(12)) | Int(1645);

    ran_mult = seed_tmp0;
  }


  //! Initialize the random number generator seed
  void setrn(const Seed& seed)
  {
    ran_seed = seed;
  }


  //! Return a copy of the random number seed
  void savern(Seed& seed)
  {
    seed = ran_seed;
  }


  //! Scalar random number generator. Done on the front end. */
  /*! 
   * It is linear congruential with modulus m = 2**47, increment c = 0,
   * and multiplier a = (2**36)*m3 + (2**24)*m2 + (2**12)*m1 + m0.  
   */
  float sranf()
  {
    bool reg = RNG::registered;

    /* Calculate the random number and update the seed according to the
     * following algorithm
     *
     * FILL(twom11,TWOM11);
     * FILL(twom12,TWOM12);
     * i3 = ran_seed(3)*ran_mult(0) + ran_seed(2)*ran_mult(1)
     *    + ran_seed(1)*ran_mult(2) + ran_seed(0)*ran_mult(3);
     * i2 = ran_seed(2)*ran_mult(0) + ran_seed(1)*ran_mult(1)
     *    + ran_seed(0)*ran_mult(2);
     * i1 = ran_seed(1)*ran_mult(0) + ran_seed(0)*ran_mult(1);
     * i0 = ran_seed(0)*ran_mult(0);
     *
     * ran_seed(0) = mod(i0, 4096);
     * i1          = i1 + i0/4096;
     * ran_seed(1) = mod(i1, 4096);
     * i2          = i2 + i1/4096;
     * ran_seed(2) = mod(i2, 4096);
     * ran_seed(3) = mod(i3 + i2/4096, 2048);
     *
     * sranf = twom11*(TO_REAL32(VALUE(ran_seed(3)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(2)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(1)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(0)))))));
     */

    Seed ran_tmp = ran_seed * ran_mult;
    ran_seed = ran_tmp;

    Real _sranf = seedToFloat(ran_seed);
    float _ssranf = toFloat(_sranf);

    return _ssranf;
  }

}

} // namespace ENSEM
