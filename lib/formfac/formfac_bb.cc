// $Id: formfac_bb.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
//
// Read building-blocks


#include "ensem/ensem.h"
#include "formfac/formfac_bb.h"

#include <sys/time.h>   // for timings
#include <assert.h>
#include <list>

#include "io/adat_byteorder.h"
#include "io/adat_io.h"

using namespace ADATIO;

namespace FF
{
  const int Nd = 4;

  //
  // Read formfactors
  //
  void read(FILE* bin, ComplexF& p)
  {
    float re,im;

    ADATUtil::bfread(&re, sizeof(float), 1, bin);
    ADATUtil::bfread(&im, sizeof(float), 1, bin);

    p = cmplx(Real(re), Real(im));
  }


  template<class T>
  void read(FILE* bin, Array<T>& p)
  {
    for(int i=0; i < p.size(); ++i)
      read(bin, p[i]);
  }


  //! Anonymous namespace
  namespace
  {
    std::ostream& operator<<(std::ostream& s, const ArrayInt& d)
    {
      if (d.size() > 0)
      {
	s << d[0];
	for(int i=1; i < d.size(); ++i)
	  s << " " << d[i];
      }

      return s;
    }


    std::ostream& operator<<(std::ostream& s, const ArrayDouble& d)
    {
      if (d.size() > 0)
      {
	s << d[0];
	for(int i=1; i < d.size(); ++i)
	  s << " " << d[i];
      }

      return s;
    }
  }


#define _DST_ 4
#define _ND_ 4

  /*#####################################################################################*/
  /*#####################################################################################*/

  void ReadBuildingBlocks_1(const std::string& BuildingBlocksFileName, BuildingBlocks_t& bar)
  {
    FILE* BBFile = NULL;
    unsigned short int Id;
    unsigned short int Version;
    unsigned short int Flavor;
    unsigned short int Contraction;
    unsigned short int NX;
    unsigned short int NY;
    unsigned short int NZ;
    unsigned short int NT;
    unsigned short int DecayDir;
    signed short int T1;
    signed short int T2;
    unsigned short int NMomPerms;
    unsigned short int MaxNLinks;
    unsigned short int NLinkPatterns;
    signed short int QX;
    signed short int QY;
    signed short int QZ;
    unsigned short int T;
    int MaxNLinkPatterns;
    char ContractionString[ 12 + 1 ]; /* 12 = strlen( "disconnected" ) > strlen( "connected" ) */
    int i;

    /*###################################################################################*/
    /* open building blocks file                                                         */
    /*###################################################################################*/

    BBFile = fopen( BuildingBlocksFileName.c_str(), "rb" );
    if( BBFile == NULL )
    {
      printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
      printf( "failed to read building blocks file %s\n", BuildingBlocksFileName.c_str() );
      exit(1);
    }

    /*###################################################################################*/
    /* read footer                                                                       */
    /*###################################################################################*/

    fseek( BBFile, - 10 * sizeof( unsigned short int ) - 5 * sizeof( signed short int ), SEEK_END );

    ADATUtil::bfread( & Flavor,        sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & Contraction,   sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & NX,            sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & NY,            sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & NZ,            sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & NT,            sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & T1,            sizeof( signed short int ),   1, BBFile );
    ADATUtil::bfread( & T2,            sizeof( signed short int ),   1, BBFile );
    ADATUtil::bfread( & MaxNLinks,     sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & NLinkPatterns, sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & QX,            sizeof( signed short int ),   1, BBFile );
    ADATUtil::bfread( & QY,            sizeof( signed short int ),   1, BBFile );
    ADATUtil::bfread( & QZ,            sizeof( signed short int ),   1, BBFile );
    ADATUtil::bfread( & Id,            sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & Version,       sizeof( unsigned short int ), 1, BBFile );

    assert( Id == 0 );
    assert( Version == 1 );

    if( MaxNLinks > 1024 )
    {
      printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
      puts( "no more than 1024 links is supported" );
      printf( "number of links = %i\n", MaxNLinks );
      exit(1);
    }

    if( Contraction > 1 )
    {
      printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
      puts( "only contraction types 0 (connected) and 1 (disconnected) are supported" );
      printf( "contraction = %i\n", Contraction );
      exit(1);
    }

    T = T2 - T1 + 1;

#if 0
    {
      if( Contraction == 0 ) strcpy( ContractionString, "connected" );
      if( Contraction == 1 ) strcpy( ContractionString, "disconnected" );

      printf( "BuildingBlocksFileName = %s\n", BuildingBlocksFileName.c_str() );
      printf( "Id            = %i\n", Id );
      printf( "Version       = %i\n", Version );
      printf( "Flavor        = %i\n", Flavor );
      printf( "Contraction   = %i = %s\n", Contraction, ContractionString );
      printf( "NX            = %i\n", NX );
      printf( "NY            = %i\n", NY );
      printf( "NZ            = %i\n", NZ );
      printf( "NT            = %i\n", NT );
      printf( "T1            = %i\n", T1 );
      printf( "T2            = %i\n", T2 );
      printf( "MaxNLinks     = %i\n", MaxNLinks );
      printf( "NLinkPatterns = %i\n", NLinkPatterns );
      printf( "QX            = %i\n", QX );
      printf( "QY            = %i\n", QY );
      printf( "QZ            = %i\n", QZ );
    }
#endif

    bar.param.latt_size.resize(Nd);
    bar.param.latt_size[0]   = NX;
    bar.param.latt_size[1]   = NY;
    bar.param.latt_size[2]   = NZ;
    bar.param.latt_size[3]   = NT;
    bar.param.t_srce    = T1;
    bar.param.t_sink    = T2;
    bar.param.DecayDir  = DecayDir;
    bar.param.links_max = MaxNLinks;
    bar.param.NMomPerms = NMomPerms;
    bar.param.Flavor    = Flavor;
    bar.param.GammaInsertion = 0;
    bar.param.Contraction    = Contraction;
    bar.param.sink_mom.resize(3);  // unfortunately, don't know the sink_mom
    bar.param.sink_mom[0]    = 0;
    bar.param.sink_mom[1]    = 0;
    bar.param.sink_mom[2]    = 0;

    Array< int > inser_mom(3);
    inser_mom[0] = QX;
    inser_mom[1] = QY;
    inser_mom[2] = QZ;


    /*###################################################################################*/
    /* read data                                                                         */
    /*###################################################################################*/

    rewind( BBFile );

    bar.links.resize(NLinkPatterns);

    for( int l = 0; l < NLinkPatterns; l ++ )
    {
      unsigned short int NLinks;
      unsigned short int mu[ 1024 ];

      ADATUtil::bfread( & NLinks, sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & mu, sizeof( unsigned short int ), NLinks, BBFile );

      bar.links[l].link_value.resize(NLinks);
      for(int i=0; i < bar.links[l].link_value.size(); ++i)
	bar.links[l].link_value[i] = mu[i];

      bar.links[l].gamma.resize(16);    // there are  Ns*Ns = 16 possible gamma matrices 
      for(int g=0; g < bar.links[l].gamma.size(); ++g)
      {
	bar.links[l].gamma[g].gamma_value = g;

	bar.links[l].gamma[g].momenta.resize(1);  // only 1 momenta here
	for(int m=0; m < bar.links[l].gamma[g].momenta.size(); ++m)
	{
	  bar.links[l].gamma[g].momenta[m].inser_mom = inser_mom;

	  bar.links[l].gamma[g].momenta[m].current.resize(NT);
	  read(BBFile, bar.links[l].gamma[g].momenta[m].current);
	}
      }
    }

    /*###################################################################################*/
    /* check that all which remains is the footer                                        */
    /*###################################################################################*/

    {
      /* "c" denotes "check" */
      unsigned short int cId;
      unsigned short int cVersion;
      unsigned short int cFlavor;
      unsigned short int cContraction;
      unsigned short int cNX;
      unsigned short int cNY;
      unsigned short int cNZ;
      unsigned short int cNT;
      signed short int cT1;
      signed short int cT2;
      unsigned short int cMaxNLinks;
      unsigned short int cNLinkPatterns;
      signed short int cQX;
      signed short int cQY;
      signed short int cQZ;

      char c;
      int n;
      int end;

      /* reproduce the footer */
      ADATUtil::bfread( & cFlavor,        sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cContraction,   sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cNX,            sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cNY,            sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cNZ,            sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cNT,            sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cT1,            sizeof( signed short int ),   1, BBFile );
      ADATUtil::bfread( & cT2,            sizeof( signed short int ),   1, BBFile );
      ADATUtil::bfread( & cMaxNLinks,     sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cNLinkPatterns, sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cQX,            sizeof( signed short int ),   1, BBFile );
      ADATUtil::bfread( & cQY,            sizeof( signed short int ),   1, BBFile );
      ADATUtil::bfread( & cQZ,            sizeof( signed short int ),   1, BBFile );
      ADATUtil::bfread( & cId,            sizeof( unsigned short int ), 1, BBFile );
      ADATUtil::bfread( & cVersion,       sizeof( unsigned short int ), 1, BBFile );

      /* verify file is at the end */
      n = ADATUtil::bfread( & c, sizeof( char ), 1, BBFile );
      end = feof( BBFile );

      if( ( ( cFlavor        == Flavor        ) &&
	    ( cContraction   == Contraction   ) &&
	    ( cNX            == NX            ) &&
	    ( cNY            == NY            ) &&
	    ( cNZ            == NZ            ) &&
	    ( cNT            == NT            ) &&
	    ( cT1            == T1            ) &&
	    ( cT2            == T2            ) &&
	    ( cMaxNLinks     == MaxNLinks     ) &&
	    ( cNLinkPatterns == NLinkPatterns ) &&
	    ( cQX            == QX            ) &&
	    ( cQY            == QY            ) &&
	    ( cQZ            == QZ            ) &&
	    ( cId            == Id            ) &&
	    ( cVersion       == Version       ) &&
	    ( n == 0 ) && ( end == 1 ) ) != 1 )
      {
	printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
	printf( "failed to read building blocks file %s correctly\n", BuildingBlocksFileName.c_str() );
	exit(1);
      }
    }

    /*###################################################################################*/
    /* finish                                                                            */
    /*###################################################################################*/

    assert( fclose( BBFile ) == 0 );

//    printf( "finished reading building blocks file %s\n", BuildingBlocksFileName.c_str() );
//    fflush( NULL );

    return;
  }



  /*#####################################################################################*/
  /*#####################################################################################*/

  void ReadBuildingBlocks_2(const std::string& BuildingBlocksFileName, BuildingBlocks_t& bar)
  {
//    std::cout << "Entering " << __func__ << std::endl;

    unsigned short int Id;
    unsigned short int Version;
    unsigned short int Flavor;
    unsigned short int Contraction;
    unsigned short int NX;
    unsigned short int NY;
    unsigned short int NZ;
    unsigned short int NT;
    unsigned short int DecayDir;
    signed short int T1;
    signed short int T2;
    unsigned short int NMomPerms;
    unsigned short int MaxNLinks;
    unsigned short int NLinkPatterns;
    signed short int GammaInsertion;
    signed short int QX;
    signed short int QY;
    signed short int QZ;
    signed short int PX;
    signed short int PY;
    signed short int PZ;
    unsigned short int T; 
    ADATUtil::n_uint32_t Checksum1, Checksum2;
    const signed short int   SeqSourceLen = 64;
    char SeqSource[SeqSourceLen];

    /*###################################################################################*/
    /* open building blocks file                                                         */
    /*###################################################################################*/
      
    BinaryFileReader BBFile(BuildingBlocksFileName);

    /*###################################################################################*/
    /* read footer                                                                       */
    /*###################################################################################*/
    
    BBFile.seekEnd(11 * sizeof( unsigned short int ) + 9 * sizeof( signed short int )
		   + 64 * sizeof( signed char ) + 1 * sizeof( ADATUtil::n_uint32_t ) );
      
    BBFile.read(Flavor);
    BBFile.read(Contraction);
    BBFile.readArray(SeqSource, sizeof( signed char ), SeqSourceLen);
    BBFile.read(GammaInsertion);
    BBFile.read(NX);
    BBFile.read(NY);
    BBFile.read(NZ);
    BBFile.read(NT);
    BBFile.read(DecayDir);
    BBFile.read(T1);
    BBFile.read(T2);
    BBFile.read(MaxNLinks);
    BBFile.read(NLinkPatterns);
    BBFile.read(QX);
    BBFile.read(QY);
    BBFile.read(QZ);
    BBFile.read(PX);
    BBFile.read(PY);
    BBFile.read(PZ);
    BBFile.read(Checksum1);
    BBFile.read(Id);
    BBFile.read(Version);

    assert( Id == 0 );
    assert( Version == 2 );

    if( Contraction > 1 )
    {
      printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
      puts( "only contraction types 0 (connected) and 1 (disconnected) are supported" );
      printf( "contraction = %i\n", Contraction );
      exit(1);
    }


#if 0
    {
      char ContractionString[ 12 + 1 ]; /* 12 = strlen( "disconnected" ) > strlen( "connected" ) */

      if( Contraction == 0 ) strcpy( ContractionString, "connected" );
      if( Contraction == 1 ) strcpy( ContractionString, "disconnected" );

      printf( "BuildingBlocksFileName = %s\n", BuildingBlocksFileName.c_str() );
      printf( "Id            = %i\n", Id );
      printf( "Version       = %i\n", Version );
      printf( "Flavor        = %i\n", Flavor );
      printf( "Contraction   = %i = %s\n", Contraction, ContractionString );
      printf( "SeqSource     = %s\n", SeqSource );
      printf( "GammaInsertion = %d\n", GammaInsertion );
      printf( "NX            = %i\n", NX );
      printf( "NY            = %i\n", NY );
      printf( "NZ            = %i\n", NZ );
      printf( "NT            = %i\n", NT );
      printf( "DecayDir      = %i\n", DecayDir );
      printf( "T1            = %i\n", T1 );
      printf( "T2            = %i\n", T2 );
      printf( "MaxNLinks     = %i\n", MaxNLinks );
      printf( "NLinkPatterns = %i\n", NLinkPatterns );
      printf( "QX            = %i\n", QX );
      printf( "QY            = %i\n", QY );
      printf( "QZ            = %i\n", QZ );
      printf( "PX            = %i\n", PX );
      printf( "PY            = %i\n", PY );
      printf( "PZ            = %i\n", PZ );
      printf( "Checksum      = %0x%x\n", Checksum1 );
    }
#endif

    bar.param.latt_size.resize(Nd);
    bar.param.latt_size[0]   = NX;
    bar.param.latt_size[1]   = NY;
    bar.param.latt_size[2]   = NZ;
    bar.param.latt_size[3]   = NT;
    bar.param.links_max = MaxNLinks;
    bar.param.t_srce    = T1;
    bar.param.t_sink    = T2;
    bar.param.DecayDir  = DecayDir;
    bar.param.links_max = MaxNLinks;
    bar.param.NMomPerms = NMomPerms;
    bar.param.Flavor    = Flavor;
    bar.param.GammaInsertion = GammaInsertion;
    bar.param.Contraction    = Contraction;
    bar.param.sink_mom.resize(3);
    bar.param.sink_mom[0]    = PX;
    bar.param.sink_mom[1]    = PY;
    bar.param.sink_mom[2]    = PZ;

    Array< int > inser_mom(3);
    inser_mom[0] = QX;
    inser_mom[1] = QY;
    inser_mom[2] = QZ;

    /*###################################################################################*/
    /* read data                                                                         */
    /*###################################################################################*/

    BBFile.rewind();

    bar.links.resize(NLinkPatterns);

    for( int l = 0; l < NLinkPatterns; l ++ )
    {
      unsigned short int NLinks;
      read(BBFile, NLinks);

      Array<unsigned short int> mu( NLinks );
      read(BBFile, mu, NLinks);

      bar.links[l].link_value.resize(NLinks);
      for(int i=0; i < bar.links[l].link_value.size(); ++i)
	bar.links[l].link_value[i] = mu[i];

      bar.links[l].gamma.resize(16);    // there are  Ns*Ns = 16 possible gamma matrices 
      for(int g=0; g < bar.links[l].gamma.size(); ++g)
      {
	bar.links[l].gamma[g].gamma_value = g;

	bar.links[l].gamma[g].momenta.resize(1);  // only 1 momenta here
	for(int m=0; m < bar.links[l].gamma[g].momenta.size(); ++m)
	{
	  bar.links[l].gamma[g].momenta[m].inser_mom = inser_mom;

	  bar.links[l].gamma[g].momenta[m].current.resize(NT);
	  read(BBFile, bar.links[l].gamma[g].momenta[m].current, NT);
	}
      }
    }

    /*###################################################################################*/
    /* check that all which remains is the footer                                        */
    /*###################################################################################*/

    {
      /* "c" denotes "check" */
      unsigned short int cId;
      unsigned short int cVersion;
      unsigned short int cFlavor;
      unsigned short int cContraction;
      unsigned short int cNX;
      unsigned short int cNY;
      unsigned short int cNZ;
      unsigned short int cNT;
      unsigned short int cDecayDir;
      signed short int cT1;
      signed short int cT2;
      unsigned short int cMaxNLinks;
      unsigned short int cNLinkPatterns;
      signed short int cGammaInsertion;
      signed short int cQX;
      signed short int cQY;
      signed short int cQZ;
      signed short int cPX;
      signed short int cPY;
      signed short int cPZ;
      ADATUtil::n_uint32_t cChecksum;
      char cSeqSource[SeqSourceLen];

      /* reproduce the footer */
      BBFile.read(cFlavor);
      BBFile.read(cContraction);
      BBFile.readArray(cSeqSource, sizeof( signed char ), SeqSourceLen );
      BBFile.read(cGammaInsertion);
      BBFile.read(cNX);
      BBFile.read(cNY);
      BBFile.read(cNZ);
      BBFile.read(cNT);
      BBFile.read(cDecayDir);
      BBFile.read(cT1);
      BBFile.read(cT2);
      BBFile.read(cMaxNLinks);
      BBFile.read(cNLinkPatterns);
      BBFile.read(cQX);
      BBFile.read(cQY);
      BBFile.read(cQZ);
      BBFile.read(cPX);
      BBFile.read(cPY);
      BBFile.read(cPZ);
      Checksum2 = BBFile.getChecksum();
      BBFile.read(cChecksum);
      BBFile.read(cId);
      BBFile.read(cVersion);

#if 0
      {
	char ContractionString[ 12 + 1 ]; /* 12 = strlen( "disconnected" ) > strlen( "connected" ) */

	if( cContraction == 0 ) strcpy( ContractionString, "connected" );
	if( cContraction == 1 ) strcpy( ContractionString, "disconnected" );

	printf( "cBuildingBlocksFileName = %s\n", BuildingBlocksFileName.c_str() );
	printf( "cId            = %i\n", cId );
	printf( "cVersion       = %i\n", cVersion );
	printf( "cFlavor        = %i\n", cFlavor );
	printf( "cContraction   = %i = %s\n", cContraction, ContractionString );
	printf( "cSeqSource     = %s\n", cSeqSource );
	printf( "cGammaInsertion = %d\n", cGammaInsertion );
	printf( "cNX            = %i\n", cNX );
	printf( "cNY            = %i\n", cNY );
	printf( "cNZ            = %i\n", cNZ );
	printf( "cNT            = %i\n", cNT );
	printf( "cDecayDir      = %i\n", cDecayDir );
	printf( "cT1            = %i\n", cT1 );
	printf( "cT2            = %i\n", cT2 );
	printf( "cMaxNLinks     = %i\n", cMaxNLinks );
	printf( "cNLinkPatterns = %i\n", cNLinkPatterns );
	printf( "cQX            = %i\n", cQX );
	printf( "cQY            = %i\n", cQY );
	printf( "cQZ            = %i\n", cQZ );
	printf( "cPX            = %i\n", cPX );
	printf( "cPY            = %i\n", cPY );
	printf( "cPZ            = %i\n", cPZ );
	printf( "cChecksum      = %0x%x\n", cChecksum );
      }
#endif

      // Check the checksum
      if (Checksum1 != cChecksum && Checksum2 != cChecksum)
      {
	printf("error at line %i in file %s\n", __LINE__, __FILE__);
	printf("computed checksum failed to match building blocks file %s\n", BuildingBlocksFileName.c_str() );
	exit(1);
      }

      if( ( ( cFlavor         == Flavor         ) &&
	    ( cContraction    == Contraction    ) &&
	    ( cGammaInsertion == GammaInsertion ) &&
	    ( cNX             == NX             ) &&
	    ( cNY             == NY             ) &&
	    ( cNZ             == NZ             ) &&
	    ( cNT             == NT             ) &&
	    ( cT1             == T1             ) &&
	    ( cT2             == T2             ) &&
	    ( cMaxNLinks      == MaxNLinks      ) &&
	    ( cNLinkPatterns  == NLinkPatterns  ) &&
	    ( cQX             == QX             ) &&
	    ( cQY             == QY             ) &&
	    ( cQZ             == QZ             ) &&
	    ( cPX             == PX             ) &&
	    ( cPY             == PY             ) &&
	    ( cPZ             == PZ             ) &&
	    ( cId             == Id             ) &&
	    ( cVersion        == Version        ) ) != 1 )
      {
	printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
	printf( "failed to read building blocks file %s correctly\n", BuildingBlocksFileName.c_str() );
	exit(1);
      }
    }

    /*###################################################################################*/
    /* finish                                                                            */
    /*###################################################################################*/

//  printf( "finished reading building blocks file %s\n", BuildingBlocksFileName.c_str() );

    return;
  }



  /*#####################################################################################*/
  /*#####################################################################################*/

  void ReadBuildingBlocks_3(const std::string& BuildingBlocksFileName, BuildingBlocks_t& bar)
  {
//    std::cout << "Entering " << __func__ << std::endl;

    unsigned short int Id;
    unsigned short int Version;
    unsigned short int Flavor;
    unsigned short int Contraction;
    unsigned short int NX;
    unsigned short int NY;
    unsigned short int NZ;
    unsigned short int NT;
    unsigned short int DecayDir;
    signed short int T1;
    signed short int T2;
    unsigned short int MaxNLinks;
    unsigned short int NLinkPatterns;
    unsigned short int NMomPerms;
    signed short int GammaInsertion;
    signed short int PX;
    signed short int PY;
    signed short int PZ;
    unsigned short int T; 
    ADATUtil::n_uint32_t Checksum1, Checksum2;
    const signed short int   SeqSourceLen = 64;
    char SeqSource[SeqSourceLen];

    /*###################################################################################*/
    /* open building blocks file                                                         */
    /*###################################################################################*/
      
    BinaryFileReader BBFile(BuildingBlocksFileName);

    /*###################################################################################*/
    /* read footer                                                                       */
    /*###################################################################################*/
    
    BBFile.seekEnd(12 * sizeof( unsigned short int ) + 6 * sizeof( signed short int )
		   + 64 * sizeof( signed char ) + 1 * sizeof( ADATUtil::n_uint32_t ) );
      
    BBFile.read(Flavor);
    BBFile.read(Contraction);
    BBFile.readArray(SeqSource, sizeof( signed char ), SeqSourceLen);
    BBFile.read(GammaInsertion);
    BBFile.read(NX);
    BBFile.read(NY);
    BBFile.read(NZ);
    BBFile.read(NT);
    BBFile.read(DecayDir);
    BBFile.read(T1);
    BBFile.read(T2);
    BBFile.read(MaxNLinks);
    BBFile.read(NLinkPatterns);
    BBFile.read(NMomPerms);
    BBFile.read(PX);
    BBFile.read(PY);
    BBFile.read(PZ);
    BBFile.read(Checksum1);
    BBFile.read(Id);
    BBFile.read(Version);

    assert( Id == 0 );
    assert( Version == 3 );

    if( Contraction > 1 )
    {
      printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
      puts( "only contraction types 0 (connected) and 1 (disconnected) are supported" );
      printf( "contraction = %i\n", Contraction );
      exit(1);
    }


#if 0
    {
      char ContractionString[ 12 + 1 ]; /* 12 = strlen( "disconnected" ) > strlen( "connected" ) */

      if( Contraction == 0 ) strcpy( ContractionString, "connected" );
      if( Contraction == 1 ) strcpy( ContractionString, "disconnected" );

      printf( "BuildingBlocksFileName = %s\n", BuildingBlocksFileName.c_str() );
      printf( "Id            = %i\n", Id );
      printf( "Version       = %i\n", Version );
      printf( "Flavor        = %i\n", Flavor );
      printf( "Contraction   = %i = %s\n", Contraction, ContractionString );
      printf( "SeqSource     = %s\n", SeqSource );
      printf( "GammaInsertion = %d\n", GammaInsertion );
      printf( "NX            = %i\n", NX );
      printf( "NY            = %i\n", NY );
      printf( "NZ            = %i\n", NZ );
      printf( "NT            = %i\n", NT );
      printf( "DecayDir      = %i\n", DecayDir );
      printf( "T1            = %i\n", T1 );
      printf( "T2            = %i\n", T2 );
      printf( "MaxNLinks     = %i\n", MaxNLinks );
      printf( "NLinkPatterns = %i\n", NLinkPatterns );
      printf( "NMomPerms     = %i\n", NMomPerms );
      printf( "PX            = %i\n", PX );
      printf( "PY            = %i\n", PY );
      printf( "PZ            = %i\n", PZ );
      printf( "Checksum      = %0x%x\n", Checksum1 );
    }
#endif

    bar.param.latt_size.resize(Nd);
    bar.param.latt_size[0]   = NX;
    bar.param.latt_size[1]   = NY;
    bar.param.latt_size[2]   = NZ;
    bar.param.latt_size[3]   = NT;
    bar.param.links_max = MaxNLinks;
    bar.param.t_srce    = T1;
    bar.param.t_sink    = T2;
    bar.param.DecayDir  = DecayDir;
    bar.param.links_max = MaxNLinks;
    bar.param.NMomPerms = NMomPerms;
    bar.param.Flavor    = Flavor;
    bar.param.GammaInsertion = GammaInsertion;
    bar.param.Contraction    = Contraction;
    bar.param.SeqSourceType  = SeqSource;
    bar.param.sink_mom.resize(3);
    bar.param.sink_mom[0]    = PX;
    bar.param.sink_mom[1]    = PY;
    bar.param.sink_mom[2]    = PZ;


    /*###################################################################################*/
    /* read data                                                                         */
    /*###################################################################################*/

    BBFile.rewind();

    bar.links.resize(NLinkPatterns);

    for( int l = 0; l < NLinkPatterns; l ++ )
    {
      unsigned short int NLinks;
      read(BBFile, NLinks);

      Array<unsigned short int> mu( NLinks );
      read(BBFile, mu, NLinks);

      bar.links[l].link_value.resize(NLinks);
      for(int i=0; i < bar.links[l].link_value.size(); ++i)
	bar.links[l].link_value[i] = mu[i];

//      std::cout << __func__ << ": NLinks = " << NLinks << "  link_value=" << bar.links[l].link_value << "\n";

      bar.links[l].gamma.resize(16);    // there are  Ns*Ns = 16 possible gamma matrices 
      for(int g=0; g < bar.links[l].gamma.size(); ++g)
      {
	bar.links[l].gamma[g].gamma_value = g;

	bar.links[l].gamma[g].momenta.resize(NMomPerms);  // varies from each BB file
	for(int m=0; m < bar.links[l].gamma[g].momenta.size(); ++m)
	{
	  Array<signed short int> Q( 3 );
	  read(BBFile, Q, Q.size());

	  bar.links[l].gamma[g].momenta[m].inser_mom.resize(Q.size());
	  for(int i=0; i < bar.links[l].gamma[g].momenta[m].inser_mom.size(); ++i)
	    bar.links[l].gamma[g].momenta[m].inser_mom[i] = Q[i];

	  bar.links[l].gamma[g].momenta[m].current.resize(NT);
	  read(BBFile, bar.links[l].gamma[g].momenta[m].current, NT);
	}
      }
    }

    /*###################################################################################*/
    /* check that all which remains is the footer                                        */
    /*###################################################################################*/

    {
      /* "c" denotes "check" */
      unsigned short int cId;
      unsigned short int cVersion;
      unsigned short int cFlavor;
      unsigned short int cContraction;
      unsigned short int cNX;
      unsigned short int cNY;
      unsigned short int cNZ;
      unsigned short int cNT;
      unsigned short int cDecayDir;
      signed short int cT1;
      signed short int cT2;
      unsigned short int cMaxNLinks;
      unsigned short int cNLinkPatterns;
      unsigned short int cNMomPerms;
      signed short int cGammaInsertion;
      signed short int cPX;
      signed short int cPY;
      signed short int cPZ;
      ADATUtil::n_uint32_t cChecksum;
      char cSeqSource[SeqSourceLen];

      /* reproduce the footer */
      BBFile.read(cFlavor);
      BBFile.read(cContraction);
      BBFile.readArray(cSeqSource, sizeof( signed char ), SeqSourceLen );
      BBFile.read(cGammaInsertion);
      BBFile.read(cNX);
      BBFile.read(cNY);
      BBFile.read(cNZ);
      BBFile.read(cNT);
      BBFile.read(cDecayDir);
      BBFile.read(cT1);
      BBFile.read(cT2);
      BBFile.read(cMaxNLinks);
      BBFile.read(cNLinkPatterns);
      BBFile.read(cNMomPerms);
      BBFile.read(cPX);
      BBFile.read(cPY);
      BBFile.read(cPZ);
      Checksum2 = BBFile.getChecksum();
      BBFile.read(cChecksum);
      BBFile.read(cId);
      BBFile.read(cVersion);

#if 0
      {
	char ContractionString[ 12 + 1 ]; /* 12 = strlen( "disconnected" ) > strlen( "connected" ) */

	if( cContraction == 0 ) strcpy( ContractionString, "connected" );
	if( cContraction == 1 ) strcpy( ContractionString, "disconnected" );

	printf( "cBuildingBlocksFileName = %s\n", BuildingBlocksFileName.c_str() );
	printf( "cId            = %i\n", cId );
	printf( "cVersion       = %i\n", cVersion );
	printf( "cFlavor        = %i\n", cFlavor );
	printf( "cContraction   = %i = %s\n", cContraction, ContractionString );
	printf( "cSeqSource     = %s\n", cSeqSource );
	printf( "cGammaInsertion = %d\n", cGammaInsertion );
	printf( "cNX            = %i\n", cNX );
	printf( "cNY            = %i\n", cNY );
	printf( "cNZ            = %i\n", cNZ );
	printf( "cNT            = %i\n", cNT );
	printf( "cDecayDir      = %i\n", cDecayDir );
	printf( "cT1            = %i\n", cT1 );
	printf( "cT2            = %i\n", cT2 );
	printf( "cMaxNLinks     = %i\n", cMaxNLinks );
	printf( "cNLinkPatterns = %i\n", cNLinkPatterns );
	printf( "cNMomPerms     = %i\n", cNMomPerms );
	printf( "cPX            = %i\n", cPX );
	printf( "cPY            = %i\n", cPY );
	printf( "cPZ            = %i\n", cPZ );
	printf( "cChecksum      = %0x%x\n", cChecksum );
      }
#endif

      // Check the checksum
      if (Checksum1 != cChecksum && Checksum2 != cChecksum)
      {
	printf("error at line %i in file %s\n", __LINE__, __FILE__);
	printf("computed checksum failed to match building blocks file %s\n", BuildingBlocksFileName.c_str() );
	exit(1);
      }

      if( ( ( cFlavor         == Flavor         ) &&
	    ( cContraction    == Contraction    ) &&
	    ( cGammaInsertion == GammaInsertion ) &&
	    ( cNX             == NX             ) &&
	    ( cNY             == NY             ) &&
	    ( cNZ             == NZ             ) &&
	    ( cNT             == NT             ) &&
	    ( cT1             == T1             ) &&
	    ( cT2             == T2             ) &&
	    ( cMaxNLinks      == MaxNLinks      ) &&
	    ( cNLinkPatterns  == NLinkPatterns  ) &&
	    ( cPX             == PX             ) &&
	    ( cPY             == PY             ) &&
	    ( cPZ             == PZ             ) &&
	    ( cId             == Id             ) &&
	    ( cVersion        == Version        ) ) != 1 )
      {
	printf( "error at line %i in file %s\n", __LINE__, __FILE__ );
	printf( "failed to read building blocks file %s correctly\n", BuildingBlocksFileName.c_str() );
	exit(1);
      }
    }

    /*###################################################################################*/
    /* finish                                                                            */
    /*###################################################################################*/

//  printf( "finished reading building blocks file %s\n", BuildingBlocksFileName.c_str() );

    return;
  }



  // Mother of all readers
  void read(const std::string& BuildingBlocksFileName, BuildingBlocks_t& bar)
  {
    FILE* BBFile = NULL;
    unsigned short int Id;
    unsigned short int Version;
    const unsigned short int MaxVersion = 3;

    /*###################################################################################*/
    /* open building blocks file                                                         */
    /*###################################################################################*/

    BBFile = fopen( BuildingBlocksFileName.c_str(), "rb" );
    if( BBFile == NULL )
    {
      printf( "\nfailed to open building blocks file %s\n\n", BuildingBlocksFileName.c_str() );
      exit( 1 );
    }

    /*###################################################################################*/
    /* check id and version                                                              */
    /*###################################################################################*/

    fseek( BBFile, - 2 * sizeof( unsigned short int ), SEEK_END );

    ADATUtil::bfread( & Id, sizeof( unsigned short int ), 1, BBFile );
    ADATUtil::bfread( & Version, sizeof( unsigned short int ), 1, BBFile );

    bar.version = Version;

    assert( fclose( BBFile ) == 0 );

    if( Id != 0 )
    {
      printf( "\nThe file does not seem to contain building blocks: %s\n", BuildingBlocksFileName.c_str() );
      printf( "The id for the file is %i and it should be 0.\n\n", Id );
      exit( 1 );
    }

    if( Version > MaxVersion )
    {
      puts( "" );
      printf( "The file %s seems to correspond to an unknown version %i.\n", BuildingBlocksFileName.c_str(), Version );
      printf( "The known versions are 0, ..., %i.\n\n", MaxVersion );
      exit( 1 );
    }

    /*###################################################################################*/
    /* switch to appropriate version                                                     */
    /*###################################################################################*/

    //  std::cout << "bb version = " << Version << std::endl;

    switch( Version )
    {
    case 0:
//      ReadBuildingBlocks_0( BuildingBlocksFileName, bar );
      std::cerr << "Version 0 not supported" << std::endl;
      exit(1);
      break;
    case 1:
      ReadBuildingBlocks_1( BuildingBlocksFileName, bar );
      break;
    case 2:
      ReadBuildingBlocks_2( BuildingBlocksFileName, bar );
      break;
    case 3:
      ReadBuildingBlocks_3( BuildingBlocksFileName, bar );
      break;
    default:
      puts( "\nconfused about version of building blocks\n" );
      exit( 1 );
    }

  }


  /*#####################################################################################*/
  /*#####################################################################################*/

  void WriteBuildingBlocks_3(const std::string& BuildingBlocksFileName, 
			     const BuildingBlocks_t& bar)
  {
//    std::cout << "Entering " << __func__ << std::endl;

    /*###################################################################################*/
    /* open building blocks file                                                         */
    /*###################################################################################*/
      
    BinaryFileWriter BBFile(BuildingBlocksFileName);

    /*###################################################################################*/
    /* write data                                                                         */
    /*###################################################################################*/

    const unsigned short int NX = bar.param.latt_size[0];
    const unsigned short int NY = bar.param.latt_size[1];
    const unsigned short int NZ = bar.param.latt_size[2];
    const unsigned short int NT = bar.param.latt_size[3];

    for( int l = 0; l < bar.links.size(); l ++ )
    {
      unsigned short int NLinks = bar.links[l].link_value.size();
      write(BBFile, NLinks);

      Array<unsigned short int> mu( NLinks );
      for(int i=0; i < bar.links[l].link_value.size(); ++i)
	mu[i] = bar.links[l].link_value[i];

      write(BBFile, mu, NLinks);

//    std::cout << __func__ << ": NLinks = " << NLinks << "  link_value=" << bar.links[l].link_value << "\n";

      for(int g=0; g < bar.links[l].gamma.size(); ++g)
      {
	for(int m=0; m < bar.links[l].gamma[g].momenta.size(); ++m)
	{
//        std::cout << __func__ << ": inser_mom = " << bar.links[l].gamma[g].momenta[m].inser_mom << std::endl;

	  Array<signed short int> Q( bar.links[l].gamma[g].momenta[m].inser_mom.size() );
	  for(int i=0; i < bar.links[l].gamma[g].momenta[m].inser_mom.size(); ++i)
	    Q[i] = bar.links[l].gamma[g].momenta[m].inser_mom[i];

	  write(BBFile, Q, Q.size());
	  write(BBFile, bar.links[l].gamma[g].momenta[m].current, NT);
	}
      }
    }

//    std::cout << __func__ << "write trailer" << std::endl;

    const unsigned short int Id = 0;  // indicates building blocks
    const unsigned short int Version = 3;  // building blocks version
    const unsigned short int Contraction = 0;  // 0 indicates connected diagram
    const signed short int   PX = bar.param.sink_mom[0];
    const signed short int   PY = bar.param.sink_mom[1];
    const signed short int   PZ = bar.param.sink_mom[2];
    const signed short int   T1 = bar.param.t_srce;
    const signed short int   T2 = bar.param.t_sink;
    const signed short int   SeqSourceLen = 64;
    std::string SeqSource = bar.param.SeqSourceType;
    SeqSource.resize(SeqSourceLen, 0);

    const signed short int Flavor = bar.param.Flavor;  // currently assumes u and d are given as f=0 and f=1
    const signed short int GammaInsertion = bar.param.GammaInsertion;
    const unsigned short int DecayDir = bar.param.DecayDir;
    const unsigned short int MaxNLinks = bar.param.links_max;
    const unsigned short int NLinkPatterns = bar.links.size();
    const unsigned short int NMomPerms = bar.param.NMomPerms;

    // possibly specific to this version
    BBFile.write( Flavor );
    BBFile.write( Contraction );
    BBFile.writeArray( SeqSource.data(), 1, SeqSourceLen );
    BBFile.write( GammaInsertion );
    BBFile.write( NX );
    BBFile.write( NY );
    BBFile.write( NZ );
    BBFile.write( NT );
    BBFile.write( DecayDir );
    BBFile.write( T1 );
    BBFile.write( T2 );
    BBFile.write( MaxNLinks );
    BBFile.write( NLinkPatterns );
    BBFile.write( NMomPerms );
    BBFile.write( PX );
    BBFile.write( PY );
    BBFile.write( PZ );
    BBFile.write( BBFile.getChecksum() );
    // generic to any building blocks file
    BBFile.write( Id );
    BBFile.write( Version );
    // close file
    BBFile.close();

    /*###################################################################################*/
    /* finish                                                                            */
    /*###################################################################################*/

//  printf( "finished writing building blocks file %s\n", BuildingBlocksFileName.c_str() );

    return;
  }


  // Writer
  void write(const std::string& BuildingBlocksFileName, const BuildingBlocks_t& bar)
  {
    /*###################################################################################*/
    /* Write only in the latest version                                                  */
    /*###################################################################################*/
    WriteBuildingBlocks_3( BuildingBlocksFileName, bar );
  }


  // Thin out a BB structure
  void thinBB(BuildingBlocks_t& new_bar, 
	      int new_links_max,
	      const BuildingBlocks_t& old_bar)
  {
    // Copying the old BB file will also copy along all those params
    new_bar = old_bar;

    // Handle oddball cases
    if (new_links_max == old_bar.param.links_max)
    {
      // NOP
      return;
    }
    else if (new_links_max > old_bar.param.links_max)
    {
      std::cerr << __func__ << ": the specified links_max is greater than in the original file" 
		<< std::endl;
      exit(1);
    }
    else if (new_links_max < 0)
    {
      std::cerr << __func__ << ": invalid links_max = " << new_links_max << std::endl;
      exit(1);
    }

    // Now thin the puppy - links_max goes down or stays the same
    std::list<BuildingBlocks_t::Links_t> links_list;

    for( int l = 0; l < old_bar.links.size(); l ++ )
    {
      // Here is where the real work happens - if the link_value list is
      // bigger than the desired size, we skip it. We accumulate desired
      // ones onto a list. Only at the end do we know the final list size
      if (old_bar.links[l].link_value.size() > new_links_max)
	continue;

      BuildingBlocks_t::Links_t new_links;
      new_links = old_bar.links[l];

      links_list.push_back(new_links);
    }

    // Rebuild the final structure
    new_bar.param.links_max = new_links_max;
    new_bar.links.resize(links_list.size());
    int l = 0;
    for(std::list<BuildingBlocks_t::Links_t>::const_iterator ll_ptr = links_list.begin();
	ll_ptr != links_list.end();
	++ll_ptr, ++l)
    {
      new_bar.links[l] = *ll_ptr;
    }

  }

} // namespace FF
