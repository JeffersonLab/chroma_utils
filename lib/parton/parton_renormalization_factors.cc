// $Id: parton_renormalization_factors.cc,v 2.0 2008/12/05 04:43:38 edwards Exp $
/*! \file
 *  \brief Renormalization factors
 */

#define Q 2 // GeV
#define _N_C 3

#define Z( me, beta, a ) 1 + \
( 2 * _N_C / beta ) / ( 16 * _PI_ * _PI_ ) * ( _N_C * _N_C - 1 ) / ( 2 * _N_C ) * \
( CONCAT2( GammaMSBar, me ) * log( Q * Q * a * a ) - ( CONCAT2( BLat, me ) - CONCAT2( BMSBar, me ) ) )

// Dmitri Dolgov did use the log term

#define Z_DD( me, beta, a ) 1 + \
( 2 * _N_C / beta ) / ( 16 * _PI_ * _PI_ ) * ( _N_C * _N_C - 1 ) / ( 2 * _N_C ) * \
( 0                                                - ( CONCAT2( BLat, me ) - CONCAT2( BMSBar, me ) ) )

//#####################################################################################
//#####################################################################################

void SetRenormalizationFactors( const float beta, const float a )
{
    Z_x_q_a = Z(   _x_q_a, beta, a );
    Z_x_q_b = Z(   _x_q_b, beta, a );
   Z_xx_q   = Z(  _xx_q,   beta, a );
  Z_xxx_q   = Z( _xxx_q,   beta, a );

   Z_1_Delta_q   = Z(  _1_Delta_q,   beta, a );
   Z_x_Delta_q_a = Z(  _x_Delta_q_a, beta, a );
   Z_x_Delta_q_b = Z(  _x_Delta_q_b, beta, a );
  Z_xx_Delta_q   = Z( _xx_Delta_q,   beta, a );

  Z_1_delta_q = Z( _1_delta_q, beta, a );
  Z_x_delta_q = Z( _x_delta_q, beta, a );

  Z_d1_q = Z( _d1_q, beta, a );
  Z_d2_q = Z( _d2_q, beta, a );
}

//#####################################################################################
//#####################################################################################

void SetRenormalizationFactors_DD( const float beta, const float a )
{
    Z_x_q_a = Z_DD(   _x_q_a, beta, a );
    Z_x_q_b = Z_DD(   _x_q_b, beta, a );
   Z_xx_q   = Z_DD(  _xx_q,   beta, a );
  Z_xxx_q   = Z_DD( _xxx_q,   beta, a );

   Z_1_Delta_q   = Z_DD(  _1_Delta_q,   beta, a );
   Z_x_Delta_q_a = Z_DD(  _x_Delta_q_a, beta, a );
   Z_x_Delta_q_b = Z_DD(  _x_Delta_q_b, beta, a );
  Z_xx_Delta_q   = Z_DD( _xx_Delta_q,   beta, a );

  Z_1_delta_q = Z_DD( _1_delta_q, beta, a );
  Z_x_delta_q = Z_DD( _x_delta_q, beta, a );

  Z_d1_q = Z_DD( _d1_q, beta, a );
  Z_d2_q = Z_DD( _d2_q, beta, a );
}

