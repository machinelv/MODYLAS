/*----------------------------------------------------------------------
*MODYLAS ver. 1.1.0
*
*Copyright (c) 2014-2019 Nagoya University
*              2020-2023 The University of Tokyo
*
*Released under the MIT license.
*see https://opensource.org/licenses/MIT
*-----------------------------------------------------------------------
*MODYLAS Developers:
*Yoshimichi Andoh, Kazushi Fujimoto, Tatsuya Sakashita, Noriyuki Yoshii, 
*Zhiye Tang, Jiachao Zhang, Yuta Asano, Ryo Urano, Tetsuro Nagai, 
*Atsushi Yamada, Hidekazu Kojima, Kensuke Iwahashi, Fumiyasu Mizutani, 
*Shin-ichi Ichikawa, and Susumu Okazaki.
*----------------------------------------------------------------------*/
/**
* @file
  @brief  Utilities to deal with MODYLAS data format
* @author  Tetsuro Nagai
*/
#ifndef MODYLAS_UTIL_H
#define MODYLAS_UTIL_H

#include "String.h"
#include <vector>
#include <math.h>

using namespace std;

// definition of utility types

template<int N>
class intn
{
public:
  intn(){
    for( int j = 0; j < N; ++j ) {
      m_i[j] = 0;
    }
  };
  intn( int i0, int i1 );
  intn( int i0, int i1, int i2 );
  intn( int i0, int i1, int i2, int i3 );

  int& operator[]( int i ) {
    return m_i[i];
  }
  int operator[]( int i ) const {
    return m_i[i];
  }
  intn<N>& operator=( const intn<N>& i ) {
    for( int j = 0; j < N; ++j ) {
      m_i[j] = i[j];
    }
    return *this;
  }
  bool operator>( const intn<N>& i ) const {
    for( int j = 0; j < N; ++j ) {
      if( m_i[j] > i[j] ) return true;
      if( m_i[j] < i[j] ) return false;
    }
    return false;
  }
  bool operator<( const intn<N>& i ) const {
    for( int j = 0; j < N; ++j ) {
      if( m_i[j] < i[j] ) return true;
      if( m_i[j] > i[j] ) return false;
    }
    return false;
  }
private:
  int m_i[N];
};

typedef intn<2> int2;
typedef intn<3> int3;
typedef intn<4> int4;

struct String3
{
  String s0;
  String s1;
  String s2;
};

struct String4
{
  String s0;
  String s1;
  String s2;
  String s3;
};

// class for utility functions
class ModylasUtil
{
public:
  static String GetElementSymbol( int i );
  static int GetAtomicNumber( const String& str );
  static String GetMassElement( double mass );
  static int GetMassAtomicNumber( double mass );
};

const double PI         = 3.1415926535897932384626;
const double KCAL_TO_KJ = 4.184;
const double KJ_TO_KCAL = 1.0 / KCAL_TO_KJ;
const double A_TO_NM    = 0.1;
const double NM_TO_A    = 1.0 / A_TO_NM;
const double SIGMA_TO_R = pow(2.0, 1.0/6.0);
const double R_TO_SIGMA = 1.0 / SIGMA_TO_R;
const double M_TO_NM    = 1.0e+9;
const double PS_TO_S    = 1.0e-12;

#endif //MODYLAS_UTIL_H


