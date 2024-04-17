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
#include "ModylasUtil.h"

#include <sstream>
#include <cmath>

using namespace std;

// partial specialization of intn<> constructor
template<>
intn<2>::intn( int i0, int i1 ) {
  m_i[0] = i0;
  m_i[1] = i1;
};

template<>
intn<3>::intn( int i0, int i1, int i2 ) {
  m_i[0] = i0;
  m_i[1] = i1;
  m_i[2] = i2;
};

template<>
intn<4>::intn( int i0, int i1, int i2, int i3 ) {
  m_i[0] = i0;
  m_i[1] = i1;
  m_i[2] = i2;
  m_i[3] = i3;
};


static const int NELEM = 112;
static const String ELEMENT[ NELEM ] = {
  "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
  "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
  "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", 
  "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
  "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
  "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
  "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
  "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
  "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
  "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
  "Ds", "Rg"
};
static const double ELEMENT_MASS[ NELEM ] = {
    0.000000,   1.007825,   4.002603,   7.016005,   9.012182,
   11.009305,  12.000000,  14.003074,  15.994915,  18.998403,
   19.992440,  22.989769,  23.985042,  26.981539,  27.976927,
   30.973762,  31.972071,  34.968853,  39.962383,  38.963707,
   39.962591,  44.955912,  47.947946,  50.943960,  55.934938,
   54.938045,  55.934938,  58.933195,  57.935343,  62.929598,
   63.929142,  68.925574,  73.921178,  74.921597,  79.916521,
   78.918337,  83.911507,  84.911790,  87.905612,  88.905848,
   89.904704,  92.906378,  97.905408,  98.906255, 101.904349,
  102.905504, 105.903486, 106.905097, 113.903359, 114.903878,
  119.902195, 120.903816, 129.906224, 126.904473, 131.904154,
  132.905452, 137.905247, 138.906353, 139.905439, 140.907653,
  143.910083, 144.912749, 151.919732, 152.921230, 157.924104,
  158.925347, 163.929175, 164.930322, 165.930293, 168.934213,
  173.938862, 174.940772, 179.946550, 180.947996, 183.950931,
  186.955753, 191.961481, 192.962926, 197.967893, 196.966569,
  201.970643, 204.974428, 207.976652, 208.980399, 209.982874,
  209.987148, 222.017578, 223.019736, 226.025410, 227.027752,
  232.038055, 231.035884, 238.050788, 237.048173, 239.052163,
  243.061381, 247.070354, 247.070307, 252.081626, 252.082979,
  257.095105, 258.098431, 259.101031, 262.109634, 267.121529,
  268.125445, 271.133472, 272.138032, 277.149841, 276.151156,
  281.162061, 280.164473
};

String ModylasUtil::GetElementSymbol( int i )
{
  if( i < 1 || NELEM < i ) {
    return "";
  }
  return ELEMENT[i];
}

int ModylasUtil::GetAtomicNumber( const String& str )
{
  if( str.size() == 0 ) {
    return 0;
  }
  String elem = str;
  if( str.size() > 1 ) {
    elem = str.substr(0,1).ToUpper() + str.substr(1,1).ToLower();
  }
  for( int i = 1; i < NELEM; ++i ) {
    if( elem == ELEMENT[i] ) return i;  
  }
  return 0;
}

String ModylasUtil::GetMassElement( double mass )
{
  int i = GetMassAtomicNumber(mass);
  if( i == 0 ) {
    return 0;
  }
  else {
    return ELEMENT[i];
  }
}

int ModylasUtil::GetMassAtomicNumber( double mass )
{
  double tol = 0.1;
  for(int i = 1; i < NELEM; ++i ) {
    if( abs( ELEMENT_MASS[i] - mass ) < tol ) {
      return i;
    }
  }
  return 0;
}

