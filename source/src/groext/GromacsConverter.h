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
  @brief  Utilities for converting from GROMACS inputs
* @author  Tetsuro Nagai
*/

#ifndef GROMACS_CONVERTER_H
#define GROMACS_CONVERTER_H

#include "ModylasData.h"
#include "GromacsData.h"
#include <map>

using namespace std;

enum EConstraints {
  CONSTRAINTS_NONE     = 1,
  CONSTRAINTS_HBONDS   = 2,
  CONSTRAINTS_ALLBONDS = 3,
};


enum class ForceFieldType {
  charmm,
  opls,
  gaff,
  unsupported
};

class GromacsConverter
{
public:
  static bool Gromacs2Modylas( const GromacsData& gmx_data, ModylasData& mod_data,
                               EConstraints eConstraints = CONSTRAINTS_HBONDS);
  static bool Modylas2Gromacs( const ModylasData& mod_data, GromacsData& gmx_data );
  static bool ConvertAngleToCell(double a, double b, double c,
                                 double alpha, double beta, double gamma,
                                 GROCELL& cell);
  static bool ConvertMdpToMddef( const String& strSessionName );

private:
  static int SearchDihedralType(
    const TOPDIHEDRAL& dihedral,
    const vector<TOPDIHEDRALTYPE> dihedraltypes,
    const vector<TOPATOM>& rAtoms,
    const map<String, String>& mapAtypeBtype);
  static vector<int> SearchDihedralTypes(
    const TOPDIHEDRAL& dihedral,
    const vector<TOPDIHEDRALTYPE> dihedraltypes,
    const vector<TOPATOM>& rAtoms,
    const map<String, String>& mapAtypeBtype);
};


#endif //GROMACS_CONVERTER_H


