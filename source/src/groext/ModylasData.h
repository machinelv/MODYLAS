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
#ifndef MODYLAS_DATA_H
#define MODYLAS_DATA_H

#include "ModylasUtil.h"
#include "ModylasTagData.h"

using namespace std;

class ModylasData
{
public:
  ModylasData();
  ModylasData(const String& f) : m_filebase(f) {}

  bool ReadMDXYZ();
  bool ReadMDFF ();

  bool WriteMDXYZ() const;
  bool WriteMDFF () const;

public:
  ModylasTagData m_mdxyz_root_tag;
  ModylasTagData m_mdff_root_tag;

private:
  String m_filebase;
};


#endif //MODYLAS_DATA_H


