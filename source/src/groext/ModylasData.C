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
#include "ModylasData.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include <assert.h>
#include <float.h>
#include <math.h>

using namespace std;

bool ModylasData::ReadMDXYZ()
{
  String filename = m_filebase + ".mdxyz";

  ifstream ifs( filename.c_str() );
  if( ifs.fail() ) {
    cerr << "Info) Could not open " + filename + ". Skip reading." << endl;
    return false;
  }
  //cout << "Info) Reading " + filename + " ..." << endl;

  String strall( string(
    ( istreambuf_iterator<char>(ifs) ),
    ( istreambuf_iterator<char>()    ) ) );

  bool b = m_mdxyz_root_tag.ParseFromString( strall );
  if( !b ) {
    cerr << "Error) Could not parse file contents : " << filename << endl;
    return false;
  }

  return true;
}
  

bool ModylasData::ReadMDFF()
{
  String filename = m_filebase + ".mdff";

  ifstream ifs( filename.c_str() );
  if( ifs.fail() ) {
    cerr << "Info) Could not open " + filename + ". Skip reading." << endl;
    return false;
  }
  //cout << "Info) Reading " + filename + " ..." << endl;

  String strall( string(
    ( istreambuf_iterator<char>(ifs) ),
    ( istreambuf_iterator<char>()    ) ) );

  bool b = m_mdff_root_tag.ParseFromString( strall );
  if( !b ) {
    cerr << "Error) Could not parse mdff file contents." << endl;
    return false;
  }

  return true;
}


bool ModylasData::WriteMDXYZ() const
{
  String filename = m_filebase + ".mdxyz";
  cout << "Info) Writing " + filename + " ..." << endl;

  ofstream ofs( filename.c_str() );
  if( ofs.fail() ) {
    cerr << "ERROR) Failed to open the file for writing : " + filename << endl;
    return false;
  }
  ofs << m_mdxyz_root_tag.WriteDataString();
  ofs.close();

  return true;
}

bool ModylasData::WriteMDFF() const
{
  String filename = m_filebase + ".mdff";
  cout << "Info) Writing " + filename + " ..." << endl;

  ofstream ofs( filename.c_str() );
  if( ofs.fail() ) {
    cerr << "ERROR) Failed to open the file for writing : " + filename << endl;
    return false;
  }
  ofs << m_mdff_root_tag.WriteDataString();
  ofs.close();

  return true;
}

