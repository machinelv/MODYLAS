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
#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "ModylasData.h"
#include "GromacsData.h"
#include "GromacsConverter.h"

using namespace std;

char* c_strMdxyz = NULL;
char* c_strMdff  = NULL;
char* c_strGro   = NULL;
char* c_strTop   = NULL;

extern "C"
char* c_get_mdxyz_ptr()
{
  return c_strMdxyz;
}

extern "C"
char* c_get_mdff_ptr()
{
  return c_strMdff;
}

extern "C"
char* c_get_gro_ptr()
{
  return c_strGro;
}

extern "C"
char* c_get_top_ptr()
{
  return c_strTop;
}

extern "C"
void c_convert_gmx_modylas(const char* pcSessionName,
                           const char* pcIncludePath,
                           const char* pcIncKeys,
                           int iReadMDP, int iWriteMddef, int iConstraints)
{
  c_strMdxyz = NULL;
  c_strMdff  = NULL;

  string strSessionName = pcSessionName;
  string strIncludePath = pcIncludePath;
  string strIncKeys     = pcIncKeys;
  GromacsData gmx_data(strSessionName, strIncludePath);

  // Read GRO
  if( !gmx_data.ReadGRO() ) {
    return;
  }

  if( iReadMDP ) {
    String mdpfile = strSessionName + ".mdp";

    ifstream ifs( mdpfile.c_str() );
    if( ifs.fail() ) {
      cerr << "Error) Could not open " + mdpfile + ". Skip reading." << endl;
      return;
    }

    String line;
    while( getline(ifs, line) ) {
      vector<String> tokens = line.Split(";");
      if( tokens.size() == 0 ) continue;
      line = tokens[0];
      tokens = line.Split();
      if( tokens.size() < 3     ) continue;
      if( tokens[0] != "define" ) continue;
      if( tokens[1] != "="      ) continue;
      for(unsigned int i = 2; i < tokens.size(); ++i) {
        String define = tokens[i];
        define = define.ReplaceAll("-D", "");
        strIncKeys += " " + define;
      }
    }

    if( iWriteMddef ) {
      GromacsConverter::ConvertMdpToMddef(strSessionName);
    }
  }

  // Read TOP
  if( !gmx_data.ReadTOP(strIncKeys) ) {
    return;
  }

  ModylasData mod_data(strSessionName);

  if( !GromacsConverter::Gromacs2Modylas( gmx_data, mod_data, (EConstraints)iConstraints ) ) {
    return;
  }
  
  string strMdxyz = mod_data.m_mdxyz_root_tag.WriteDataString();
  string strMdff  = mod_data.m_mdff_root_tag .WriteDataString();
  string strGro   = String::ReadFileString(strSessionName + ".gro");
  string strTop   = String::ReadFileString(strSessionName + ".top");

  int n1 = strMdxyz.length();
  int n2 = strMdff .length();
  int n3 = strGro  .length();
  int n4 = strTop  .length();
  c_strMdxyz = (char*)malloc((n1+1) * sizeof(char));
  c_strMdff  = (char*)malloc((n2+1) * sizeof(char));
  c_strGro   = (char*)malloc((n3+1) * sizeof(char));
  c_strTop   = (char*)malloc((n4+1) * sizeof(char));
  strcpy(c_strMdxyz, strMdxyz.c_str());
  strcpy(c_strMdff , strMdff .c_str());
  strcpy(c_strGro  , strGro  .c_str());
  strcpy(c_strTop  , strTop  .c_str());
}

extern "C"
void c_convert_modylas_gmx(const char* pcSessionName)
{
  c_strMdxyz = NULL;
  c_strMdff  = NULL;

  string strSessionName = pcSessionName;
  ModylasData mod_data(strSessionName);

  // Read GRO
  if( !mod_data.ReadMDXYZ() ) {
    return;
  }

  // Read TOP
  if( !mod_data.ReadMDFF() ) {
    return;
  }

  GromacsData gmx_data(strSessionName);

  if( !GromacsConverter::Modylas2Gromacs( mod_data, gmx_data ) ) {
    return;
  }

  string strMdxyz = String::ReadFileString(strSessionName + ".mdxyz");
  string strMdff  = String::ReadFileString(strSessionName + ".mdff");
  string strGro   = gmx_data.WriteGROString();
  string strTop   = gmx_data.WriteTOPString();

  int n1 = strMdxyz.length();
  int n2 = strMdff .length();
  int n3 = strGro  .length();
  int n4 = strTop  .length();
  c_strMdxyz = (char*)malloc((n1+1) * sizeof(char));
  c_strMdff  = (char*)malloc((n2+1) * sizeof(char));
  c_strGro   = (char*)malloc((n3+1) * sizeof(char));
  c_strTop   = (char*)malloc((n4+1) * sizeof(char));
  strcpy(c_strMdxyz, strMdxyz.c_str());
  strcpy(c_strMdff , strMdff .c_str());
  strcpy(c_strGro  , strGro  .c_str());
  strcpy(c_strTop  , strTop  .c_str());
}

extern "C"
void c_update_gro_positions(double xyz[], double vel[], double cell[3], size_t nSize)
{
  double a     = cell[0] * M_TO_NM;
  double b     = cell[1] * M_TO_NM;
  double c     = cell[2] * M_TO_NM;
  double alpha = cell[3];
  double beta  = cell[4];
  double gamma = cell[5];
  GROCELL grocell;
  
  if( !GromacsConverter::ConvertAngleToCell(a, b, c, alpha, beta, gamma, grocell) )
    return;

  double move[3];
  move[0] = ( grocell.v1[0] + grocell.v2[0] + grocell.v3[0] ) * 0.5;
  move[1] = ( grocell.v1[1] + grocell.v2[1] + grocell.v3[1] ) * 0.5;
  move[2] = ( grocell.v1[2] + grocell.v2[2] + grocell.v3[2] ) * 0.5;

  int nAtom = nSize / 3;

  String strGro = c_strGro;
  vector<String> lines = strGro.Split("\n", false);
  for(int i = 0; i < nAtom; ++i) {
    String& line = lines[i+2];
    double x  = xyz[3*i+0] * M_TO_NM + move[0];
    double y  = xyz[3*i+1] * M_TO_NM + move[1];
    double z  = xyz[3*i+2] * M_TO_NM + move[2];
    double vx = vel[3*i+0] * M_TO_NM * PS_TO_S; // m/s -> nm/ps
    double vy = vel[3*i+1] * M_TO_NM * PS_TO_S;
    double vz = vel[3*i+2] * M_TO_NM * PS_TO_S;

    line = line.substr(0, 15) + 
      String::SPrintf("%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f",
                      x, y, z, vx, vy, vz);
  }
  
  vector<String> tokens = lines[nAtom+2].Split();
  GROCELL cell_old;
  if( 2 < tokens.size() ) {
    cell_old.v1[0] = tokens[0].ToDouble();
    cell_old.v1[1] = 0.0;
    cell_old.v1[2] = 0.0;
    cell_old.v2[0] = 0.0;
    cell_old.v2[1] = tokens[1].ToDouble();
    cell_old.v2[2] = 0.0;
    cell_old.v3[0] = 0.0;
    cell_old.v3[1] = 0.0;
    cell_old.v3[2] = tokens[2].ToDouble();
    if( 8 < tokens.size() ) {
      cell_old.v1[1] = tokens[3].ToDouble();
      cell_old.v1[2] = tokens[4].ToDouble();
      cell_old.v2[0] = tokens[5].ToDouble();
      cell_old.v2[2] = tokens[6].ToDouble();
      cell_old.v3[0] = tokens[7].ToDouble();
      cell_old.v3[1] = tokens[8].ToDouble();
    }
  }

  bool bChange = false;
  if( fabs(cell_old.v1[0] - grocell.v1[0]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v1[1] - grocell.v1[1]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v1[2] - grocell.v1[2]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v2[0] - grocell.v2[0]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v2[1] - grocell.v2[1]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v2[2] - grocell.v2[2]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v3[0] - grocell.v3[0]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v3[1] - grocell.v3[1]) > 1e-5 ) bChange = true;
  if( fabs(cell_old.v3[2] - grocell.v3[2]) > 1e-5 ) bChange = true;
  if( bChange ) {
    bool bCuboid = true;
    if( fabs(grocell.v1[1]) > 1e-5 ) bCuboid = false;
    if( fabs(grocell.v1[2]) > 1e-5 ) bCuboid = false;
    if( fabs(grocell.v2[0]) > 1e-5 ) bCuboid = false;
    if( fabs(grocell.v2[2]) > 1e-5 ) bCuboid = false;
    if( fabs(grocell.v3[0]) > 1e-5 ) bCuboid = false;
    if( fabs(grocell.v3[1]) > 1e-5 ) bCuboid = false;
    String line = String::SPrintf(" %9.5f %9.5f %9.5f",
                                  grocell.v1[0], grocell.v2[1], grocell.v3[2]);
    if( !bCuboid )  
      line += String::SPrintf(" %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f",
                              grocell.v1[1], grocell.v1[2], 
                              grocell.v2[0], grocell.v2[2],
                              grocell.v3[0], grocell.v3[1]);
    lines[nAtom+2] = line;
  }

  strGro = "";
  for(int i = 0; i < lines.size(); ++i) {
    strGro += lines[i] + "\n";
  }

  free(c_strGro);
  int n = strGro.length();
  c_strGro = (char*)malloc((n+1) * sizeof(char));
  strcpy(c_strGro, strGro.c_str());
}

extern "C"
void c_update_gro_names(const char* pcSessionName)
{
  string strSessionName = pcSessionName;
  GromacsData gmx_data(strSessionName);
  GromacsData gmx_name_data(strSessionName);

  String strGro = c_strGro;

  // Read GRO from string
  if( !gmx_data.ParseGroFromString(strGro) ) {
    return;
  }

  // Read GRO
  if( !gmx_name_data.ReadGRO() ) {
    return;
  }

  // copy names and residus
  int nAtom = gmx_data.m_groatoms.size();
  if( nAtom != gmx_name_data.m_groatoms.size() ) {
    return;
  }
  
  for( int i = 0; i < nAtom; ++i ) {
    gmx_data.m_groatoms[i].resname = gmx_name_data.m_groatoms[i].resname;
    gmx_data.m_groatoms[i].name    = gmx_name_data.m_groatoms[i].name;
  }

  strGro = gmx_data.WriteGROString();

  free(c_strGro);
  int n = strGro.length();
  c_strGro = (char*)malloc((n+1) * sizeof(char));
  strcpy(c_strGro, strGro.c_str());
}



