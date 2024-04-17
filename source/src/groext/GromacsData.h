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
  @brief  Utilities to deal with GROMACS data format
* @author  Tetsuro Nagai
*/
#ifndef GROMACS_DATA_H
#define GROMACS_DATA_H

#include "String.h"
#include "ModylasUtil.h"
#include <map>
#include <set>

using namespace std;

class GROATOM
{
public:
  GROATOM()
    : resid(0)
  {
    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }

public:
  int    resid;
  String resname;
  String name;
  double pos[3];
  double vel[3];
};

class GROCELL
{
public:
  GROCELL() {
    v1[0] = 0.0;
    v1[1] = 0.0;
    v1[2] = 0.0;
    v2[0] = 0.0;
    v2[1] = 0.0;
    v2[2] = 0.0;
    v3[0] = 0.0;
    v3[1] = 0.0;
    v3[2] = 0.0;
  }
public:
  double v1[3];
  double v2[3];
  double v3[3];
};

class TOPDEFAULTS
{
public:
  TOPDEFAULTS()
    : nbfunc(1), comb_rule(1), gen_pairs(true), fudgeLJ(1.0), fudgeQQ(1.0)
  {
  }

public:
  int    nbfunc;
  int    comb_rule;
  bool   gen_pairs;
  double fudgeLJ;
  double fudgeQQ;
};

class TOPATOMTYPE
{
public:
  TOPATOMTYPE()
    : mass(0.0), charge(0.0), sigma(0.0), epsilon(0.0)
  {
  }

public:
  String name;
  String btype;
  double mass;
  double charge;
  double sigma;
  double epsilon;
};

class TOPBONDTYPE
{
public:
  TOPBONDTYPE()
    : funct(1), r(0.0), k(0.0)
  {
  }

public:
  String namei;
  String namej;
  int funct;
  double r;
  double k;
};

class TOPANGLETYPE
{
public:
  TOPANGLETYPE()
    : funct(1), theta(0.0), k(0.0), ub0(0.0), cub(0.0)
  {
  }

public:
  String namei;
  String namej;
  String namek;
  int    funct;
  double theta;
  double k;
  double ub0;
  double cub;
};

class TOPDIHEDRALTYPE
{
public:
  TOPDIHEDRALTYPE()
    : funct(9), phase(0.0), kd(0.0), pn(1),
      c0(0.0), c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0)
  {
  }

public:
  String namei;
  String namej;
  String namek;
  String namel;
  int    funct;
  double phase;
  double kd;
  int    pn;
  // for Ryckaert-Bellemans
  double c0, c1, c2, c3, c4, c5;
};

class TOPCMAPTYPE
{
public:
  TOPCMAPTYPE()
  {
  }

public:
  String namei;
  String namej;
  String namek;
  String namel;
  String namem;
};

class TOPPAIRTYPE
{
public:
  TOPPAIRTYPE()
    : sigma(0.0), epsilon(0.0)
  {
  }

public:
  String namei;
  String namej;
  double sigma;
  double epsilon;
};

class TOPATOM
{
public:
  TOPATOM()
    : resid(0), charge(0.0), mass(0.0), cgnr(1),
      bHasCharge(false), bHasMass(false)
  {
  }

public:
  String type;
  int    resid;
  String resname;
  String name;
  int    cgnr; // charge group number
  bool   bHasCharge;
  double charge;
  bool   bHasMass;
  double mass;
};

class TOPBOND
{
public:
  TOPBOND()
    : ai(0), aj(0), funct(1), bHasR(false), r(0.0), k(0.0)
  {
  }
  TOPBOND(int _ai, int _aj, int _funct, bool _bHasR, double _r, double _k)
    : ai(_ai), aj(_aj), funct(_funct), bHasR(_bHasR), r(_r), k(_k)
  {
  }

public:
  int    ai;
  int    aj;
  int    funct;
  bool   bHasR;
  double r;
  double k;
};

class TOPPAIR
{
public:
  TOPPAIR()
    : ai(0), aj(0), bHasSigma(false), sigma(0.0), epsilon(0.0)
  {
  }

public:
  int    ai;
  int    aj;
  bool   bHasSigma;
  double sigma;
  double epsilon;
};

class TOPCMAP
{
public:
  TOPCMAP()
    : ai(0), aj(0), ak(0), al(0), am(0)
  {
  }

public:
  int ai;
  int aj;
  int ak;
  int al;
  int am;
};

class TOPANGLE
{
public:
  TOPANGLE()
    : ai(0), aj(0), ak(0), funct(1),
      bHasThetaK(false), theta(0.0), k(0.0),
      bHasUB(false), ub0(0.0), cub(0.0)
  {
  }

public:
  int    ai;
  int    aj;
  int    ak;
  int    funct;
  bool   bHasThetaK;
  double theta;
  double k;
  bool   bHasUB;
  double ub0;
  double cub;
};

class TOPDIHEDRAL
{
public:
  TOPDIHEDRAL()
    : ai(0), aj(0), ak(0), al(0), funct(9),
      bHasPhase(false), phase(0.0), kd(0.0), bHasN(false), pn(1),
      bHasC(false), c0(0.0), c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0)
  {
  }

public:
  int    ai;
  int    aj;
  int    ak;
  int    al;
  int    funct;
  bool   bHasPhase;
  double phase;
  double kd;
  bool   bHasN;
  int    pn;
  // for Ryckaert-Bellemans
  bool bHasC;
  double c0, c1, c2, c3, c4, c5;
};

class TOPSETTLE
{
public:
  TOPSETTLE()
    : iow(1), funct(1), doh(0.0), dhh(0.0)
  {
  }
public:
  int     iow;
  int     funct;
  double  doh;
  double  dhh;
};

class TOPEXCLUSION
{
public:
  TOPEXCLUSION()
    : ai(0)
  {
  }
public:
  int         ai;
  vector<int> vaj;
};

class TOPMOLECULETYPE
{
public:
  TOPMOLECULETYPE() : nrexcl(3) {}
  
public:
  String               name;
  int                  nrexcl;
  vector<TOPATOM>      atoms;
  vector<TOPBOND>      bonds;
  vector<TOPPAIR>      pairs;
  vector<TOPANGLE>     angles;
  vector<TOPCMAP>      cmaps;
  vector<TOPDIHEDRAL>  dihedrals;
  vector<TOPSETTLE>    settles;
  vector<TOPEXCLUSION> exclusions;
};

class TOPMOLECULE
{
public:
  TOPMOLECULE() : nmol(0) {}
  
public:
  String molname;
  int    nmol;
};

class GromacsData
{
public:
  GromacsData();
  GromacsData(const String& f) : m_filebase(f) {}
  GromacsData(const String& f, const String& includePath) : m_filebase(f) {
    m_includeDirs = includePath.Split(":");
  }

  bool ReadGRO();
  bool ReadTOP(String strIncKeys = "");

  bool WriteGRO() const;
  bool WriteTOP() const;
  String WriteGROString() const;
  String WriteTOPString() const;
  
  bool ParseGroFromString( const String& str );

private:
  bool GetTopString(const String& filename, map<String, String>& defined,
                    vector<String>& lines) const;
  bool SearchItpFile(String& itpfile, const String& topfile) const;
  bool IsWriteLine(const map<String, String>& defined,
		   const vector<pair<String, bool> >& in_ifdef) const;
  String ReplaceDefined(const String& line, 
                        const map<String, String>& defined) const;
  bool ReplaceTokenString(String& line, const String& from, const String& to) const;

public:
  vector<GROATOM> m_groatoms;
  GROCELL         m_grocell;

  TOPDEFAULTS             m_defaults;
  vector<TOPATOMTYPE>     m_atomtypes;
  vector<TOPPAIRTYPE>     m_pairtypes;
  vector<TOPBONDTYPE>     m_bondtypes;
  vector<TOPANGLETYPE>    m_angletypes;
  vector<TOPDIHEDRALTYPE> m_dihedraltypes;
  vector<TOPCMAPTYPE>     m_cmaptypes;
  vector<TOPMOLECULETYPE> m_moltypes;
  vector<TOPMOLECULE>     m_mols;
  
private:
  String m_filebase;
  vector<String> m_includeDirs;
};


#endif //GROMACS_DATA_H


