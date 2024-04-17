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
#include "GromacsData.h"

#include <fstream>
#include <iostream>
#include <stack>
#include <algorithm>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <iterator>
#include <sstream>
#include "cmap_itp.h"

using namespace std;

bool GromacsData::ReadGRO()
{
  String filename = m_filebase + ".gro";

  ifstream ifs( filename.c_str() );
  if( ifs.fail() ) {
    cerr << "Info) Could not open " + filename + ". Skip reading." << endl;
    return false;
  }
  //cout << "Info) Reading " + filename + " ..." << endl;

  String strall( string(
    ( istreambuf_iterator<char>(ifs) ),
    ( istreambuf_iterator<char>()    ) ) );

  bool b = ParseGroFromString( strall );
  if( !b ) {
    cerr << "Error) Could not parse gro file contents." << endl;
    return false;
  }

  return true;
}

bool GromacsData::ParseGroFromString( const String& str )
{
  stringstream ss(str);

  String line;
  try {
    getline(ss, line);
    getline(ss, line);
    int natoms = line.ToInt();
    m_groatoms.resize( natoms );
    for( int i = 0; i < natoms; ++i ) {
      getline(ss, line);
      GROATOM& atom = m_groatoms[i];
      atom.resid   = line.substr(0 ,5).ToInt();
      atom.resname = line.substr(5 ,5).Trim();
      atom.name    = line.substr(10,5).Trim();
      atom.pos[0]  = line.substr(20,8).ToDouble();
      atom.pos[1]  = line.substr(28,8).ToDouble();
      atom.pos[2]  = line.substr(36,8).ToDouble();
      if( 68 <= line.length() ) {
        atom.vel[0]  = line.substr(44,8).ToDouble();
        atom.vel[1]  = line.substr(52,8).ToDouble();
        atom.vel[2]  = line.substr(60,8).ToDouble();
      }
    }
    getline(ss, line);
    vector<String> tokens = line.Split();
    if( tokens.size() < 3 ) throw "";
    m_grocell.v1[0] = tokens[0].ToDouble();
    m_grocell.v1[1] = 0.0;
    m_grocell.v1[2] = 0.0;
    m_grocell.v2[0] = 0.0;
    m_grocell.v2[1] = tokens[1].ToDouble();
    m_grocell.v2[2] = 0.0;
    m_grocell.v3[0] = 0.0;
    m_grocell.v3[1] = 0.0;
    m_grocell.v3[2] = tokens[2].ToDouble();
    if( 8 < tokens.size() ) {
      m_grocell.v1[1] = tokens[3].ToDouble();
      m_grocell.v1[2] = tokens[4].ToDouble();
      m_grocell.v2[0] = tokens[5].ToDouble();
      m_grocell.v2[2] = tokens[6].ToDouble();
      m_grocell.v3[0] = tokens[7].ToDouble();
      m_grocell.v3[1] = tokens[8].ToDouble();
    }
  }
  catch(...) {
    cerr << "ERROR) Failed to read gro file" << endl;
    cerr << "ERROR) Error while reading line : " << line << endl;
    return false;
  }

  return true;
}

bool GromacsData::ReadTOP(String strIncKeys)
{
  String filename = m_filebase + ".top";

  vector<String> lines;
  map<String, String> defined;
  vector<String> tokens = strIncKeys.Split();
  for( unsigned int i = 0; i < tokens.size(); ++i ) {
    defined[tokens[i]] = "";
  }
  if( !GetTopString(filename, defined, lines) ) {
    cerr << "Info) Could not read " + filename + ". Skip reading." << endl;
    return false;
  }
  //cout << "Info) Reading " + filename + " ..." << endl;

  String line;
  try {
    String mode;

    for(unsigned int iLine = 0; iLine < lines.size(); ++iLine) {
      line = lines[iLine];

      string::size_type i = line.find(";");
      if( i != NPOS ) {
        line = line.substr(0, i);
      }
      i = line.find("#");
      if( i != NPOS ) {
        line = line.substr(0, i);
      }
      if( 0 < line.length() && line.substr(0,1) == "*" ) {
        continue;
      }
      line = line.Trim();
      if( line == "" ) continue;
      if( line[0] == '[' ) {
        if( line.back() != ']' ) throw "";
        mode = line.substr(1, line.size() - 2).Trim();
        continue;
      }
      vector<String> tokens = line.Split();
      int ntokens = tokens.size();
      if( mode == "defaults" ) {
        if( ntokens < 5 ) throw "";
        m_defaults.nbfunc    = tokens[0].ToInt();
        m_defaults.comb_rule = tokens[1].ToInt();
        m_defaults.gen_pairs = (tokens[2] == "yes");
        m_defaults.fudgeLJ   = tokens[3].ToDouble();
        m_defaults.fudgeQQ   = tokens[4].ToDouble();
        continue;
      }
      if( mode == "pairtypes" ) {
        if( ntokens < 5 ) throw "";
        TOPPAIRTYPE ptype;
        ptype.namei   = tokens[0];
        ptype.namej   = tokens[1];
        ptype.sigma   = tokens[3].ToDouble();
        ptype.epsilon = tokens[4].ToDouble();
        m_pairtypes.push_back(ptype);
        continue;
      }
      if( mode == "atomtypes" ) {
        if( ntokens < 7 ) throw "";
        // atomtypes might include bonded type and/or atomic number
        // see push_at() function in toppush.c of GROMACS to find how to identify the format
        TOPATOMTYPE atype;
        if( tokens[3].length() == 1 && isalpha(tokens[3][0]) ) {
          atype.name    = tokens[0];
          atype.btype   = tokens[0];
          atype.mass    = tokens[1].ToDouble();
          atype.charge  = tokens[2].ToDouble();
          atype.sigma   = tokens[4].ToDouble();
          atype.epsilon = tokens[5].ToDouble();
        }
        else if( tokens[5].length() == 1 && isalpha(tokens[5][0]) ) {
          atype.name    = tokens[0];
          atype.btype   = tokens[1];
          atype.mass    = tokens[3].ToDouble();
          atype.charge  = tokens[4].ToDouble();
          atype.sigma   = tokens[6].ToDouble();
          atype.epsilon = tokens[7].ToDouble();
        }
        else if( isalpha(tokens[1][0]) ) {
          atype.name    = tokens[0];
          atype.btype   = tokens[1];
          atype.mass    = tokens[2].ToDouble();
          atype.charge  = tokens[3].ToDouble();
          atype.sigma   = tokens[5].ToDouble();
          atype.epsilon = tokens[6].ToDouble();
        }
        else {
          atype.name    = tokens[0];
          atype.btype   = tokens[0];
          atype.mass    = tokens[2].ToDouble();
          atype.charge  = tokens[3].ToDouble();
          atype.sigma   = tokens[5].ToDouble();
          atype.epsilon = tokens[6].ToDouble();
        }
        m_atomtypes.push_back(atype);
        continue;
      }
      if( mode == "bondtypes" ) {
        if( ntokens < 5 ) throw "";
        TOPBONDTYPE btype;
        btype.namei = tokens[0];
        btype.namej = tokens[1];
        btype.funct = tokens[2].ToInt();
        btype.r     = tokens[3].ToDouble();
        btype.k     = tokens[4].ToDouble();
        m_bondtypes.push_back(btype);
        continue;
      }
      if( mode == "angletypes" ) {
        if( ntokens < 6 ) throw "";
        TOPANGLETYPE angtype;
        angtype.namei = tokens[0];
        angtype.namej = tokens[1];
        angtype.namek = tokens[2];
        angtype.funct = tokens[3].ToInt();
        angtype.theta = tokens[4].ToDouble();
        angtype.k     = tokens[5].ToDouble();
        if( angtype.funct == 5 ) {
          if( ntokens < 8 ) throw "";
          angtype.ub0 = tokens[6].ToDouble();
          angtype.cub = tokens[7].ToDouble();
        }
        m_angletypes.push_back(angtype);
        continue;
      }
      if( mode == "dihedraltypes" ) {
        if( ntokens < 4 ) throw "";
        bool b4atoms = false;
        try {
          tokens[2].ToInt();
        }
        catch(...) {
          b4atoms = true;
        }
        if( !b4atoms ) {
          tokens.insert(tokens.begin()    , "X");
          tokens.insert(tokens.begin() + 3, "X");
          ntokens = tokens.size();
        }
        if( ntokens < 7 ) throw "";
        TOPDIHEDRALTYPE dihedtype;
        dihedtype.namei = tokens[0];
        dihedtype.namej = tokens[1];
        dihedtype.namek = tokens[2];
        dihedtype.namel = tokens[3];
        dihedtype.funct = tokens[4].ToInt();
        if( dihedtype.funct == 3 ) {
          if( ntokens < 11 ) throw "";
          dihedtype.c0 = tokens[5].ToDouble();
          dihedtype.c1 = tokens[6].ToDouble();
          dihedtype.c2 = tokens[7].ToDouble();
          dihedtype.c3 = tokens[8].ToDouble();
          dihedtype.c4 = tokens[9].ToDouble();
          dihedtype.c5 = tokens[10].ToDouble();
        }
        else if( dihedtype.funct == 2 ) {
          if( ntokens < 7 ) throw "";
          dihedtype.phase = tokens[5].ToDouble();
          dihedtype.kd    = tokens[6].ToDouble();
        }
        else if( dihedtype.funct == 1 || dihedtype.funct == 4 || dihedtype.funct == 9 ) {
          if( ntokens < 8 ) throw "";
          dihedtype.phase = tokens[5].ToDouble();
          dihedtype.kd    = tokens[6].ToDouble();
          dihedtype.pn    = tokens[7].ToInt();
        }
        m_dihedraltypes.push_back(dihedtype);
        continue;
      }
      if( mode == "cmaptypes" ) {
        if( tokens[0].substr(0, 1) == "-" || isdigit(tokens[0][0])) continue;
        if( ntokens < 5 ) throw "";
        TOPCMAPTYPE ctype;
        ctype.namei = tokens[0];
        ctype.namej = tokens[1];
        ctype.namek = tokens[2];
        ctype.namel = tokens[3];
        ctype.namem = tokens[4];
        m_cmaptypes.push_back(ctype);
        continue;
      }
      if( mode == "moleculetype" ) {
        if( ntokens < 1 ) throw "";
        m_moltypes.push_back(TOPMOLECULETYPE());
        m_moltypes.back().name = tokens[0];
        if( 1 < ntokens ) {
          m_moltypes.back().nrexcl = tokens[1].ToInt();
        }
      }
      if( mode == "atoms" ) {
        if( ntokens < 6 ) throw "";
        TOPATOM atom;
        atom.type    = tokens[1];
        atom.resid   = tokens[2].ToInt();
        atom.resname = tokens[3];
        atom.name    = tokens[4];
        atom.cgnr    = tokens[5].ToInt();
        if( 6 < ntokens ) {
          atom.charge     = tokens[6].ToDouble();
          atom.bHasCharge = true;
        }
        if( 7 < ntokens ) {
          atom.mass       = tokens[7].ToDouble();
          atom.bHasMass   = true;
        }
        m_moltypes.back().atoms.push_back(atom);
        continue;
      }
      if( mode == "bonds" ) {
        if( ntokens < 3 ) throw "";
        TOPBOND bond;
        bond.ai    = tokens[0].ToInt();
        bond.aj    = tokens[1].ToInt();
        bond.funct = tokens[2].ToInt();
        if( 4 < ntokens ) {
          bond.bHasR = true;
          bond.r     = tokens[3].ToDouble();
          bond.k     = tokens[4].ToDouble();
        }
        m_moltypes.back().bonds.push_back(bond);
        continue;
      }
      if( mode == "pairs" ) {
        if( ntokens < 2 ) throw "";
        TOPPAIR pair;
        pair.ai = tokens[0].ToInt();
        pair.aj = tokens[1].ToInt();
        if( 4 < ntokens ) {
          pair.bHasSigma = true;
          pair.sigma   = tokens[3].ToDouble();
          pair.epsilon = tokens[4].ToDouble();
        }
        m_moltypes.back().pairs.push_back(pair);
        continue;
      }
      if( mode == "angles" ) {
        if( ntokens < 4 ) throw "";
        TOPANGLE angle;
        angle.ai    = tokens[0].ToInt();
        angle.aj    = tokens[1].ToInt();
        angle.ak    = tokens[2].ToInt();
        angle.funct = tokens[3].ToInt();
        if( 5 < ntokens ) {
          angle.bHasThetaK = true;
          angle.theta = tokens[4].ToDouble();
          angle.k     = tokens[5].ToDouble();
        }
        if( angle.funct == 5 && 7 < ntokens ) {
          angle.bHasUB = true;
          angle.ub0    = tokens[6].ToDouble();
          angle.cub    = tokens[7].ToDouble();
        }
        m_moltypes.back().angles.push_back(angle);
        continue;
      }
      if( mode == "cmap" ) {
        if( ntokens < 5 ) throw "";
        TOPCMAP cmap;
        cmap.ai = tokens[0].ToInt();
        cmap.aj = tokens[1].ToInt();
        cmap.ak = tokens[2].ToInt();
        cmap.al = tokens[3].ToInt();
        cmap.am = tokens[4].ToInt();
        m_moltypes.back().cmaps.push_back(cmap);
        continue;
      }
      if( mode == "dihedrals" ) {
        if( ntokens < 5 ) throw "";
        TOPDIHEDRAL dihedral;
        dihedral.ai    = tokens[0].ToInt();
        dihedral.aj    = tokens[1].ToInt();
        dihedral.ak    = tokens[2].ToInt();
        dihedral.al    = tokens[3].ToInt();
        dihedral.funct = tokens[4].ToInt();
        if( dihedral.funct == 3 ) {
          if( 10 < ntokens ) { 
            dihedral.bHasC = true;
            dihedral.c0 = tokens[5].ToDouble();
            dihedral.c1 = tokens[6].ToDouble();
            dihedral.c2 = tokens[7].ToDouble();
            dihedral.c3 = tokens[8].ToDouble();
            dihedral.c4 = tokens[9].ToDouble();
            dihedral.c5 = tokens[10].ToDouble();
          }
        }
        else {
          if( 6 < ntokens ) {
            dihedral.bHasPhase = true;
            dihedral.phase = tokens[5].ToDouble();
            dihedral.kd    = tokens[6].ToDouble();
          }
          if( 7 < ntokens ) {
            dihedral.bHasN = true;
            dihedral.pn    = tokens[7].ToInt();
          }
        }
        m_moltypes.back().dihedrals.push_back(dihedral);
        continue;
      }
      if( mode == "settles" ) {
        if( ntokens < 4 ) throw "";
        TOPSETTLE settle;
        settle.iow   = tokens[0].ToInt();
        settle.funct = tokens[1].ToInt();
        settle.doh   = tokens[2].ToDouble();
        settle.dhh   = tokens[3].ToDouble();
        m_moltypes.back().settles.push_back(settle);
        continue;
      }
      if( mode == "exclusions" ) {
        if( ntokens < 2 ) throw "";
        TOPEXCLUSION exclusion;
        exclusion.ai = tokens[0].ToInt();
        for( int j = 1; j < tokens.size(); ++j ) {
          exclusion.vaj.push_back(tokens[j].ToInt());
        }
        m_moltypes.back().exclusions.push_back(exclusion);
        continue;
      }
      if( mode == "molecules" ) {
        if( ntokens < 2 ) throw "";
        TOPMOLECULE mol;
        mol.molname = tokens[0];
        mol.nmol    = tokens[1].ToInt();
        m_mols.push_back(mol);
        continue;
      }
    }
  }
  catch(...) {
    cerr << "Error while reading the line : " << line << endl;
    cerr << "ERROR) Failed to read top file" << endl;
    return false;
  }

  return true;
}

String GromacsData::WriteGROString() const
{
  int natoms = m_groatoms.size();

  String str;
  str += "gro file written by MODYLAS\n";
  str += String(natoms) + "\n";

  for( int i = 0; i < natoms; ++i ) {
    const GROATOM& atom = m_groatoms[i];
    int    resid   = atom.resid % 100000;
    String resname = atom.resname.substr(0,5);
    String name    = atom.name.substr(0,5);
    int    index   = ( i + 1 ) % 100000;
    double x       = atom.pos[0];
    double y       = atom.pos[1];
    double z       = atom.pos[2];
    double vx      = atom.vel[0];
    double vy      = atom.vel[1];
    double vz      = atom.vel[2];
    str += String::SPrintf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                           resid, resname.c_str(), name.c_str(), index,
                           x, y, z, vx, vy, vz);
  }

  double ax = m_grocell.v1[0];
  double ay = m_grocell.v1[1];
  double az = m_grocell.v1[2];
  double bx = m_grocell.v2[0];
  double by = m_grocell.v2[1];
  double bz = m_grocell.v2[2];
  double cx = m_grocell.v3[0];
  double cy = m_grocell.v3[1];
  double cz = m_grocell.v3[2];
  str += String(ax) + " " + String(by) + " " + String(cz);
  if( 1e-10 < fabs(bx) || 1e-10 < fabs(cx) || 1e-10 < fabs(cy) ) {
    str +=
      " " + String(ay) + " " + String(az) +
      " " + String(bx) + " " + String(bz) +
      " " + String(cx) + " " + String(cy);
  }
  str += "\n";

  return str;
}

bool GromacsData::WriteGRO() const
{
  String str = WriteGROString();

  String filename = m_filebase + ".gro";

  ofstream ofs( filename.c_str() );
  if( ofs.fail() ) {
    cerr << "ERROR) Failed to open the file for writing : " + filename << endl;
    return false;
  }
  ofs << str;
  ofs.close();

  return true;
}

String GromacsData::WriteTOPString() const
{
  String str;

  str += "[ defaults ]\n";
  str += "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n";
  str += "1 2 yes 1.0 1.0\n";
  str += "\n";

  // atomtypes
  if( !m_atomtypes.empty() ) {
    str += "[ atomtypes ]\n";
    str += ";name at.num  mass  charge  ptype sigma epsilon\n";
    for( unsigned int i = 0; i < m_atomtypes.size(); ++i ) {
      const TOPATOMTYPE& atype = m_atomtypes[i];
      String mass   = String(atype.mass   );
      String charge = String(atype.charge );
      String sigma  = String(atype.sigma  );
      String eps    = String(atype.epsilon);
      str += atype.name + " " + atype.btype + " " + mass + " " + charge
        + " A " + sigma + " " + eps + "\n";
    }
    str += "\n";
  }
  
  if( !m_pairtypes.empty() ) {
    str += "[ pairtypes ]\n";
    str += "; i j func sigma1-4 epsilon1-4\n";
    for( unsigned int i = 0; i < m_pairtypes.size(); ++i ) {
      const TOPPAIRTYPE& ptype = m_pairtypes[i];
      String sigma = String(ptype.sigma  );
      String eps   = String(ptype.epsilon);
      str += ptype.namei + " " + ptype.namej + " 1 "
        + sigma + " " + eps + "\n";
    }
    str += "\n";
  }

  if( !m_bondtypes.empty() ) {
    str += "[ bondtypes ]\n";
    str += "; i j func b0 kb\n";
    for( unsigned int i = 0; i < m_bondtypes.size(); ++i ) {
      const TOPBONDTYPE& btype = m_bondtypes[i];
      String f  = String(btype.funct);
      String b0 = String(btype.r    );
      String kb = String(btype.k    );
      str += btype.namei + " " + btype.namej + " " + f + " " + b0 + " " + kb + "\n";
    }
    str += "\n";
  }

  if( !m_angletypes.empty() ) {
    str += "[ angletypes ]\n";
    str += "; i j k func th0 cth ub0 cub\n";
    for( unsigned int i = 0; i < m_angletypes.size(); ++i ) {
      const TOPANGLETYPE& angtype = m_angletypes[i];
      String f   = String(angtype.funct);
      String th0 = String(angtype.theta);
      String cth = String(angtype.k    );
      str += angtype.namei + " " + angtype.namej + " " + angtype.namek
        + " " + f + th0 + " " + cth;
      if( angtype.funct == 5 ) { 
        String ub0 = String(angtype.ub0);
        String cub = String(angtype.cub);
        str += " " + ub0 + " " + cub + "\n";
      }
      str += "\n";
    }
    str += "\n";
  }

  if( !m_dihedraltypes.empty() ) {
    str += "[ dihedraltypes ]\n";
    str += "; i j k l func phi0 cp mult\n";
    for( unsigned int i = 0; i < m_dihedraltypes.size(); ++i ) {
      const TOPDIHEDRALTYPE& dtype = m_dihedraltypes[i];
      String f     = String(dtype.funct);
      String phi0  = String(dtype.phase);
      String cp    = String(dtype.kd   );
      str += dtype.namei + " " + dtype.namej + " " + dtype.namek + " " + dtype.namel
        + " " + f + " " + phi0 + " " + cp;
      if( dtype.funct == 1 || dtype.funct == 9 ) {
        String multi = String(dtype.pn);
        str += " " + multi;
      }
      str += "\n";
    }
    str += "\n";
  }

  bool bWriteCmapTypes = false;
  for( unsigned int i = 0; i < m_moltypes.size(); ++i ) {
    if( !m_moltypes[i].cmaps.empty() ) {
      bWriteCmapTypes = true;
      break;
    }
  }
  if( bWriteCmapTypes ) {
    str += CMAP_ITP;
  }

  
  for( unsigned int imol = 0; imol < m_moltypes.size(); ++imol ) {
    const TOPMOLECULETYPE& moltype = m_moltypes[imol];
    str += "[ moleculetype ]\n";
    str += "; name  nrexcl\n";
    str += moltype.name + " " + String(moltype.nrexcl) +  "\n";
    str += "\n";

    if( !moltype.atoms.empty() ) {
      str += "[ atoms ]\n";
      str += "; nr type resi res atom cgnr charge mass\n";
      for( unsigned int i = 0; i < moltype.atoms.size(); ++i ) {
        const TOPATOM& atom = moltype.atoms[i];
        String nr    = String(int(i + 1));
        String resid = String(atom.resid);
        String cgnr  = String(atom.cgnr );
        str += nr + " " + atom.type + " " + resid + " " + atom.resname + " "
          + atom.name + " " + cgnr;
        if( atom.bHasCharge ) {
          str += " " + String(atom.charge);
          if( atom.bHasMass ) {
            str += " " + String(atom.mass);
          }
        }
        str += "\n";
      }
      str += "\n";
    }
    

    if( !moltype.bonds.empty() ) {
      str += "[ bonds ]\n";
      str += "; ai aj funct r k\n";
      for( unsigned int i = 0; i < moltype.bonds.size(); ++i ) {
        const TOPBOND& bond = moltype.bonds[i];
        String ai = String(bond.ai);
        String aj = String(bond.aj);
        String f  = String(bond.funct);
        str += ai + " " + aj + " " + f;
        if( bond.bHasR ) {
          String r  = String(bond.r);
          String k  = String(bond.k);
          str += " " + r + " " + k;
        }
        str += "\n";
      }
        str += "\n";
    }

    if( !moltype.pairs.empty() ) {
      str += "[ pairs ]\n";
      str += "; ai aj funct r e\n";
      for( unsigned int i = 0; i < moltype.pairs.size(); ++i ) {
        const TOPPAIR& pair = moltype.pairs[i];
        String ai = String(pair.ai);
        String aj = String(pair.aj);
        str += ai + " " + aj + " 1";
        if( pair.bHasSigma ) {
          String sigma(pair.sigma);
          String eps  (pair.epsilon);
          str += " " + sigma + " " + eps;
        }
        str += "\n";
      }
      str += "\n";
    }

    if( !moltype.angles.empty() ) {
      str += "[ angles ]\n";
      str += "; ai aj ak funct theta k\n";
      for( unsigned int i = 0; i < moltype.angles.size(); ++i ) {
        const TOPANGLE& angle = moltype.angles[i];
        String ai = String(angle.ai);
        String aj = String(angle.aj);
        String ak = String(angle.ak);
        String f  = String(angle.funct);
        str += ai + " " + aj + " " + ak + " " + f;
        if( angle.bHasThetaK ) {
          String theta = String(angle.theta);
          String k     = String(angle.k    );
          str += " " + theta + " " + k;
        }
        if( angle.bHasUB ) {
          String ub0 = String(angle.ub0);
          String cub = String(angle.cub);
          str += " " + ub0 + " " + cub;
        }
        str += "\n";
      }
      str += "\n";
    }

    if( !moltype.dihedrals.empty() ) {
      str += "[ dihedrals ]\n";
      str += "; ai aj ak al funct gamma k n\n";
      for( unsigned int i = 0; i < moltype.dihedrals.size(); ++i ) {
        const TOPDIHEDRAL& dihedral = moltype.dihedrals[i];
        String ai = String(dihedral.ai);
        String aj = String(dihedral.aj);
        String ak = String(dihedral.ak);
        String al = String(dihedral.al);
        String f  = String(dihedral.funct);
        str += ai + " " + aj + " " + ak + " " + al + " " + f;
        if( dihedral.bHasPhase ) {
          String gamma = String(dihedral.phase);
          String k     = String(dihedral.kd   );
          str += " " + gamma + " " + k;
        }
        if( dihedral.bHasN ) {
          String n = String(dihedral.pn   );
          str += " " + n;
        }
        str += "\n";
      }
      str += "\n";
    }

    if( !moltype.cmaps.empty() ) {
      str += "[ cmap ]\n";
      str += "; ai aj ak al funct gamma k n\n";
      for( unsigned int i = 0; i < moltype.cmaps.size(); ++i ) {
        const TOPCMAP& cmap = moltype.cmaps[i];
        String ai = String(cmap.ai);
        String aj = String(cmap.aj);
        String ak = String(cmap.ak);
        String al = String(cmap.al);
        String am = String(cmap.am);
        str += ai + " " + aj + " " + ak + " " + al + " " + am + " 1\n";
      }
      str += "\n";
    }
  } // species

  str += "[system]\n";
  str += "converted by modylas2gmx\n";
  str += "\n";
  
  str += "[molecules]\n";
  for( unsigned int i = 0; i < m_mols.size(); ++i ) {
    String nmol = String(m_mols[i].nmol);
    str +=  m_mols[i].molname + " " + nmol + "\n";
  }
  str += "\n";
  
  return str;
}

bool GromacsData::WriteTOP() const
{
  String filename = m_filebase + ".top";

  String str = WriteTOPString();

  ofstream ofs( filename.c_str() );
  if( ofs.fail() ) {
    cerr << "ERROR) Failed to open the file for writing : " + filename << endl;
    return false;
  }
  ofs << str;
  ofs.close();

  return true;
}

bool GromacsData::GetTopString(const String& filename, map<String, String>& defined,
                               vector<String>& lines) const
{
  ifstream ifs( filename.c_str() );
  if( ifs.fail() ) {
    cerr << "ERROR) Failed to open the file for reading : " + filename << endl;
    return false;
  }

  String line_org;
  vector<pair<String, bool> > in_ifdef;
  try {
    while( getline(ifs, line_org) ) {
      String line = line_org;
      // Remove comment
      string::size_type i = line.find(";");
      if( i != NPOS ) {
        line = line.substr(0, i);
      }
      line = line.Trim();

      vector<String> tokens = line.Split();
      if( 0 < tokens.size() ) {
        if( tokens[0] == "#else" ) {
          in_ifdef.back().second = !in_ifdef.back().second;
          continue;
        }
        if( tokens[0] == "#endif" ) {
          in_ifdef.pop_back();
          continue;
        }
        if( 1 < tokens.size() ) {
          if( tokens[0] ==  "#define" ) {
            String str = line.substr(tokens[0].length()).Trim();
            str = str.substr(tokens[1].length()).Trim();
            defined[tokens[1]] = str;
            continue;
          }
          if( tokens[0] ==  "#undef" ) {
            defined.erase(tokens[1]);
            continue;
          }
          if( tokens[0] == "#ifdef" ) {
            in_ifdef.push_back(pair<String, bool>(tokens[1], true));
            continue;
          }
          if( tokens[0] == "#ifndef" ) {
            in_ifdef.push_back(pair<String, bool>(tokens[1], false));
            continue;
          }
          if( tokens[0] == "#include" ) {
            if( !IsWriteLine(defined, in_ifdef) ) continue;
            
            String itpfile = line.substr(9).Trim(); // file name might contain space
            int n = itpfile.length();
            if( 1 < n && itpfile.substr(0, 1) == "\"" && itpfile.substr(n-1, 1) == "\"")
              itpfile = itpfile.substr(1, n-2); // remove double quatation
            if( !SearchItpFile(itpfile, filename) ) {
              cerr << "ERROR) itp file not found : " + itpfile << endl;
              return false;
            }
            
            vector<String> lines_itp;
            bool b = GetTopString(itpfile, defined, lines_itp);
            if( !b ) return false;
            copy(lines_itp.begin(), lines_itp.end(), back_inserter(lines));
            
            continue;
          }
        }
        
        if( IsWriteLine(defined, in_ifdef) ) {
          line = ReplaceDefined(line, defined);
          lines.push_back(line);
        }
      }
    }
  }
  catch(...) {
    cerr << "ERROR) Failed to read top/itp file : " << filename << endl;
    cerr << "ERROR) Error while reading line : " << line_org << endl;
    ifs.close();
    return false;
  }
  ifs.close();

  return true;
}

String GromacsData::ReplaceDefined(const String& line, const map<String, String>& defined) const
{
  set<String> used; //to avoid circular reference

  String line0;
  String line1 = line;
  while(line0 != line1) {
    line0 = line1;
    map<String, String>::const_iterator it = defined.begin();
    for(; it != defined.end(); ++it) {
      if( used.find(it->first) != used.end() ) continue;
      bool b = ReplaceTokenString(line1, it->first, it->second);
      if( b ) {
        used.insert(it->first);
      }
    }
  }
  
  return line0;
}

bool GromacsData::ReplaceTokenString(String& line, const String& from, const String& to) const
{
  bool bReplace = false;

  vector<String> tokens = line.Split();
  for(unsigned int i = 0; i < tokens.size(); ++i) {
    if( tokens[i] == from ) {
      tokens[i] = to;
      bReplace = true;
    }
  }

  if( bReplace ) {
    line = "";
    for(unsigned int i = 0; i < tokens.size(); ++i) {
      line += tokens[i] + " ";
    }
    line = line.Trim();
  }

  return bReplace;
}

// Search itp file from
//  - top file directory
//  - current directory
//  - specfied include path
bool GromacsData::SearchItpFile(String& itpfile, const String& topfile) const
{
  String topdir = topfile.GetDirectory();
  
  {
    String itpfullpath = topdir + "/" + itpfile;
    ifstream ifs( itpfullpath.c_str() );
    if( ifs.is_open() ) {
      itpfile = itpfullpath;
      return true;
    }
  }
  
  {
    ifstream ifs( itpfile.c_str() );
    if( ifs.is_open() ) return true;
  }

  for(unsigned int i = 0; i < m_includeDirs.size(); ++i) {
    String dir = m_includeDirs[i];
    String s = dir.substr(dir.length()-1, 1);
    if( (s != "/") && (s != "\\") ) {
      dir += "/";
    }
    String itpfullpath = dir + itpfile;
    ifstream ifs( itpfullpath.c_str() );
    if( ifs.is_open() ) {
      itpfile = itpfullpath;
      return true;
    }
  }
  
  return false;
}

bool GromacsData::IsWriteLine(const map<String, String>& defined,
			      const vector<pair<String, bool> >& in_ifdef) const
{
  for(unsigned int i = 0; i < in_ifdef.size(); ++i) {
    bool bDefined = defined.find(in_ifdef[i].first) != defined.end();
    if( bDefined != in_ifdef[i].second ) {
      return false;
    }
  }
  return true;
}

