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
#include "String.h"

#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdarg.h>
#include <algorithm>
#include <fstream>
#include <iterator>

using namespace std;

String::String( int n )
{
  stringstream ss;
  ss << n;
  *this = ss.str();
}

String::String( double d )
{
  stringstream ss;
  ss << scientific;
  ss << setprecision(6);
  ss << setw(14);
  ss << uppercase;
  ss << right;
  ss << d;
  *this = ss.str();
}

String String::RemoveExt() const
{
  String filename = *this;
  SIZE_TYPE pos = filename.find_last_of( "." );
  if( pos == NPOS ) {
    return filename;
  }
  return filename.substr(0, pos );
}

String String::GetExt() const
{
  String filename = *this;
  SIZE_TYPE pos = filename.find_last_of( "." );
  if( pos == NPOS ) {
    return "";
  }
  return filename.substr( pos + 1 );
}

String String::RemoveDir() const
{
  String filename = *this;
  SIZE_TYPE pos = filename.find_last_of( "/\\" );
  if( pos == NPOS ) {
    return filename;
  }
  return filename.substr( pos + 1 );
}

String String::GetDirectory() const
{
  String filename = *this;
  SIZE_TYPE pos = filename.find_last_of( "/\\" );
  if( pos == NPOS ) {
    return filename;
  }
  return filename.substr( 0, pos );
}


vector<String> String::Split( const String& delimiter,
  bool ignore_blank ) const
{
  vector<String> vstr;
  SIZE_TYPE start = 0;
  SIZE_TYPE pos( find_first_of( delimiter ) );
  while( pos != NPOS ) {
    vstr.push_back( substr(start, pos-start) );
    start = pos + 1;
    pos = find_first_of( delimiter, start );
  }
  vstr.push_back( substr(start, length()-start) );
    
  if( ! ignore_blank ) {
    return vstr;
  }

  vector<String> vstr_;
  for( unsigned int i = 0; i < vstr.size(); ++i ) {
    if( vstr[i] != "" ) {
      vstr_.push_back( vstr[i] );
    }
  }
  return vstr_;
}


double String::ToDouble() const
{
  double d;
  stringstream ss( *this );
  ss >> d;
  if( ss.fail() ) throw "Failed to convert string to double";
  return d;
}

int String::ToInt() const
{
  int n;
  stringstream ss( *this );
  ss >> n;
  if( ss.fail() ) throw "Failed to convert string to int";
  return n;
}

String String::ToUpper() const
{
  String ret = *this;
  transform( ret.begin(), ret.end(), ret.begin(), ::toupper );
  return ret;
}

String String::ToLower() const
{
  String ret = *this;
  transform( ret.begin(), ret.end(), ret.begin(), ::tolower );
  return ret;
}

String String::Trim() const
{
  SIZE_TYPE b = find_first_not_of(" \r\n\t");
  if( b == NPOS ) {
    return "";
  }
  SIZE_TYPE e = find_last_not_of(" \r\n\t");
  return substr( b, e - b + 1 );
}

String String::ReplaceAll( const String& key1,
  const String& key2 ) const
{
  String ret = *this;
  SIZE_TYPE pos( ret.find( key1 ) );
  while( pos != NPOS ) {
    ret.replace( pos, key1.length(), key2 );
    SIZE_TYPE start = pos + key2.length();
    pos = ret.find( key1, start );
  }
  return ret;
}

String String::SPrintf( const char* format, ... )
{
  int nBuf = 256;
  vector<char> buf(nBuf);
  while(1) {
    va_list argp;
    va_start(argp, format);
    int nSize = vsnprintf(&(buf[0]), nBuf, format, argp);
    if( 0 <= nSize && nSize < nBuf ) {
      break;
    }
    nBuf *= 2;
    buf.resize( nBuf );
    va_end(argp);
  }
  return String( &(buf[0]) );
}

String String::ReadFileString( const String& filename )
{
  ifstream ifs( filename.c_str() );
  if( ifs.fail() ) {
    return String();
  }
  istreambuf_iterator<char> it(ifs);
  istreambuf_iterator<char> last;
  string str(it, last);
  return str;
}

