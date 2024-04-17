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
  @brief  Utilities to deal with string
* @author  Tetsuro Nagai
*/
#ifndef MODYLAS_STRING_H
#define MODYLAS_STRING_H

#include <string>
#include <vector>

using namespace std;

typedef string::size_type SIZE_TYPE;
const SIZE_TYPE NPOS = string::npos;

// derived class from std::string with useful functions 
class String : public string
{
public:
  String() {}
  String( const string& str ) : string(str) {}
  String( const char* p ) : string(p) {}
  String( int n );
  String( double d );
  ~String() {}

  String& operator=( const char* p )
  {
      *this = String( p );
      return *this;
  }
  String substr( SIZE_TYPE index, SIZE_TYPE length = NPOS ) const
  {
    return String( string::substr( index, length ) );
  }
  
  char& back() { return *rbegin(); };

  String RemoveExt() const;
  String GetExt() const;
  String RemoveDir() const;
  String GetDirectory() const;
  vector<String> Split( const String& delimiter = " \n\r\t",
    bool ignore_blank = true ) const;
  double ToDouble() const ;
  int ToInt() const;
  String ToUpper() const;
  String ToLower() const;
  String Trim() const;
  String ReplaceAll( const String& key1, const String& key2 ) const;

  static String SPrintf( const char* format, ... );
  static String ReadFileString( const String& filename );
};

#endif //MODYLAS_STRING_H


