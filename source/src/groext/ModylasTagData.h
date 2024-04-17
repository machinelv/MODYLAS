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
  @brief  Utilities for dealing with MODYLAS files
* @author  Tetsuro Nagai
*/
#ifndef MODYLAS_TAG_DATA_H
#define MODYLAS_TAG_DATA_H

#include "String.h"
#include <vector>

using namespace std;

// Modylas Data Class
class ModylasTagData
{  
  typedef enum {
    ET_STRING,
    ET_KEYVALUE,
    ET_TAG,
    ET_COMMENT,
  } ElementType;
    
  typedef struct {
    String str;
    String key;
    String value;
    ModylasTagData* tag;
    ElementType type;
  } Element;

public:
  ModylasTagData() : m_ncolumn(1) {};
  ~ModylasTagData();
    
public:
  const String& GetName() const { return m_name; }
  void SetName( const String& name ) { m_name = name; }
  int GetNColumn() const { return m_ncolumn; }
  void SetNColumn( int ncolumn ) { m_ncolumn = ncolumn; }
  vector<String> GetStrings() const ;
  String GetValue( const String& key ) const ;
  void SetValue( const String& key, const String& value );
  ModylasTagData* GetTagData( const String& name ) const ;
  vector<ModylasTagData*> GetTagDatas( const String& name ) const ;
  bool ParseFromString( const String& str );
  String WriteDataString( int indent = 0 ) const ;
  ModylasTagData* AddTagData( const String& name, int ncolumn = 1 );
  void RemoveTagData( const String& name );
  void AddKeyValue( const String& key, const String& value );
  void AddKeyValue( const String& key, double d );
  void AddKeyValue( const String& key, int n );
  void RemoveKeyValue( const String& key );
  void AddString( const String& str );
  void AddString( double d );
  void AddString( int n );
  void AddComment( const String& str );
  void Clear() { m_elems.clear(); }
  bool Empty() { return m_elems.empty(); }
  void RemoveEmptyTag();

private:
  bool AddFromNoTagString( const String& str );

private:
  String m_name;
  int m_ncolumn;
  vector<Element> m_elems;
};


#endif //MODYLAS_TAG_DATA_H


