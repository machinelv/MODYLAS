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
#include "ModylasTagData.h"
#include "ModylasUtil.h"

#include <iostream>

using namespace std;

// implementation

ModylasTagData::~ModylasTagData()
{
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type == ET_TAG ) {
      delete m_elems[i].tag;
    }
  }
}

vector<String> ModylasTagData::GetStrings() const
{
  vector<String> strs;
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type == ET_STRING ) {
      strs.push_back( m_elems[i].str );
    }
  }
  return strs;
}

String ModylasTagData::GetValue( const String& key ) const
{
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type != ET_KEYVALUE ) continue;
    if( m_elems[i].key  != key         ) continue;
    return m_elems[i].value;
  }
  return String("");
}

void ModylasTagData::SetValue( const String& key, const String& value )
{
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type != ET_KEYVALUE ) continue;
    if( m_elems[i].key  != key         ) continue;
    m_elems[i].value = value;
    return;
  }
  AddKeyValue( key, value );
}

ModylasTagData* ModylasTagData::GetTagData( const String& name ) const
{
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type           != ET_TAG ) continue;
    if( m_elems[i].tag->GetName() != name   ) continue;
    return m_elems[i].tag;
  }
  return NULL;
}

vector<ModylasTagData*> ModylasTagData::GetTagDatas( const String& name ) const
{
  vector<ModylasTagData*> datas;
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type           != ET_TAG ) continue;
    if( m_elems[i].tag->GetName() != name   ) continue;
    datas.push_back( m_elems[i].tag );
  }
  return datas;
}

bool ModylasTagData::AddFromNoTagString( const String& str )
{
  // remove comments and blank line
  String no_tag;
  SIZE_TYPE last_pos = 0;
  while( true ) {
    SIZE_TYPE pos = str.find( "#", last_pos );
    SIZE_TYPE len = ( pos == NPOS ) ? NPOS : pos - last_pos;
    no_tag += str.substr( last_pos, len );
    if( pos == NPOS ) break;
    last_pos = pos + 1;
    pos = str.find( "\n", last_pos );
    if( pos == NPOS ) break;
    last_pos = pos;
  }

  // Strings case 
  if( no_tag.find( "=" ) == NPOS ) {
    // count columns
    vector<String> lines = no_tag.Split( "\n", true );
    for( unsigned int i = 0; i < lines.size(); ++i ) {
      vector<String> tokens = lines[i].Split();
      if( (int)tokens.size() > m_ncolumn ) m_ncolumn = tokens.size();
    }

    vector<String> strs = no_tag.Split();
    for( unsigned int i = 0; i < strs.size(); ++i ) {
      Element elem = { strs[i], "", "", NULL, ET_STRING };
      m_elems.push_back( elem );
    }
    return true;
  }

  // key & value case
  no_tag = no_tag.ReplaceAll( "=" , " = " );
  vector<String> strs = no_tag.Split();
  if( strs.size() % 3 != 0 ) {
    cerr << "ERROR) invalid expression of \"key = value\"." << endl;
    return false;
  }
  int nkeys = strs.size() / 3; 
  for( int i = 0; i < nkeys; ++i ) {
    if( strs[3*i+1] != "=" ) {
      cerr << "ERROR) invalid expression of \"key = value\"." << endl;
      return false;
    }
  }
  for( int i = 0; i < nkeys; ++i ) {
    const String& key   = strs[3*i  ];
    const String& value = strs[3*i+2];
    Element elem = { "", key, value, NULL, ET_KEYVALUE };
    m_elems.push_back( elem );
  }

  return true;
}

bool ModylasTagData::ParseFromString( const String& str )
{
  m_elems.clear();

  SIZE_TYPE last_pos = 0;
  while(true) {
    // search "<"
    SIZE_TYPE pos = str.find("<", last_pos );
    SIZE_TYPE len = ( pos == NPOS ) ? NPOS : pos - last_pos;
    String no_tag = str.substr( last_pos, len );
    bool b = AddFromNoTagString( no_tag );
    if( !b ) {
      return false;
    }
    if( pos == NPOS ) {
      break;
    }
    last_pos = pos + 1;

    // search ">"
    pos = str.find(">", last_pos );
    if( pos == NPOS ) {
      cerr << "ERROR) \">\" not found after \"<\"." << endl;
      return false;
    }

    // "<tag_name>" found, search "</tag_name>"
    String tag_name = str.substr(last_pos, pos - last_pos);
    String open_tag  = String("<" ) + tag_name + ">";
    String close_tag = String("</") + tag_name + ">";

    SIZE_TYPE open_pos = pos + 1;

    pos = str.find(close_tag, open_pos );
    if( pos == NPOS ) {
      cerr << "ERROR) close tag not found for " + open_tag + "." << endl;
      return false;
    }

    // Make String in the tag, and make ModylasTagData recursively
    String str_tag = str.substr( open_pos, pos - open_pos );
    ModylasTagData* tag_data = new ModylasTagData();
    {
      bool b = tag_data->ParseFromString( str_tag );
      if( !b ) {
        delete tag_data;
        return false;
      }
      tag_data->SetName( tag_name );
    }
    Element elem = { "", "", "", tag_data, ET_TAG };
    m_elems.push_back( elem );

    last_pos = pos + close_tag.length();
  }

  return true;
}

String ModylasTagData::WriteDataString( int indent ) const
{
  String ret;
  const String& name = GetName();
  String space;
  for( int i = 0; i < indent; ++i ) {
    space += " ";
  }
  String elemspace = space;

  // write open tag if need
  if( name != "" ) {
    ret += space + "<" + name + ">\n";
    indent += 2;
    elemspace += "  ";
  }

  // write elemets
  unsigned int j = 0;
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    const Element& elem = m_elems[i];
    switch( elem.type ) {
    case ET_STRING:
      if( j % m_ncolumn == 0 ) {
        ret += elemspace;
      }
      else {
        ret += "   ";
      }
      ret += elem.str;
      if( ( j + 1 ) % m_ncolumn == 0 ) {
        ret += "\n";
      }
      ++j;
      break;
    case ET_KEYVALUE:
      ret += elemspace + elem.key + " = " + elem.value + "\n";
      break;
    case ET_TAG:
      ret += elem.tag->WriteDataString( indent );
      break;
    case ET_COMMENT:
      ret += elemspace + "# " + elem.str + "\n";
      break;
    default:
      break;
    }
  }

  // write close tag if need
  if( name != "" ) {
    ret += space + "</" + name + ">\n";
  }
  return ret;
}

ModylasTagData* ModylasTagData::AddTagData( const String& name, int ncolumn )
{
  ModylasTagData* tag = new ModylasTagData();
  tag->SetName( name );
  tag->SetNColumn( ncolumn );

  Element elem = { "", "", "", tag, ET_TAG };
  m_elems.push_back( elem );

  return tag;
}

void ModylasTagData::RemoveTagData( const String& name )
{
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type           != ET_TAG ) continue;
    if( m_elems[i].tag->GetName() != name   ) continue;
    m_elems.erase( m_elems.begin() + i );
    return;
  }
}

void ModylasTagData::AddKeyValue( const String& key, const String& value )
{
  Element elem = { "", key, value, NULL, ET_KEYVALUE };
  m_elems.push_back( elem );
}

void ModylasTagData::AddKeyValue( const String& key, double d )
{
  AddKeyValue( key, String( d ) );
}

void ModylasTagData::AddKeyValue( const String& key, int n )
{
  AddKeyValue( key, String( n ) );
}

void ModylasTagData::RemoveKeyValue( const String& key )
{
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type != ET_KEYVALUE ) continue;
    if( m_elems[i].key  != key         ) continue;
    m_elems.erase( m_elems.begin() + i );
    return;
  }
}

void ModylasTagData::AddString( const String& str )
{
  Element elem = { str, "", "", NULL, ET_STRING };
  m_elems.push_back( elem );
}

void ModylasTagData::AddString( double d )
{
  AddString( String(d) );
}

void ModylasTagData::AddString( int n )
{
  AddString( String(n) );
}

void ModylasTagData::AddComment( const String& str )
{
  Element elem = { str, "", "", NULL, ET_COMMENT };
  m_elems.push_back( elem );
}

void ModylasTagData::RemoveEmptyTag()
{
  for( unsigned int i = 0; i < m_elems.size(); ++i ) {
    if( m_elems[i].type != ET_TAG ) continue;
    if( m_elems[i].tag->Empty() ) {
      delete m_elems[i].tag;
      m_elems.erase( m_elems.begin() + i );
      --i;
    }
    else {
      m_elems[i].tag->RemoveEmptyTag();
    }
  }
}

