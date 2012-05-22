#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __PhyloTree_h__
#define __PhyloTree_h__

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stack>

/*PhyloTree.h - Tree parsing and navigation routines
Copyright (c) 2005-2010 Aaron Darling 
This file based on PhyloTree.h, from the libMems library
See http://sourceforge.net/projects/mauve for more information
about libMems

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see LICENSE.txt).  If not, see 
<http://www.gnu.org/licenses/>.
*/

				    //typedef unsigned int node_id_t;
typedef size_t node_id_t;
class TreeNode 
{
 public:
 TreeNode() : distance(0) {};
  std::string name;/**< node name */
  double distance;/**< distance to parent */
  std::vector< std::list<TreeNode>::iterator > parents;   /**< if parents.size() == 0 this is a root node */
  std::vector< std::list<TreeNode>::iterator > children;  /**< if children.size() == 0 this is a leaf node */
  int number;
};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

template< class T >
class PhyloTree 
{
 public:
  PhyloTree();
  PhyloTree( const PhyloTree<T>& pt );
  PhyloTree<T>& operator=( const PhyloTree<T>& pt );
  double weight;/**< Overall tree weight */
  typename std::list<T>::iterator root;/**< root of the tree */
  std::list< T > nodes;/**< nodes of the tree */
  void clear();
  /**
   * Reads a tree in Newick format.  WARNING:  only reads rooted trees correctly
   */
  void readTree( std::istream& tree_file );
  /**
   * Writes a tree in Newick format
   */
  void writeTree( std::ostream& os ) const;
  /**
   * Determines the height of the tree along the path from the root to the left-most leaf node
   */
  //  double getHeight() const;
  /**
   * Determines the height of the tree along the path from nodeI to its left-most descendant leaf node
   */
  //  double getHeight( node_id_t nodeI ) const;
  /**
   * Counts the number of nodes that are leaves
   */
  int getNleaves();
  /**
   * Count the number of leaves
   */
  int deleteLeaf( const char*  node_name );
  /**
   * Delete the (first) leaf with the name node_name; 
   * if node_name="", delete all leaves with no name 
   * (internal nodes that have turned into leaves)
   */
  int smooth();
  /**
   * Removed nodes with only one child and add distance to the child
   */
  int check_root();
  /**
   * Check that the root has two children, if not, change the root
   */
  //T& operator[]( const unsigned i ){ return nodes[i]; }
  //const T& operator[]( const unsigned i ) const{ return nodes[i]; }

  void renumber();
  void writeAllNodes( std::ostream& os ) const;
  
  size_t size() const{ return nodes.size(); }
  typename std::list<T>::iterator begin(){ return nodes.begin(); }
  typename std::list<T>::iterator end(){   return nodes.end(); }
  void push_back( T& t ){ nodes.push_back(t); }
  T& back() { return nodes.back(); }
  const T& back() const{ return nodes.back(); }
  void resize( const unsigned s ){ nodes.resize(s); }
  
  void swap( PhyloTree<T>& other ){
    std::swap( weight, other.weight );
    std::swap( root, other.root );
    nodes.swap( other.nodes );
  }

 protected:

};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

template< class T >
PhyloTree<T>::PhyloTree()
{
  weight = 0;
  root = nodes.end();
}

///////////////////////////////////////////////////////////////////////

template< class T >
PhyloTree<T>::PhyloTree( const PhyloTree<T>& pt ) :
nodes( pt.nodes ),
  weight( pt.weight ),
  root( pt.root )
{}

///////////////////////////////////////////////////////////////////////

template< class T >
PhyloTree<T>& PhyloTree<T>::operator=( const PhyloTree<T>& pt )
{
  nodes = pt.nodes;
  weight = pt.weight;
  root = pt.root;
  return *this;
}

///////////////////////////////////////////////////////////////////////

template< class T >
void PhyloTree<T>::clear()
{
  nodes.clear();
  weight = 0;
  root = nodes.end();
}

///////////////////////////////////////////////////////////////////////
/**
 *  readTree version 2.0: read in a phylogenetic tree in the Newick file format.
 *
 */
template< class T >
void PhyloTree<T>::readTree( std::istream& tree_file )
{
  std::string line;
  clear();
  if( !std::getline( tree_file, line ) )
    return;
  // look for either a ; or a matched number of parenthesis, if
  // not found then read another line
  while(true){
    int paren_count = 0;
    for( size_t charI = 0; charI < line.size(); charI++ )
      {
	if( line[charI] == '(' )
	  paren_count++;
	if( line[charI] == ')' )
	  paren_count--;
      }
    if( paren_count == 0 )
      break;
    if( paren_count != 0 ){
      std::string another_line;
      if( !std::getline( tree_file, another_line ) )
	return;
      line += another_line;
    }
  }
  
  std::stringstream line_str( line );
  
  // look for a weight
  std::string::size_type open_bracket_pos = line.find( "[" );
  std::string::size_type bracket_pos = line.find( "]" );
  if( open_bracket_pos != std::string::npos && bracket_pos != std::string::npos && 
      open_bracket_pos < bracket_pos && bracket_pos < line.find( "(" ) ){
    // read in a weight
    getline( line_str, line, '[' );
    getline( line_str, line, ']' );
    std::stringstream weight_str( line );
    weight_str >> weight;
  }
  
  // ready to begin parsing the tree data.
  std::list<TreeNode>::iterator new_it;
  std::string tree_line;
  std::getline( line_str, tree_line, ';' );
  size_t read_state = 0;/**< read_state of 0 indicates nothing has been parsed yet */
  size_t section_start = 0;
  std::stack<typename std::list<T>::iterator> node_stack;
  std::stringstream blen_str;
  T new_node;
  new_node.distance = 0;// default the distance to 0
  bool already_read_name = false;
  bool blen_found = false;
  for( size_t charI = 0; charI < tree_line.size(); charI++ ){
    switch( tree_line[ charI ] ){
      // if this is an open parens then simply create a new
      // parent node and push it on the parent stack
    case '(':
      if( node_stack.size() > 0 ){
	new_node.parents.clear();
	new_node.parents.push_back( node_stack.top() );
      }
      nodes.push_back( new_node );
      new_it = nodes.end(); new_it--;
      if( node_stack.size() > 0 ){
	node_stack.top()->children.push_back( new_it );
      } else {
	// This is the first node, set it as root
	std::cout << "Setting root\n";
	root = new_it;
      }
      node_stack.push( new_it );
      read_state = 1;
      section_start = charI + 1;
      break;
    case ')':
      if( blen_found )
	{
	  // read off a branch length
	  blen_str.clear();
	  blen_str.str( tree_line.substr( section_start, charI - section_start ) );
	  blen_str >> node_stack.top()->distance;
	}else{
	// read off a name, if possible
	if( read_state == 1 ){
	  new_node.parents.clear();
	  new_node.parents.push_back( node_stack.top() );
	  nodes.push_back( new_node );
	  new_it = nodes.end(); new_it--;
	  node_stack.top()->children.push_back( new_it );
	  node_stack.push( new_it );
	  read_state = 2;// pop this node afstart, charI - section_start ) );
	  blen_str >> node_stack.top()->distance;
	}else{
	  // read off a name, if possible
	  if( read_state == 1 ){
	    new_node.parents.clear();
	    new_node.parents.push_back( node_stack.top() );
	    nodes.push_back( new_node );
	    new_it = nodes.end(); new_it--;
	    node_stack.top()->children.push_ba;
	    read_state = 2;// pop this node after reading its branch length
	  }
	  if( tree_line[section_start-1] == ')' ){
	    // This is a FastTree score (do nothing for now)
	  } else {
	    // This is a node name
	    node_stack.top()->name = tree_line.su( it->parents.size() > 0 ){
	      os << it->parents[0]->number  << "\t";
	    } else {
	      os << "\t";
	    }
	    // children
	    os << it->children.size()    << "\t";
	    if( it->children.size() > 0 ){
	      os << it->children[0]->number << "\t";
	    } else {
	      os <<< "]";
	    }
	    os << "(";

	    while( node_stack.size() > 0 ) {
	      if(  node_stack.top()->children.size() != 0 ){
		// this is a parent node
		// if we have scanned all its children then pop it
		if( child_stack.top() == node_stack.top()->children.size() ){
		  os << ")";
		  if( node_stack.size() > 1 && write_branch_lengths ){
		    os << ":" << node_stack.top()->distance;
		  }
		  node_stack.pop();
		  child_< node_stack.top()->distance;
		}
		// pop the child
		node_stack.pop();
	      }
	      os << ";" << std::endl;
	    }

	    ///////////////////////////////////////////////////////////////////////
	    /*
template< class T >
double PhyloTree<T>::getHeight() const
{
  return getHeight( root );
}
	    */
	       ///////////////////////////////////////////////////////////////////////
	       /*
template< class T >
double Phylo/////////////////////

template< class T >
void PhyloTree<T>::renumber(){
  int i=0;
  for( typename std::list<T>::iterator it=nodes.begin(); it!=nodes.end(); it++ ){
    it->number = i;
    i++;
  }
  std::cout << "Tree sucessfully renumbered\n";
}

///////////////////////////////////////////////////////////////////////
/// WARNING!!! Deleting nodes made convert an unrooted tree into a rooted tree ///
template< c+){
    if( Cit->children.size() == 0 && Cit->name.compare(node_name)==0 ){
      // Remove all references to this node from the parents
      for (  Pit = Cit->parents.begin();  Pit != Cit->parents.end();   Pit++){
      for (PCit = (*Pit)->children.begin(); PCit != (*Pit)->children.end(); PCit++){
        if( (*PCit) == Cit ){
	    // Delete the reference to this child
	        (*Pit)->children.erase(PCit);
		    brt<T>::iterator >::iterator Parent;
  typename std::vector< typename std::list<T>::iterator >::iterator Child;
  typename std::vector< typename std::list<T>::iterator >::iterator Pit;
  typename std::vector< typename std::list<T>::iterator >::iterator Cit;
  for (it=nodes.begin(); it !=nodes.end(); it++){
    if( it->children.size() == 1 && it->parents.size() == 1 ){
      // This is a dummy node, change references and remove it
      // std::cout << "Removing node with branch length " << it->distance << std::endl;
      Child  = it->children.begin();
      Parent = it->parents.begin();

      /e reference that needs to be changed
      *Pit = *Parent;
        break;
	}
      }
      // Make the child's distance also into the node we're removing
      (*Child)->distance += it->distance;
      
      // Delete the node
      it_todele
  }
  typename std::list<T>::iterator newroot;
  // Set the new root to be the child of the old root
  newroot = *(root->children.begin());
  // The new root should have no parents
  newroot->parents.clear();
  // Erase the old root, reset the root variable
  nodes.erase(root);
  root = newroot;
  std::cout << "Reset the root\n";

  return 0;
}

///////////////////////////////////////////////////////////////////////
/** determine which nodesen.size() > 0 )
{
  node_stack.push(alignment_tree[cur_node].children[0]);
    node_stack.push(alignment_tree[cur_node].children[1]);
    }
      descendants.push_back(cur_node);
    }
}


///////////////////////////////////////////////////////////////////////
/** determine which nodes are leaf nodes below a given node */
///////////////////////////////////////

	    namespace std {
  
  template< class T > inline
    void swap( PhyloTree<T>& a, PhyloTree<T>& b )
    {
      a.swap(b);
    }
  
  template<> inline void swap( PhyloTree<TreeNode>& a, PhyloTree<TreeNode>& b){ a.swap(b); }
	    }

#endif // __PhyloTree_h__
