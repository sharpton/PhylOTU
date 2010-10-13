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

//typedef unsigned int node_id_t;
typedef size_t node_id_t;
class TreeNode 
{
 public:
  TreeNode() : distance(0) {};
    std::string name;	/**< node name */
    double distance;	/**< distance to parent */
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
  double weight;	/**< Overall tree weight */
  typename std::list<T>::iterator root;	/**< root of the tree */
  std::list< T > nodes;	/**< nodes of the tree */
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
  //	T& operator[]( const unsigned i ){ return nodes[i]; }
  //	const T& operator[]( const unsigned i ) const{ return nodes[i]; }

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
  size_t read_state = 0;	/**< read_state of 0 indicates nothing has been parsed yet */
  size_t section_start = 0;
  std::stack<typename std::list<T>::iterator> node_stack;
  std::stringstream blen_str;
  T new_node;
  new_node.distance = 0;	// default the distance to 0
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
	  read_state = 2;	// pop this node after reading its branch length
	}
	node_stack.top()->name = tree_line.substr( section_start, charI - section_start );
      }
      if( read_state == 2 ){
	node_stack.pop();
      }
      section_start = charI + 1;
      blen_found = false;
      
      // pop off the top of the node stack
      read_state = 2;
      break;
    case ',':
      if( blen_found ){
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
	  read_state = 2;	// pop this node after reading its name
	}
	node_stack.top()->name = tree_line.substr( section_start, charI - section_start );
      }
      if( read_state == 2 )
	node_stack.pop();
      section_start = charI + 1;
      read_state = 1;	// indicates that we'll be creating a new node when we hit :
      blen_found = false;
      break;
    case ':':
      // read off a name, if possible
      if( read_state == 1 ){
	new_node.parents.clear();
	new_node.parents.push_back( node_stack.top() );
	nodes.push_back( new_node );
	new_it = nodes.end(); new_it--;
	node_stack.top()->children.push_back( new_it );
	node_stack.push( new_it );
	read_state = 2;	// pop this node after reading its branch length
      }
      if( tree_line[section_start-1] == ')' ){
	// This is a FastTree score (do nothing for now)
      } else {
	// This is a node name
	node_stack.top()->name = tree_line.substr( section_start, charI - section_start );
      }
      section_start = charI + 1;
      blen_found = true;
      break;
    default:
      break;
    }
  }
  std::cout << "Done reading in tree, WARNING: you may want to run renumber() to number the tree\n";
}

///////////////////////////////////////////////////////////////////////
template< class T >
void PhyloTree<T>::writeAllNodes( std::ostream& os ) const{
  for( typename std::list<T>::const_iterator it=nodes.begin(); it!=nodes.end(); it++ ){
    os << it->number             << "\t" << it->distance << "\t";
    // parents
    os << it->parents.size() << "\t";
    if( it->parents.size() > 0 ){
      os << it->parents[0]->number  << "\t";
    } else {
      os << "\t";
    }
    // children
    os << it->children.size()    << "\t";
    if( it->children.size() > 0 ){
      os << it->children[0]->number << "\t";
    } else {
      os << "\t";
    }
    if( it->children.size() > 1 ){
      os << it->children[1]->number << "\t";
    } else {
      os << "\t";
    }
    os << it->name               << "\t" << std::endl;
  }

}
///////////////////////////////////////////////////////////////////////

template< class T >
void PhyloTree<T>::writeTree( std::ostream& os ) const{
  std::stack<typename std::list<T>::iterator> node_stack;
  std::stack< size_t > child_stack;
  node_stack.push( root );
  child_stack.push( 0 );
  bool write_branch_lengths = false;
  for( typename std::list<T>::const_iterator it=nodes.begin(); it!=nodes.end(); it++ ){
    if( it->distance != 0 ){
      write_branch_lengths = true;
      break;
    }
  }
  if( (*this).weight != 0 ){
    os << "[" << weight << "]";
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
	child_stack.pop();
	continue;
      }

      // try to recurse to its children
      // if the child is a parent as well spit out a parent
      typename std::list<T>::iterator child = node_stack.top()->children[ child_stack.top() ];
      node_stack.push( child );
      child_stack.top()++;
      // print a comma to separate multiple children
      if( child_stack.top() > 1 ){
	os << ",";
      }

      if( child->children.size() > 0 ){
	child_stack.push( 0 );
	os << "(";
      }
      continue;
    }
    
    // this is a leaf node
    os << node_stack.top()->name;
    if( write_branch_lengths ){
      os << ":" << node_stack.top()->distance;
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
double PhyloTree<T>::getHeight( node_id_t nodeI ) const
{
  // Not properly implemented with new version of tree
  if( (*this)[ nodeI ].children.size() == 0 )
    return (*this)[ nodeI ].distance;
  return (*this)[ nodeI ].distance + getHeight( (*this)[ nodeI ].children[ 0 ] );
}
 */
///////////////////////////////////////////////////////////////////////

template< class T >
int PhyloTree<T>::getNleaves()
{
  int nleaves=0;
  for( typename std::list<T>::iterator it=nodes.begin(); it!=nodes.end(); it++ ){
    if( it->children.size() == 0 ){
      nleaves++;
    }
  }
  
  return nleaves;
}

///////////////////////////////////////////////////////////////////////

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

template< class T >
int PhyloTree<T>::deleteLeaf( const char* node_name ){

  int count = 0;
  int named = strlen(node_name);
  if( named ){
    //    std::cout << "\nDeleting node with name " << node_name << std::endl;
  } else {
    //    std::cout << "Deleting all unnamed nodes\n";
  }

  // Loop over leaves until I find this one
  typename std::list<T>::iterator Cit;
  typename std::list<T>::iterator Cit_todelete;
  typename std::vector< typename std::list<T>::iterator >::iterator Pit;
  typename std::vector< typename std::list<T>::iterator >::iterator PCit;
  for (Cit=nodes.begin(); Cit !=nodes.end(); Cit++){
    if( Cit->children.size() == 0 && Cit->name.compare(node_name)==0 ){
      // Remove all references to this node from the parents
      for (  Pit = Cit->parents.begin();  Pit != Cit->parents.end();   Pit++){
	for (PCit = (*Pit)->children.begin(); PCit != (*Pit)->children.end(); PCit++){
	  if( (*PCit) == Cit ){
	    // Delete the reference to this child
	    (*Pit)->children.erase(PCit);
	    break;
	  }
	}
      }
      // Delete the node
      Cit_todelete = Cit;
      Cit--;
      nodes.erase(Cit_todelete);
      count++;
      // If this is a named node we're done
      if( named ){
	return count;
      }
    }
  }
	  
  std::cout << "Removed " << count << " unnamed nodes, WARNING: You may want to renumber() the tree\n";
  return count;
}

///////////////////////////////////////////////////////////////////////

template< class T >
int PhyloTree<T>::smooth(){
  int count = 0;

  typename std::list<T>::iterator it;
  typename std::list<T>::iterator it_todelete;
  typename std::vector< typename std::list<T>::iterator >::iterator Parent;
  typename std::vector< typename std::list<T>::iterator >::iterator Child;
  typename std::vector< typename std::list<T>::iterator >::iterator Pit;
  typename std::vector< typename std::list<T>::iterator >::iterator Cit;
  for (it=nodes.begin(); it !=nodes.end(); it++){
    if( it->children.size() == 1 && it->parents.size() == 1 ){
      // This is a dummy node, change references and remove it
      // std::cout << "Removing node with branch length " << it->distance << std::endl;
      Child  = it->children.begin();
      Parent = it->parents.begin();

      // Set the grandparent's child as it's grandchild
      for (Cit = (*Parent)->children.begin(); Cit != (*Parent)->children.end(); Cit++){
	if( (*Cit) == it ){
	  // This is the reference that needs to be changed
	  *Cit = *Child;
	  break;
	}
      }
      // Set the child's grandparent as it's parent
      for (Pit = (*Child)->parents.begin(); Pit != (*Child)->parents.end(); Pit++){
	if( (*Pit) == it ){
	  // This is the reference that needs to be changed
	  *Pit = *Parent;
	  break;
	}
      }
      // Make the child's distance also into the node we're removing
      (*Child)->distance += it->distance;
      
      // Delete the node
      it_todelete = it;
      it--;
      nodes.erase(it_todelete);
      count++;
    }
  }
	  
  std::cout << "Removed " << count << " 1-child nodes, WARNING: You may want to renumber() the tree\n";
  return count;
}

///////////////////////////////////////////////////////////////////////
/** determine which nodes are descendants of a given node */
template< class TreeType >
void getDescendants( TreeType& alignment_tree, node_id_t node, std::vector< node_id_t >& descendants )
{
  // do a depth first search
  std::stack< node_id_t > node_stack;
  node_stack.push( node );
  descendants.clear();
  while( node_stack.size() > 0 )
    {
      node_id_t cur_node = node_stack.top();
      node_stack.pop();
      if( alignment_tree[cur_node].children.size() > 0 )
	{
	  node_stack.push(alignment_tree[cur_node].children[0]);
	  node_stack.push(alignment_tree[cur_node].children[1]);
	}
      descendants.push_back(cur_node);
    }
}


///////////////////////////////////////////////////////////////////////
/** determine which nodes are leaf nodes below a given node */
template< class TreeType >
void getLeaves( TreeType& tree, node_id_t node, std::vector< node_id_t >& leaves )
{
  // do a depth first search
  std::stack< node_id_t > node_stack;
  node_stack.push( node );
  leaves.clear();
  while( node_stack.size() > 0 )
    {
      node_id_t cur_node = node_stack.top();
      node_stack.pop();
      if( tree[cur_node].children.size() > 0 )
	{
	  node_stack.push(tree[cur_node].children[0]);
	  node_stack.push(tree[cur_node].children[1]);
	}else
	leaves.push_back(cur_node);
    }
}

///////////////////////////////////////////////////////////////////////

namespace std {
  
  template< class T > inline
    void swap( PhyloTree<T>& a, PhyloTree<T>& b )
    {
      a.swap(b);
    }
  
  template<> inline void swap( PhyloTree<TreeNode>& a, PhyloTree<TreeNode>& b){ a.swap(b); }
}

#endif // __PhyloTree_h__
