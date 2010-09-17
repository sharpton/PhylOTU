#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include "PhyloTree.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

  using namespace boost;

int main(int argc, char *argv[])
{

  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  //  std::cout << "Start " << asctime (timeinfo);

  /////////////////////////////////////////
  // READ IN TREE FROM FILE
  char* infilename  = argv[1];
  char* outfilename = argv[2];
  int startrow = atoi(argv[3]);
  int endrow   = atoi(argv[4]);
  std::cout << "Reading in " << infilename << std::endl;
  PhyloTree<TreeNode>* tr = new PhyloTree<TreeNode>();
  std::ifstream infile;
  infile.open(infilename);
  tr->readTree(infile);
  std::ofstream outfile;
  outfile.open(outfilename);
  outfile.precision(5);

  /////////////////////////////////////////
  // PRINT OUT THE HEADER (only before top row)
  if( startrow==0 ){
    for( int j=0; j < tr->size(); j++ ){
      if( (*tr)[j].children.size()==0 ){
	outfile << "\"" << (*tr)[j].name.c_str() << "\" ";
      }
    }
    outfile << std::endl;
  }

  /////////////////////////////////////////
  // DEFINE SOME STUFF
  typedef adjacency_list < listS, vecS, undirectedS,
    no_property, property < edge_weight_t, int > > graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
  typedef std::pair<int, int> Edge;


  /////////////////////////////////////////
  // FILL IN THE EDGE ARRAY
  const int num_nodes = tr->size();
  Edge edge_array[num_nodes-1];   // every node has a parent, except the root
  int weights[num_nodes-1];
  const char* name[num_nodes];
  std::cout << "Startrow is " << startrow << std::endl; 
  if( endrow >= tr->getNleaves() || endrow==0 ){ 
    endrow = tr->getNleaves()-1; 
    std::cout << "Endrow set to " << endrow << std::endl; 
  } else {
    std::cout << "Endrow is " << endrow << std::endl; 
  }
  double scale_factor = 1;
  if( (*tr)[1].distance < 1 ){
    scale_factor = 1000000;
  }
  size_t eI = 0;
  int nleaf = 0;
  for( size_t vI = 0; vI < num_nodes; ++vI ){
    name[vI] = (*tr)[vI].name.c_str();
    if( (*tr)[vI].parents.size() != 0 ){
      edge_array[eI] = Edge( vI, (*tr)[vI].parents[0] );
      weights[eI] = (int)( (*tr)[vI].distance*scale_factor );
      eI++;
    }
    if( (*tr)[vI].children.size() == 0 ){
      if( startrow == nleaf ){ startrow = -1*vI;  std::cout << "Effective start row " << -1*startrow << std::endl; }
      if( endrow   == nleaf ){   endrow = -1*vI;  std::cout << "Effective end   row " << -1*endrow << std::endl; }
      nleaf++;
    }
  }
  startrow *= -1;
  endrow   *= -1;


  /////////////////////////////////////////
  // MAKE A GRAPH
  int num_arcs = sizeof(edge_array) / sizeof(Edge);
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
  graph_t g(num_nodes);
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
  for (std::size_t j = 0; j < num_arcs; ++j) {
    edge_descriptor e; bool inserted;
    boost::tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
    weightmap[e] = weights[j];
  }
#else
  graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
#endif
  std::vector<vertex_descriptor> p(num_vertices(g));
  std::vector<int> d(num_vertices(g));

  
  /////////////////////////////////////////
  // LOOP OVER SELECTED ROWS
  for( int i=startrow; i<=endrow; i++ ){
    if( (*tr)[i].children.size()==0 ){

      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      //      std::cout << "Row " << i << " start at " << asctime (timeinfo);
      /////////////////////////////////////////
      // RUN DIJKSTRA
      vertex_descriptor s = vertex(i, g);
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
      // VC++ has trouble with the named parameters mechanism
      property_map<graph_t, vertex_index_t>::type indexmap = get(vertex_index, g);
      dijkstra_shortest_paths(g, s, &p[0], &d[0], weightmap, indexmap, 
			      std::less<int>(), closed_plus<int>(), 
			      (std::numeric_limits<int>::max)(), 0,
			      default_dijkstra_visitor());
#else
      dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));
#endif

      /////////////////////////////////////////
      // PRINT OUT THE ROW
      outfile << "\"" << name[i] << "\" ";
      for( int j=0; j<num_vertices(g); j++ ){
	if( (*tr)[j].children.size()==0 ){
	  outfile << d[j]/scale_factor << " ";
	}
      }
      outfile << std::endl;
    }
  }

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  //  std::cout << "Done at " << asctime (timeinfo);
  
  return EXIT_SUCCESS;
}
