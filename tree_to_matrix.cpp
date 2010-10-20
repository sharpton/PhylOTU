#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <string>
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
  // Setup
  if( argc < 10 || argc > 12){
    std::cout << "Wrong number of input arguments (" << argc << "), should have format:\n";
    std::cout << "\ttree_to_matrix <infile> <tmpfile> <prunedfile> <refalignment> <outfile> <starting_row> <ending_row> <format M=matrix E=esprit> <Do_Pruning 0=no 1=yes 2=only prune> [outfile_freq] [maxdistance(E format only)]\n";
  }
  char* infilename     = argv[1];
  char* tempfilename   = argv[2];
  char* prunedfilename = argv[3];
  char* refalignname   = argv[4];
  char* outfilename    = argv[5];
  int startrow    = atoi(argv[6]);
  int endrow      = atoi(argv[7]);
  char format          = argv[8][0];
  int do_pruning  = atoi(argv[9]);
  //    M = matrix format, used by mothur
  //    E = ESPRIT list format
  char* frqfilename;
  float maxdist=0.1;
  if( argc == 12 ){
    frqfilename      = argv[10];
    maxdist     = atof(argv[11]);
    std::cout << frqfilename << " " << maxdist << std::endl;
  } else {
    if( format == 'E' ){
      std::cerr << "maximum distance required for ESPRIT printout; quitting\n";
      return EXIT_FAILURE;
    }
  }
  int srow = startrow;
  char* inname;
  if( do_pruning>0 ){
    // Read in raw file, then prune it
    inname = infilename;
  } else {
    // Read in pruned file directly
    inname = prunedfilename;
  }
  if( format == 'E' ){
    std::cout << "Printing output in ESPRIT list format\n";
  } else if( format == 'M' ){
    std::cout << "Printing output in Mothur matrix format\n";
  } else {
    std::cerr << "Unknown format " << format << ". Quitting\n";
    return EXIT_FAILURE;
  }
  std::list<TreeNode>::iterator startit;
  std::list<TreeNode>::iterator endit;

  /////////////////////////////////////////
  // READ IN TREE FROM FILE
  std::cout << "Reading in " << inname << std::endl;
  PhyloTree<TreeNode>* tr = new PhyloTree<TreeNode>();
  std::ifstream infile;
  infile.open(inname);
  tr->readTree(infile);

  /////////////////////////////////////////
  // Prune tree (if necessary)
  if( do_pruning>0 ){
    std::cout << "Pruning tree\n";
    // Read in reference alignment file and grab reference file names
    std::ifstream reffile;
    reffile.open(refalignname);
    char line[100];
    reffile >> line;
    while( !reffile.eof() ){
      if( line[0] == '>' ){
	// Clean-up the file name
	std::string name(line);
	int slash = (int)name.find("/");
	name = name.substr(1, slash-1);
	// Remove this leaf from the tree
	tr->deleteLeaf(name.c_str());
      }
      reffile >> line;
    }
    reffile.close();

    // Print to tmp file, just in case 
    std::ofstream treeout;
    treeout.open( tempfilename );
    treeout.precision(5);
    treeout.setf(std::ios::fixed,std::ios::floatfield);
    tr->writeTree( treeout );
    treeout.close();
    std::cout << "Printed to file " << tempfilename << std::endl;
    
    // Remove internal nodes that are now leaves
    while( tr->deleteLeaf("") > 0 );

    // Smooth to remove single child nodes
    while( tr->smooth() > 0 );


    // Print pruned file, for use by parallel jobs
    treeout.open( prunedfilename );
    treeout.precision(6);
    treeout.setf(std::ios::fixed,std::ios::floatfield);
    if( !treeout.is_open() ){ std::cout << "Unable to open file " << prunedfilename << std::endl; }
    tr->writeTree( treeout );
    treeout.close();
    std::cout << "Printed to file " << prunedfilename << std::endl;

    // If I only needed to prune then I'm done
    if( do_pruning>1 ){
      std::cout << "Done pruning tips, ready to launch parallel tree_to_matrix jobs\n";
      return EXIT_SUCCESS;
    }
  }

  /////////////////////////////////////////
  // PRINT OUT THE HEADER (only before top row)
  std::ofstream outfile;
  outfile.open(outfilename);
  outfile.precision(5);
  outfile.setf(std::ios::fixed,std::ios::floatfield);

  if( startrow==0 ){
    if( format == 'M' ){
      // phylip format
      outfile << tr->getNleaves() << std::endl;
      //R format (not used)
      //      for( std::list<TreeNode>::iterator it=tr->begin(); it!=tr->end(); it++ ){
      //	if( it->children.size()==0 ){
      //	  outfile << "\"" << it->name.c_str() << "\" ";
      //	}
      //      }
      //outfile << std::endl;
    } else if( format == 'E' ){
      // Print out frequency file, which also provide a map from ID name to index #
      std::ofstream frqfile;
      frqfile.open(frqfilename);
      for( std::list<TreeNode>::iterator it=tr->begin(); it!=tr->end(); it++ ){
	if( it->children.size()==0 ){
	  frqfile << it->name.c_str() << " 1" << std::endl;
	}
      }
      frqfile.close();
    }
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
  if( maxdist < 1 ){
    scale_factor = 1000000;
  }
  size_t eI = 0;
  int nleaf = 0;
  size_t vI = 0;
  // I'm about to use the parent numbers so they better be correct!
  tr->renumber();
  for( std::list<TreeNode>::iterator it=tr->begin(); it!=tr->end(); it++ ){
    name[vI] = it->name.c_str();
    if( it->parents.size() != 0 ){
      edge_array[eI] = Edge( vI, it->parents[0]->number );
      weights[eI] = (int)( it->distance*scale_factor );
      eI++;
    }
    if( it->children.size() == 0 ){
      if( startrow == nleaf ){ startrow = -1*vI;  startit = it;       std::cout << "Effective start row " << -1*startrow << std::endl; }
      if( endrow   == nleaf ){   endrow = -1*vI;  endit=it; endit++;  std::cout << "Effective end   row " << -1*endrow << std::endl; }
      nleaf++;
    }
    vI++;
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
  graph_t h();
  graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
#endif
  std::vector<vertex_descriptor> p(num_vertices(g));
  std::vector<int> d(num_vertices(g));

  
  /////////////////////////////////////////
  // LOOP OVER SELECTED ROWS
  int leafrow=srow;
  int i=startrow;
  for( std::list<TreeNode>::iterator it=startit; it!=endit; it++ ){
    if( it->children.size()==0 ){    // Only care about leaf nodes
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
      if( format == 'M' ){ outfile << "\"" << name[i] << "\" "; }
      int leafcol = 0;
      int j=0;
      for( std::list<TreeNode>::iterator jt=tr->begin(); jt!=tr->end(); jt++ ){
	if( jt->children.size()==0 ){
	  if( format == 'M' ){
	    // Matrix format, print all the leaves in this row
	    outfile << d[j]/scale_factor << " ";
	  } else if( format == 'E' ){
	    // ESPRIT list format, print only index numbers and 
	    if( j>i && (d[j]/scale_factor)<maxdist ){
	      outfile << leafrow << " " << leafcol << " " << d[j]/scale_factor << "\n";
	      }
	  }
	  leafcol++;
	}
	j++;
      }
      if( format == 'M' ){
	outfile << std::endl;
      }
      leafrow++;
    }
    i++;
  }

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  //  std::cout << "Done at " << asctime (timeinfo);
  
  return EXIT_SUCCESS;
}
