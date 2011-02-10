#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

/*align2profile_qc_Col.cpp - alignment column cleanup
Copyright (C) 2011  Thomas J. Sharpton 
author contact: thomas.sharpton@gladstone.ucsf.edu

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

int width  = 8000;    // max number of columns in alignment
int MAXNUM = width*8; // max number of characters in %coverage line
int main(int argc, char *argv[]){

  char* infaname   = argv[1];
  char* ingapname  = argv[2];
  char* outfaname  = argv[3];
  char* outgapname = argv[4];
  int min     = atoi(argv[5]);

  cout << "Reading  input alignment " << infaname <<  "\tand  input gap file " << ingapname << endl;
  cout << "Writing output alignment " << outfaname << "\tand output gap file " << outgapname << endl;

  int printedchars=0;
  int badind=0;
  int badchar=0;
  int i=0;
  int num;
  int MAXCHAR=60; // Number of characters per line in the fasta file
  int printed=0;
  float cov=0;
  char line[MAXNUM];
  char achar;
  int bad[width];
  ifstream infile;
  ofstream outfile;

  ////////////////////////////////////////////////////////////////////////////
  // Read in and print out gap file and collect bad clounms along the way
  infile.open(ingapname);
  outfile.open(outgapname);
  if( !infile.is_open() ){  cout << "Failed to open input gap file, quitting\n";  return 1; }
  if( !outfile.is_open() ){ cout << "Failed to open output gap file, quitting\n"; return 1; }
  ////////////////////////////////////////////////////
  // Keep the firstline as is, save as the denominator
  infile.getline(line, MAXNUM);
  outfile << line << endl;
  /////////////////////////////////////////////////////////////////
  // The second line is the coverage, save the list of bad indicies
  i=0;
  badind=0;
  printed=0;
  //  infile >> line;
  infile.getline(line, MAXNUM);
  istringstream in( line );
  while( in >> num ){
    if( num >= min ){
      //This is a sufficiently filled column, preserve it
      outfile << num << " ";
      printed++;
    } else {
      //Mark this column to be removed
      bad[badind] = i;
      badind++;
    }
    i++;
  }
  outfile << endl;
  bad[badind] = i;
  cout << "read in gap file with " << i << " columns, " << badind << " bad columns and printed " << printed << " cloumns\n";
  //////////////////////////////////////////////////
  // The third line is the % coverage, print selected
  i=0;
  badind = 0;
  printed=0;
  badchar = bad[badind];
  outfile.precision(5);
  outfile.setf(ios::fixed);
  infile >> cov;
  while( !infile.eof() ){
    if( i < badchar ){
      outfile << cov << " ";
      printed++;
    } else {
      badind++;
      badchar = bad[badind];
    }
    infile >> cov;
    i++;
  }
  cout << "read in gap file with " << i << " columns, " << badind << " bad columns and printed " << printed << " cloumns\n";
  ////////////
  // Finish up
  infile.close();
  outfile.close();


  ////////////////////////////////////////////////////////////////////////////
  // Now read though the alignment and print only good columns
  i = -1;
  infile.open(infaname);
  outfile.open(outfaname);
  if( !infile.is_open() ){  cout << "Failed to open input alignment file, quitting\n";  return 1; }
  if( !outfile.is_open() ){ cout << "Failed to open output alignment file, quitting\n"; return 1; }
  //  infile.getline(line, MAXNUM);
  infile >> line;
  while( !infile.eof() ){
    if( line[0] == '>' ){
      if( i != -1 ){
	// Output a newline if this isn't the first sequence
	outfile << endl;
      }
      //Print the new sequence identifier
      outfile << line << endl;
      i = 0;
      printedchars = 0;
      printed = 0;
      badind = 0;
      badchar = bad[badind];
    } else {
      // This is a part of the current sequence
      istringstream in( line );
      while( in >> achar ){
	if( i < badchar ){
	  // Good char, print it
	  outfile << achar;
	  printed++;
	  printedchars++;
	  if(printedchars == MAXCHAR) { 
	    // Done printing this line, start a new one
	    outfile << endl;
	    printedchars = 0;
	  }
	} else {
	  // This is the badchar, skip it
	  badind++;
	  badchar = bad[badind];
	}
	i++;
      }
    }
    infile >> line;
  }
  infile.close();
  outfile.close();
  
}
