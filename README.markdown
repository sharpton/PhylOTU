Welcome to PhylOTU
==================

##ANNOUNCEMENTS##

###February 2011###
The PhylOTU manuscript has been published at PLoS Computational Biology (January 2011 edition) and the source code is now officially live. The manuscript describes many aspects of the workflow logic and decisions and demonstrates the method's use, applicability, and validity. I highly recommend reading the manuscript prior to using PhylOTU to ensure its proper implementation. The citation follows:

> Sharpton TJ, Riesenfeld SJ, Kembel SW, Ladau J, O'Dwyer JP, et al. 2011 PhylOTU: A High-Throughput Procedure Quantifies Microbial Community Diversity and Resolves Novel Taxa from Metagenomic Data. PLoS Comput Biol 7(1): e1001061. doi:10.1371/journal.pcbi.1001061

If you used PhylOTU in your study, please cite the above manuscript. If you have any questions that the manuscript and this documentation fail to address, please contact the author (Thomas Sharpton - thomas.sharpton@gladstone.ucsf.edu). Additional help may be found via the comments embedded in the source code, principally otu_handler.pl and OTU.pm.

TJS

###August 2010###
This software is still in beta and is subject to change. More specific documentation is being worked up. 
In the meantime, please see the source code for more information (otu_handler.pl and OTU.pm) and 
contact the author (Thomas Sharpton - thomas.sharpton@gladstone.ucsf.edu) with specific inquiries.

Keep an eye out for the PhylOTU manuscript, which is currently in review.

TJS

##INTRODUCTION##

PhylOTU is a software workflow that identifies OTUs directly from metagenomic reads. Because of the 
particular methodology employed in PhylOTU, reads that share no sequence overlap (i.e., cannot be 
aligned to one another) can still be clustered into the same OTU. In addition, reads are not recruited into
top BLAST hits, but are effectively clustered de novo.

Briefly, PhylOTU uses a CMmodel of the 16S locus (via INFERNAL) to align metagenomic reads. The CMmodel
is contrusted using a reference alignment of high-quality, full-length 16S sequences. A phylogenetic tree
is constructed from the multiple sequence alignment of reads and reference sequences via FastTree. The
phylogenetic distance spanning every pair of reads is then calculated from the resulting tree and used 
to create a distance matrix describing pairwise read relationships. This distance matrix is then fed into 
MOTHUR, which implements average-neighbor clustering to identify OTUs. The final output (a list of OTU clusters)
is a file in the matrix subdirectory of the metagenomic sample set with the following typed name:

> <path_to_database>/samples/<SAMPLE_NAME>/matrix/<SAMPLE_NAME>.an.list

##REQUIREMENTS##

There are lots of dependencies. Please read carefully to ensure your system is properly configured. I have included links to the software in question. Please contact me if you find these are dead. Don't forget your local open source developer for contributing to this infrastructure.

PhylOTU is primarily written in Perl 5 and R. Run time dependencies include the following Perl packages:

-  [Bioperl-live libraries](https://github.com/bioperl/bioperl-live)
-  [IPC::System::Simple](http://search.cpan.org/~pjf/IPC-System-Simple-1.21/lib/IPC/System/Simple.pm)
-  [Bio::Phylo::IO](http://search.cpan.org/~rvosa/Bio-Phylo-0.34/lib/Bio/Phylo/IO.pm) 
-  [File::Basename](http://search.cpan.org/~rjbs/perl-5.12.3/lib/File/Basename.pm)
-  [File::Path](http://search.cpan.org/~dland/File-Path-2.08/Path.pm)
-  [File::Copy](http://search.cpan.org/~rjbs/perl-5.12.3/lib/File/Copy.pm)

In addition, you must have the following R library installed:

-  [APE](http://cran.r-project.org/web/packages/ape/index.html) Note: easiest install is via the R command install.packages

As well as the following C++ library:

-  [Boost](http://www.boost.org/)

PhylOTU stitches together various software packages written by other authors. While the software can easily be modified to accomodate various software suites, it is currently designed to implement the following tools, which you will also need to install on your system. The use of methods alternative to those listed here should be coupled with an independent validation test as described in the PhylOTU manuscript

-  [INFERNAL](http://infernal.janelia.org/)
-  [FastTree](http://www.microbesonline.org/fasttree/)
-  [BLAST legacy version](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/)
-  [MOTHUR](http://www.mothur.org/wiki/Main_Page)

Please see the authors' websites for instructions regarding the installation of any of the above software.

##INSTALLATION AND SETUP##
###A. C++ compilation###

PhylOTU uses two c++ programs to efficiently process large data sets (align2profile_qc_Col.cpp and tree_to_matrix.cpp). Compile these programs by navigating to the directory that contains the source code and running the following commands at the command line:

> g++ align2profile_qc_Col.cpp -o align2profile_qc_Col

> g++ -I /usr/include/boost/ tree_to_matrix.cpp -o tree_to_matrix

###B. Reference Data and CMmodel building###

PhylOTU leverages reference 16S sequence data (high-quality, full-length) to construct an evolutionary model via the cmbuild function in INFERNAL. 

Because of the particular licensing agreement, we cannot currently distribute the reference data we use in our study with this software (though we are working with the publishers of the reference data to change this). We recommend contacting the Ribosomal Database Project to obtain the reference library described in the RDP v 10 manuscript:

> Cole et. al. The Ribosomal Database Project: improved alignments and new tools for rRNA analysis. Nucleic Acids Res. 2009 Jan;37(Database issue):D141-5. Epub 2008 Nov 12.

You can alternatively create your own 16S reference alignment (that includes secondary structural information), though we would recommend independengly evaluating the model you build from this alignment. 

The reference alignment is used to build a CMmodel using INFERNAL's cmbuild program. This model enables sensitive alignment of fragmented reads that are phylogenetically diverse from the full-length reference sequences within the context of the reference multiple sequence alignment. The specific cmbuild run-time parameters depend on various characteristics of the reference alignment. If you use the RDP reference alignment mentioned above, the following command will build a CMmodel (the same we use in the manuscript): 

> cmbuild --rf --ere 1.4 <name of model> <reference alignment file>

For the above command to work, your reference alignment file must be in stockholm format. Alternative alignments may require different settings, but will generally look like the above command. PhylOTU will build the model for you (using the above parameters as defaults, see section C below for more), but it only needs to occur once. Since this is a time intensive process, ensure that you not telling PhylOTU to unnecessairly rebuild a model. When you use the *-first* parameter (see section C, above), PhylOTU will automatically try to build the model for you, so odds are good you'll never need to worry about this. If you do this by hand, any models you construct should be placed in 

> <database>/reference/profiles/ 

and the corresponding reference alignment should be placed in 

> <database>/reference/aligns
 
###C. Build the flat file database###

PhylOTU organizes the sample, reference, and workflow output data via the use of a flat file database. While you may specify the root location of the database, the subdirectory structure is controlled by PhylOTU. Once you have designated a location for your database, you need to initialize it. This involves creating the organizational hierarchy, building the STAP blast databases, and training the 16S CMmodel from the reference alignment (see section B, above). Some of these steps take a few moments (dependent on reference sequence database size), but only need to be conducted once. With the infrastructure established, PhylOTU can run relatively quickly on very large metagenomic datasets. 

Prior to initialization, you'll need to set a few variables: the database path, the PhylOTU code path, and the path to your reference alignment (optional). To set the database path, open otu_handler.pl and point the variable $masterdir to the location you want to install your database. Alternatively, run PhylOTU at the command line using the *-db <path_to_database>* parameter. You will also need to hardcode the location of the PhylOTU source code on your computer in otu_hander.pl. Specifically, point the variable $scripts_path to your code. Alternatively, use the *-sd <path_to_code>* parameter at run time. It is recommended that you hardcode these location within otu_handler.pl since this is unlikely to change once installed. To point to the reference alignment, use the *-ra <path_to_reference_alignment>* parameter at run-time.

To initialize the database and have PhylOTU build your CMmodel from a reference alignment of your choosing, run the following command at the command line in the directory where you have the PhylOTU code installed (this builds the bacteria model):
 
> perl otu_handler.pl -first -ra <path_to_bacteria_reference_alignment> -bac

to build the archaeal model, use the following command:

> perl otu_handler.pl -first -ra <path_to_archaea_reference_alignment> -arc

If you didn't hardcode the database and code paths in otu_handler.pl, you'd amend the above statements with the appropriate command line options. For the bacteria example:

> perl otu_handler.pl -first -ra <path_to_reference_alignment> -db <path_to_database> -sd <path_to_code>

This command will initialize the aforementioned database and reference files, including the CMmodel. If you want to build multiple CMmodels (e.g., one for bacteria and one for archaea), simply run the first of the two above commands multiple times, pointing *-ra* to a different reference alignment each time. If you want to initialize the database and then build a CMmodel by hand (see section B, above), run the following command in the same directory:

> perl otu_handler.pl -first

Note that you will be unable to use PhylOTU until you build a cmmodel.

Generally speaking, you are now ready to process a metagenomic library with PhylOTU.

##IMPLEMENTATION##

PhylOTU is run directly at the command line and is implemented via a single command:

> perl otu_handler.pl -bac -arc -i <metagenomic sequence file> 

That said, there are many settings that need to be controlled for successful implementation of PhylOTU. These settings are controlled through the handler script otu_handler.pl under the section USER RUNTIME OPTIONS. Please see the source code of otu_hander.pl for detailed information about the various settings and their use. The output of PhylOTU is a database of files in addition to
run-time logs. You may elect to pipe these logs to a file to reduce screen clutter with the follwing command:

> perl otu_handler.pl -i <metagenomic sequence file> > standard.out 2> error.log

Please send all bug reports and inquiries to the author (Thomas Sharpton - thomas.sharpton@gladstone.ucsf.edu).
