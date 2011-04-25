// Written by Avinash using Andre's Library Feb 11 2011
// sample usage - ./ancestral_parsimony -t input.txt -o output.txt
// This is a program to find all possible rooted trees.
extern const char *builddate;
#include "common.h"
#include "argument.h"
#include "tree.h"
#include "tree_IO.h"
#include "tree_traversal.h"
#include "tree_subtree_info.h"
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/foreach.hpp>
#include <math.h>
using namespace std;

// -----------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------

int main(int ac, char* av[]) 
{
  {
    std::ostringstream os; os << "command:";
    for (int i = 0; i < ac; i++) os << ' ' << av[i];
    MSG(os.str());
  }
  std::string tree_filename;
  std::string seqtree_filename;
  std::string tab_filename;
  std::string output_filename;
  {
    Argument a; a.add(ac, av);
    // help
    if (a.existArg2("-h","--help"))  
      {
	MSG("options:");
	MSG("  -t [ --tree ] filename       tree (file in plain NEWICK format)");
	MSG("  -ct [ --chartree ] filename  characters as taxa labels in tree form");
	MSG("  -cl [ --chartab ] filename   characters in table form");
	MSG("  -o [ --output ] filename     output to file instead of screen");
	MSG("  -v [ --version ]             date of compilation");
	MSG("  -h [ --help ]                produce help message");
	MSG("");
	MSG("example:");
	MSG("  " << av[0] << " -i trees.newick -c 'acient_fox'");
	exit(0);
      }
    // version
    if (a.existArg2("-v","--version")) 
      {
	MSG("date of compilation: " << builddate);
	exit(0);
      }
    // input tree
    if (a.getArgVal2("-t", "--tree", tree_filename)) MSG("tree file: " << tree_filename);
    // characters as taxa labels in tree form
    if (a.existArgVal2("-ct", "--chartree", seqtree_filename)) MSG("character tree file: " << seqtree_filename);
    // characters in table form
    if (a.existArgVal2("-cl", "--chartab", tab_filename)) MSG("character table file: " << tab_filename);
    // output file
    if (a.existArgVal2("-o", "--output", output_filename)) MSG("output file: " << output_filename);
    // unknown arguments?
    a.unusedArgsError();
    //cout<<"1\n";
  }

  // -----------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------
  //MSG("2");//avi_test
  // read trees
  std::ifstream ifs,ifs2;
  const std::string filename1 = tree_filename;
  ifs.open (filename1.c_str());
  std::istream &is = filename1.empty() ? std::cin : ifs;
  ifs2.open(filename1.c_str());
  std::istream &is2 = filename1.empty() ? std::cin : ifs2;
  const std::string &filename2 = output_filename;
  std::ofstream ofs(filename2.c_str());
  std::ostringstream os,os2;	
  //int k=0;
  //string tree;
  //getline(ifs,tree); 
  while(ifs.good())
  {
    aw::Tree t,t2;
    aw::idx2name t_name,t_name2;
    aw::idx2weight_double t_weight,t_weight2; 
    //cout<<"\nk = "<<++k;		     
    if (!aw::stream2tree(is, t, t_name, t_weight))     break;			      
    if (!aw::stream2tree(is2, t2, t_name2, t_weight2)) break; 
    //getline(ifs,tree); 
    		       
    aw::tree2newick(os, t, t_name, t_weight);
    os << std::endl;  
    
			      

    //int tree_count = 0;
    TREE_POSTORDER2(k2,t2) //if(t_temp.is_leaf(k.idx)) 
      {
				 
	//cout<<"\nNode Count"<<++tree_count<<endl;
	cout<<"\nRoot is"<<t2.root<<endl;      	      
	unsigned int root_node = t.root;
			      	      
	if(k2.idx != t.root)
	  {
	    //cout<<"\nidx"<<k2.idx;
	    unsigned int l[3],k[3],i = 0,j =0;
			      	        
	    BOOST_FOREACH(j, t.adjacent(k2.idx)) 
	      {
		//cout<<"1";
		l[i] = j;
		i++;
	      }
	    i=0;
	    if((l[0]!=0 && l[1]!=0) || (l[0]!=0 && t.is_leaf(k2.idx)))
	      {
		BOOST_FOREACH(j, t.adjacent(t.root)) 
		  {
		    //cout<<"\nROOT ADJACENT"<<j;
		    t.remove_edge(j,t.root);
		    k[i] = j;
		    i++;
		  }
		
		t.add_edge(k[0],k[1]);      	        
		//cout<<"\nPARENT"<<l[0];
		t.remove_edge(k2.idx,l[0]);
		t.add_edge(l[0],root_node);
		t.add_edge(root_node,k2.idx);
		//BOOST_FOREACH(j, t.adjacent(t.root)) 
		//{
		//cout<<"\nROOT ADJACENT AFTER"<<j;
		//}
					   
		// output taxa			       
		aw::tree2newick(os, t, t_name, t_weight);
		os << std::endl;		      
	      }
	    //else cout<<"\nNOt if"<<k2.idx;
			      
	  } 
	  //else {
	  
	  //}cout<<"\nROOTNAME = "<<t_name[k2.idx];
	t.contract_all_chains();   
	t = t2;
	t_name = t_name2;
	t_weight = t_weight2;
				 
      }
			  
 
  }
  ifs.close();
  ifs2.close();
  ofs << os.str()<< std::endl;
  if (!ofs) ERROR_exit("cannot write file '" << filename2 << "'");			      
  ofs.close();
  cout<<"DONE !\n"; 
  return 0;
}
  
