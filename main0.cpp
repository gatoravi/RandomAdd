/* TOTALLY RANDOM ADDITION Version 1*/
// Written by Avinash using Andre Wehe's Library Feb 11 2011
// sample usage - ./exec tree_file leaves_file  
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

#define NEWLEAFNO 1000 /* Max Number of new leaves that can be added */

#define ADDED 1
#define NOTADDED 0
using namespace std;

// -----------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------
int binarysearch(double *larray, int size, double key){
  int first = 0;
  int last = size-1;
  cout<<"\nSize of array "<<size<<" key "<<key<<" last "<<last;
  while (first <= last) {
    int mid = (first + last) / 2;  // compute mid point.
    cout<<"\nlarraymid "<<mid<<" "<<larray[mid]<<"\nlarraymid+1 "<<mid+1<<" "<<larray[mid+1];
    if (key > larray[mid] && key > larray[mid+1]) 
      first = mid + 1;  // repeat search in top half.
    else if (key < larray[mid] && key < larray[mid+1]) 
      last = mid - 1; // repeat search in bottom half.
    else{
      cout<<"\nReturned value: "<<mid+1;
      return mid+1;     // found it. return position 
    }
           
  }
  return 0;
}


int main(int argc, char* argv[]) { 
 
  
  if(argc!=3){
    cout<<"Incorrect Arguments ! Exiting! \n\t usage - ./exec tree_file leaves_file\n";
    exit(1);
  }
 
  
  srand(time(NULL));
  char* treefile = argv[1];
  char* leaves_file = argv[2];
  
  
  /* Read in tree from treefile */
  std::ifstream ifs;
  ifs.open (treefile);
  aw::Tree t;
  aw::idx2name t_name;
  aw::idx2weight_double t_weight;   
  if(ifs.good()){       
    if (!aw::stream2tree(ifs, t, t_name, t_weight)){
      cout<<"Unable to read tree from file! Exiting!";
      exit(1);
    }    	
  }
  else {
    cout<<"Unable to open tree file! Exiting!\n";
    exit(1);
  }  
  
  /* Print out tree */
  ostringstream os;
  aw::tree2newick(os, t, t_name, t_weight);
  cout<<"\nOutput : "<<os.str()<<endl;
  
  
  /* Read in leaves array from file */
  string leaves_array[NEWLEAFNO];
  int leafcount = 0;  
  std::ifstream ifs2(leaves_file);
  string current_leaf;
  if(ifs2.is_open()){
    getline(ifs2 ,current_leaf);
    while(ifs2.good()){  
      leaves_array[leafcount++] = current_leaf;    
      getline(ifs2 ,current_leaf);    
    }
  }
  
  /* Traverse tree and store branch lengths in an array */
  
  int total_nodes = 2*(leafcount + NEWLEAFNO -1);
  double* bl_array = new double[total_nodes];
  int* parent_array = new int[total_nodes];
  int* translate_index = new int[total_nodes];
  double total_blength = 0;
  int leafn = 0;
  int currnodecount = 0;
  
  TREE_POSTORDER2(k,t){
    int currnode = k.idx;
    cout<<"\nWeight for node "<<currnode<<" "<<t_name[currnode];
    cout<<" is "<<t_weight[currnode][0];    
    total_blength += t_weight[currnode][0]; 
    translate_index[currnodecount] = currnode;
    bl_array[currnodecount++] = total_blength;//this array stores cumulative branch lengths & hence is sorted.
    if(t.is_leaf(k.idx)){
      leafn++;
    }
    parent_array[currnode] = k.parent;     
  } 
  cout<<"\n\n";
  for(int i=0; i<currnodecount; i++){
    int translate = translate_index[i]; 
    cout<<"\nbrancharray "<<translate<<" "<<bl_array[i];
  }
  cout<<"\nThe total branch length is "<<total_blength; 
  
  int currleafcount = leafn;  
  int* added_leaf = new int[leafn + leafcount]; 
  int* leafindex  = new int[leafn + leafcount]; 
  for(int j=0; j<leafcount; j++) {
    leafindex[j] = j;
    added_leaf[j] = NOTADDED;
  }
  
  int edgec = t.edge_size();
  int nodec = t.node_size();
  cout<<"\nNumber of edges is "<<edgec;
  cout<<"\nNumber of nodes is "<<nodec;
  int templc = leafcount;
  cout<<"\nNumber of leaves to be added"<<templc;
  
  for(int j=0; j<leafcount; j++) {
    /* Pick a random leaf */
    int random_leaf;
    int rnum;
    do{
      rnum = rand() % templc;
      random_leaf = leafindex[rnum];
      cout<<"\nj"<<j;
      cout<<"\nTEMPLC"<<templc;
      cout<<"\nleafcount"<<leafcount;
      cout<<endl<<random_leaf<<endl<<" "<<added_leaf[random_leaf];
      //cout<<"chk";
    }while(added_leaf[random_leaf]==ADDED);
    added_leaf[random_leaf] = ADDED;
    
    /* reduce length of leafindex by 1 and swap chosen leaf with last element */
    leafindex[rnum] = leafindex[templc-1];
    templc--;    
    string newleaf = leaves_array[random_leaf];
    
    /* Pick a random length to select branch*/
    cout<<"\nTotal blength"<<total_blength;
    double randomblength = (double)rand() * (double)total_blength / (double)RAND_MAX;
    //double randomblength = fmod((double)rand() , (double)total_blength );
    
    /* Find the node to insert leaf based on branch lengths */
    unsigned int untranslated_node = binarysearch( bl_array, currnodecount, randomblength);  
    unsigned int insert_node = translate_index[untranslated_node]; 
    cout<<"\nNode to be inserted is: "<<insert_node; 
    
    
    /* insert new internal node to attach new leaf */
    unsigned int i_n = t.new_node();
    if(insert_node != t.root)
      bl_array[currnodecount] = bl_array[currnodecount-1] + reduce_length;
    translate_index[currnodecount] = i_n;
    currnodecount++;
  
    /* insert the new leaf */
    unsigned int l_n = t.new_node();  
    //bl_array[currnodecount] = bl_array[currnodecount-1] + randomblength;
    
    translate_index[currnodecount] = l_n;
    currnodecount++;
    currleafcount++;
    
    /* Reduce the branch length of the node where you insert by the original length - new random length*/
    double original_length = bl_array[untranslated_node];
    double reduce_length = original_length - randomblength;
    bl_array[untranslated_node] = bl_array[untranslated_node-1] + randomblength;
    
    
    /* insert the branch lengths of two new nodes after the insertion node */
    double temp1 = bl_array[untranslated_node+1];
    bl_array[currnodecount] = bl_array[currnodecount-1] + reduce_length;
    double temp1 = bl_array[untranslated_node+2];
    bl_array[currnodecount] = bl_array[currnodecount-1] + randomblength;
    
    
    
    
    /* Reduce the branch lengths of remaining nodes */
    for(int i= untranslated_node+2; i<currnodecount; i++) {
      bl_array[i] -= reduce_length;
    }
    
    total_blength =  bl_array[currnodecount];
    
    
      
    
  
    cout<<"\nReduce Length "<<reduce_length<<" Randomblength "<<randomblength<<" original length "<<original_length;
    cout<<"\nnew curr node count "<<currnodecount;
  
    unsigned int current_parent; 
    if(insert_node != t.root){
      //cout<<"\nnot root selected";
      current_parent = parent_array[insert_node];
      //cout<<"\nThe parent is : "<<parent;	      
      t.remove_edge(current_parent, insert_node);  
    }  
  
    t_name[l_n] = newleaf;
    t.add_edge(i_n, l_n);
    parent_array[l_n] = i_n;
    t.add_edge(i_n, insert_node);
    parent_array[insert_node] = i_n;
    if(insert_node != t.root){
      t.add_edge(i_n, current_parent);
      parent_array[i_n] = current_parent;      
    }      
    else
      t.root = i_n;
    
    edgec = t.edge_size();
    nodec = t.node_size();
  
    cout<<"\nNumber of edges is "<<edgec;
    cout<<"\nNumber of nodes is "<<nodec;
  
  
    //for(int i=0; i<currnodecount; i++){
      //int translate = translate_index[i]; 
      //cout<<"\nbrancharray "<<translate<<" "<<bl_array[i];
    //}
  
    //cout<<"\n";
    //TREE_POSTORDER2(k,t) {
    //cout<<k.idx<<" ";
    //} 
    /* Print out tree */
    /*ostringstream os;
    aw::tree2newick(os, t, t_name, t_weight);
    cout<<"\nOutput : "<<os.str()<<endl;*/
  } 
  cout<<"\n";
  return 0;
  
  
}  





/*
 
  srand(time(NULL));
  char leaves_file[100];
  int parent_node[LEAFNO*2-1];
  int parent_node_temp[LEAFNO*2-1];
  int inout[LEAFNO*2-1];
  if(argc!=6){
  cout<<"Incorrect Arguments ! Exiting! \n\t Usage - ./exec leaves_file  mast_size treenum mastpresentpc subtreenoise\n";
  exit(1);
  }
  strcpy(leaves_file, argv[1]);
  string current_leaf;
  string leaves_array[LEAFNO];
  int leafcount = 0;  
  std::ifstream ifs(leaves_file);
  if(ifs.is_open()){
  getline(ifs ,current_leaf);
  while(ifs.good()){  
  leaves_array[leafcount++] = current_leaf;    
  getline(ifs ,current_leaf);    
  }
  }
  else{
  cout<<"\nUnable to open taxa file !\n"; 
  exit(1);
  }
  ifs.close();
  string output_file = "";
  output_file.append(argv[2]);
  output_file.append("_");
  output_file.append(argv[4]);
  output_file.append("_");
  output_file.append(argv[5]);
  output_file.append(".txt");
  char opfile[100];
  strcpy(opfile, output_file.c_str());
  std::ofstream fout(opfile);
  cout<<"\nthe number of leaves is "<<leafcount;  
  //cout<<"\ncheck";  
  cout<<"\n";
  aw::Tree MAST;
  aw::idx2name t_name;
  //aw::idx2weight_double t_weight; 
  
  istringstream test_ipstream;  
  string test_tree;
  test_tree = "(";
  test_tree.append(leaves_array[0]);
  test_tree.append(",");
  test_tree.append(leaves_array[1]);
  test_tree.append(");");
  
  test_ipstream.str(test_tree);
   
  int edgec, nodec;
  if (!aw::stream2tree(test_ipstream, MAST, t_name)) cout<<"Unable to read tree  !!";
  edgec = MAST.edge_size();
  nodec = MAST.node_size();
  //cout<<"\nNumber of edges is "<<edgec;
  //cout<<"\nNumber of nodes is "<<nodec;
  
  int mast_lc;
  //cout<<"\nEnter size of MAST: ";
  mast_lc = atoi(argv[2]);
  //cin>>mast_lc;
  
  if(mast_lc<2){
  cout<<"\nInvalid MAST size exiting.\n";
  exit(1);
  }
  
  int mast_l = 2;
  parent_node[0] = 0;
  parent_node[1] = 0;
  parent_node[2] = 0;
  
  for(int i=0; i<LEAFNO*2-1; i++) {
  inout[i] = OUTER;
  }
  
  
  for(int i=0; i<mast_lc-2; i++){
  string newleaf = leaves_array[mast_l];
  //cout<<"\n\nNew leaf is "<<newleaf;
  int nodes = 2*(mast_l++) - 1;   
  unsigned int node_sel = rand()%nodes;  
  //cout<<"\nthe selected node is "<<node_sel;
  int parent; 
  if(node_sel != MAST.root){
  //cout<<"\nnot root selected";
  parent = parent_node[node_sel];
  //cout<<"\nThe parent is : "<<parent;	      
  MAST.remove_edge(parent, node_sel);  
  }  
  //else 
  //cout<<"\n root selected";
  int i_n = MAST.new_node();
  int l_n = MAST.new_node();
  t_name[l_n] = newleaf;
  MAST.add_edge(i_n, l_n);
  parent_node[l_n] = i_n;
  MAST.add_edge(i_n, node_sel);
  parent_node[node_sel] = i_n;
  if(node_sel != MAST.root){
  MAST.add_edge(i_n, parent);
  parent_node[i_n] = parent;      
  }      
  else
  MAST.root = i_n;
  edgec = MAST.edge_size();
  nodec = MAST.node_size();
  //cout<<"\nNumber of edges is "<<edgec;
  //cout<<"\nNumber of nodes is "<<nodec; 
    
  //cout<<"\n";
  //TREE_POSTORDER2(k,MAST) {
  //cout<<k.idx<<" ";
  //}  
     
  }
  ostringstream os;
  aw::tree2newick(os, MAST, t_name);
  //cout<<"\nMAST : "<<os.str()<<endl;
  fout<<os.str()<<endl;
  int treenum;
  int MAST_presentpc;
  //cout<<"\nEnter the total number of trees to be created: ";
  //cin>>treenum;  
  treenum = atoi(argv[3]);
  //cout<<"\nEnter MAST present percentage (from 0 to 100 %) : ";
  //cin>>MAST_presentpc;
  MAST_presentpc = atoi(argv[4]);
  //cout<<"\nEnter the subtree noise (from 0 to 100 %)";
  int subtreen;
  subtreen = atoi(argv[5]);
  //cin>>subtreen;
  
  int present_trees = (MAST_presentpc * treenum) / 100;
  int absent_trees = treenum - present_trees;
*/
  
/* Store the parents of the MAST subtree nodes */
/*  for(int i=0; i<2*mast_lc-1; i++) {
    parent_node_temp[i] = parent_node[i];
    }*/
  
/* Trees Containing the MAST */  
/* for(int i=0; i<present_trees; i++){  
//cout<<"\n";
for(int i=0; i<2*mast_lc-1; i++) {
inout[i] = INNER;
}
inout[MAST.root] = OUTER;
ostringstream os;
aw::Tree tree_i = MAST;
aw::idx2name name_i = t_name;    
for(int j=mast_lc; j<leafcount; j++){      
string newleaf = leaves_array[j];
//cout<<"\n\nNew leaf is "<<newleaf;
int nodes = 2*(j) - 1;   
unsigned int node_sel;
int random_n = rand()%100;
      
//if number is between 0 and subtreeprob then put it out.
if(random_n<subtreen){
do{
node_sel = rand()%nodes; 
}while(inout[node_sel] != OUTER);
//cout<<"OUTER";
}
else{
//if number is between  subtreeprob and 100 then put it in.
do{
node_sel = rand()%nodes; 
}while(inout[node_sel] != INNER);
//cout<<"INNER";
}
//cout<<"\nthe selected node is "<<node_sel;
int parent; 
if(node_sel != tree_i.root){
//cout<<"\nnot root selected";
parent = parent_node[node_sel];
//cout<<"\nThe parent is : "<<parent;	      
tree_i.remove_edge(parent, node_sel);  
}  
//else 
//cout<<"\n root selected";
int i_n = tree_i.new_node();
int l_n = tree_i.new_node();
if(inout[node_sel] == INNER){
inout[i_n] = INNER;
inout[l_n] = INNER;
}
else{
inout[i_n] = OUTER;
inout[l_n] = OUTER;      
}      
name_i[l_n] = newleaf;
tree_i.add_edge(i_n, l_n);
parent_node[l_n] = i_n;
tree_i.add_edge(i_n, node_sel);
parent_node[node_sel] = i_n;
if(node_sel != tree_i.root){
tree_i.add_edge(i_n, parent);
parent_node[i_n] = parent;      
}      
else
tree_i.root = i_n;
//edgec = tree_i.edge_size();
//nodec = tree_i.node_size();
//cout<<"\nNumber of edges is "<<edgec;
//cout<<"\nNumber of nodes is "<<nodec;     
cout<<"\n";
TREE_POSTORDER2(k,tree_i) {
cout<<k.idx<<" ";
}*/     
/*  }  
    aw::tree2newick(os, tree_i, name_i);
    //cout<<"\nOutput : "<<os.str(); 
    fout<<os.str();
    fout<<endl;
    Now copy the parents of the MAST from the temp array to the main array*/
/*  for(int i=0; i<2*mast_lc-1; i++){
    parent_node[i] = parent_node_temp[i];
    }
    }
  
    aw::Tree MAST2;
    istringstream test_ipstream2;  
    test_tree = "(";
    test_tree.append(leaves_array[0]);
    test_tree.append(",");
    test_tree.append(leaves_array[1]);
    test_tree.append(");");  
  
    //cout<<"\ntest_tree"<<test_tree;
    test_ipstream2.str(test_tree);  
    if (!aw::stream2tree(test_ipstream2, MAST2, t_name)) cout<<"Unable to read tree  !!";
    parent_node[0] = 0;
    parent_node[1] = 0;
    parent_node[2] = 0;
  
    //cout<<"\nPresent trees: "<<present_trees;
    //cout<<"\nAbsent  trees: "<<absent_trees;
  
    Totally Random Trees for the remaining trees */
/*  for(int i=0; i<absent_trees; i++){    
    ostringstream os;
    aw::Tree tree_i = MAST2;
    aw::idx2name name_i = t_name;    
    for(int j=2; j<leafcount; j++){
    string newleaf = leaves_array[j];
    //cout<<"\n\nNew leaf is "<<newleaf;
    int nodes = 2*(j) - 1;   
    unsigned int node_sel = rand()%nodes;  
    //cout<<"\nthe selected node is "<<node_sel;
    int parent; 
    if(node_sel != tree_i.root){
    //cout<<"\nnot root selected";
    //cout<<"\nj "<<j<<"node sel"<<node_sel;
    parent = parent_node[node_sel];
    //cout<<"\nThe parent is : "<<parent;	      
    tree_i.remove_edge(parent, node_sel);  
    }  
    //else 
    //cout<<"\n root selected";
    int i_n = tree_i.new_node();
    int l_n = tree_i.new_node();
    name_i[l_n] = newleaf;
    tree_i.add_edge(i_n, l_n);
    parent_node[l_n] = i_n;
    tree_i.add_edge(i_n, node_sel);
    parent_node[node_sel] = i_n;
    if(node_sel != tree_i.root){
    tree_i.add_edge(i_n, parent);
    parent_node[i_n] = parent;      
    }      
    else
    tree_i.root = i_n;
    edgec = tree_i.edge_size();
    nodec = tree_i.node_size();
    //cout<<"\nNumber of edges is "<<edgec;
    //cout<<"\nNumber of nodes is "<<nodec; 
    
    cout<<"\n";
    TREE_POSTORDER2(k,tree_i) {
    cout<<k.idx<<" ";
    }*/     
/*   }  
     aw::tree2newick(os, tree_i, name_i);
     //cout<<"\nOutput : "<<os.str(); 
     fout<<os.str();
     fout<<endl;
     parent_node[0] = 0;
     parent_node[1] = 0;
     parent_node[2] = 0;
     }
  
     //aw::tree2newick(os, MAST, t_name, t_weight);
     //cout<<"\nOutput : "<<os.str();
     cout<<"\n";
     return 0;
     }










  
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
     }*/



// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
//MSG("2");//avi_test
// read trees
/*std::ifstream ifs,ifs2;
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
  ofs << os.str() << std::endl;
  if (!ofs) ERROR_exit("cannot write file '" << filename2 << "'");			      
  ofs.close();
  cout<<"DONE !\n"; 
  return 0;
  }*/
  
