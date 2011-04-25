/* TOTALLY RANDOM ADDITION Version 1
*  Written by Avinash using Andre Wehe's Library Feb 11 2011
*  usage - ./executable tree_file leaves_file replicates opfile 
*  Program to add new taxa to an existing phylogeny based on branch lengths of the initial tree
*/

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

#define MAXADDLEAF 2000 /* Max Number of new leaves that can be added */
#define ADDED 1
#define NOTADDED 0

using namespace std;

// -----------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------
int binarysearch(double *larray, int size, double key){
  int first = 0;
  int last = size-1;
  //cout<<"\nSize of array "<<size<<" key "<<key<<" last "<<last;
  while (first <= last) {
    int mid = (first + last) / 2;  // compute mid point.
    //cout<<"\nlarraymid "<<mid<<" "<<larray[mid]<<"\nlarraymid+1 "<<mid+1<<" "<<larray[mid+1];
    if (key > larray[mid] && key > larray[mid+1]) 
      first = mid + 1;  // repeat search in top half.
    else if (key < larray[mid] && key < larray[mid+1]) 
      last = mid - 1; // repeat search in bottom half.
    else{
      //cout<<"\nReturned value: "<<mid+1;
      return mid+1;     // found it. return position 
    }
           
  }
  return 0;
}


int main(int argc, char* argv[]) {  
  
  if(argc!=5){
    cout<<"Incorrect Arguments ! Exiting! \n\t usage - ./exec tree_file leaves_file replicates opfile\n";
    exit(1);
  }
 
  /* set seed of Rand number generator */
  srand(time(NULL));
  
  /* Read in command line arguments */
  char* treefile = argv[1];
  char* leaves_file = argv[2];  
  int replicates = atoi(argv[3]);
  char* opfile = argv[4];
  
  /* Read initial tree from treefile */
  std::ifstream ifs;
  ifs.open (treefile);
  aw::Tree initial_t;
  aw::idx2name initial_t_name;
  aw::idx2weight_double initial_t_weight;   
  if(ifs.good()){       
    if (!aw::stream2tree(ifs, initial_t, initial_t_name, initial_t_weight)){
      cout<<"Unable to read tree from file! Exiting!";
      exit(1);
    }    	
  }
  else {
    cout<<"Unable to open tree file! Exiting!\n";
    exit(1);
  }  
  
  
  /* Open output file */
  ofstream ofs(opfile);
  if(!ofs.good()) {
    cout<<"unable to open output file!";
    exit(1);
  }   
   
  /* Read in leaves to be added from file */
  string leaves_array[MAXADDLEAF];
  int leafcount = 0;  //number of leaves to be added.
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
  int total_nodes = 2*(leafcount + MAXADDLEAF -1);
  double* bl_array_initial = new double[total_nodes];
  int* parent_array_initial = new int[total_nodes];
  unsigned int* translate_index_initial = new unsigned int[total_nodes];
  double total_blength_initial = 0;
  int leafn = 0;
  int initial_currnodecount = 0;
  
  TREE_POSTORDER2(k,initial_t){
    unsigned int currnode = k.idx;
    //cout<<"\nWeight for node "<<currnode<<" "<<initial_t_name[currnode];
    //cout<<" is "<<initial_t_weight[currnode][0];    
    total_blength_initial += initial_t_weight[currnode][0]; 
    translate_index_initial[initial_currnodecount] = currnode;
    if(currnode != initial_t.root)
      bl_array_initial[initial_currnodecount] = total_blength_initial;//this array stores cumulative branch lengths & hence is sorted.
    initial_currnodecount++;
    parent_array_initial[currnode] = k.parent;     
  } 
  
  /*cout<<"\n\nBEFORE ADDITION\n\n";
  for(int i=0; i<initial_currnodecount-1; i++){
    int translate = translate_index_initial[i]; 
    cout<<"\nbranch number "<<translate;
    if(i!=0)cout<<" length "<<bl_array_initial[i] - bl_array_initial[i-1];
    else cout<<" length "<<bl_array_initial[i];
    cout<<" cumul length "<<bl_array_initial[i];
  }*/
     
  //int currleafcount = leafn;  
  int* added_leaf = new int[leafn + leafcount]; 
  int* leafindex  = new int[leafn + leafcount]; 
  
  
  /* CREATE REPLICATES OF NEW  TREE */
  for(int k=0; k<replicates; k++) {
    
    /* Copy the initial node count and branch length */
    int currnodecount = initial_currnodecount;
    //cout<<"Initial curr node count is: "<<currnodecount;
    double total_blength = total_blength_initial;
    cout<<"\n\nReplicate number : "<<k+1;
    cout<<"\nThe total branch length is "<<total_blength;
    
    /* Copy the initial branch lengths and node parents to new arrays */
    double* bl_array = new double[total_nodes];
    int* parent_array = new int[total_nodes];
    unsigned int* translate_index = new unsigned int[total_nodes];
    for(int i=0; i<currnodecount; i++) {
      translate_index[i] = translate_index_initial[i];
      parent_array[i]    = parent_array_initial[i];
      bl_array[i]        = bl_array_initial[i];
    }
  
    /* Declare tree, labels & weights */
    aw::Tree t = initial_t;
    aw::idx2name t_name = initial_t_name;
    aw::idx2weight_double t_weight = initial_t_weight;
       
    /* Initialise all leaves as not added and copy leaf indices */
    for(int j=0; j<leafcount; j++) {
      leafindex[j] = j;
      added_leaf[j] = NOTADDED;
    }
  
    //int edgec = t.edge_size();
    //int nodec = t.node_size();
    //cout<<"\nBEfore addition Number of edges is "<<edgec;
    //cout<<"\nBefore addition Number of nodes is "<<nodec;
    
    /* Temporary count of leaves remaining to be added */
    int templc = leafcount;
    //cout<<"\nNumber of leaves to be added "<<templc;
  
    for(int j=0; j<leafcount; j++) {       
      /* Check the initial number of nodes and edges */
      int edgec = t.edge_size();
      int nodec = t.node_size();
      
      /* Pick a random leaf */
      int random_leaf;
      int rnum;
      do {
	rnum = rand() % templc;
	random_leaf = leafindex[rnum];      
      }while(added_leaf[random_leaf]==ADDED);   
      added_leaf[random_leaf] = ADDED;
      //cout<<"\nLeaf Number "<<random_leaf+1; // leaf numbers 1 to n
     
      /* Reduce length of leafindex by 1 and swap chosen leaf with last element */
      leafindex[rnum] = leafindex[templc-1];
      templc--;     
      string newleaf = leaves_array[random_leaf];    
    
    
      /* Pick a random length to select branch */
      double randomblength;
          
      /* Find the node to insert leaf based on branch lengths */
      unsigned int untranslated_node;
      do {
        randomblength = (double)rand() * (double)total_blength / (double)RAND_MAX; 
        //cout<<"\nTemp total branch length"<<total_blength;
        //cout<<" temp random branch length"<<randomblength;
        untranslated_node = binarysearch( bl_array, currnodecount-1, randomblength);  
      } while(translate_index[untranslated_node] == t.root);
      
      
      /* obtain the individual lengths by subtracting the cumulative lengths */
      double original_length = bl_array[untranslated_node] - bl_array[untranslated_node-1];
      randomblength = randomblength - bl_array[untranslated_node-1];
    
      double reduce_length =   original_length - randomblength; 
      unsigned int insert_branch = translate_index[untranslated_node]; 
    
      /* change length of branch where inserted */
      t_weight[insert_branch][0] = randomblength;
   
      /* insert new internal node to attach new leaf */
      unsigned int i_n = t.new_node();
      t_weight[i_n][0] = reduce_length;    
  
      /* insert the new leaf */
      unsigned int l_n = t.new_node(); 
      t_weight[l_n][0] = randomblength; 
      
      /* remove the edge b/w selected node and parent */
      unsigned int current_parent = parent_array[insert_branch];
      //cout<<"\nCurrent parent "<<current_parent<<" insert branch "<<insert_branch;
      t.remove_edge(current_parent, insert_branch);       
      
      //cout<<"\nnot root selected";
      //cout<<"\nThe parent is : "<<parent;	      
     
      /* add an edge b/w new leaf and new internal */ 
      t_name[l_n] = newleaf;
      t.add_edge(i_n, l_n);
      parent_array[l_n] = i_n;
    
      /* add an edge b/w new internal and selected node */
      t.add_edge(i_n, insert_branch);
      parent_array[insert_branch] = i_n;
      t.add_edge(i_n, current_parent);
    
      /* assign parent of new internal to old parent of selected */
      parent_array[i_n] = current_parent; 
      edgec = t.edge_size();
      nodec = t.node_size();
      
      /* Recalculate total branch length and number of nodes */
      total_blength = 0;
      currnodecount = 0;    
      TREE_POSTORDER2(k,t){
	unsigned int currnode = k.idx;
	//cout<<"\nWeight for node "<<currnode<<" "<<t_name[currnode];
	//cout<<" is "<< t_weight[currnode][0];    
	total_blength += t_weight[currnode][0]; 
	translate_index[currnodecount] = currnode;
	if(currnode!=t.root)
	  bl_array[currnodecount++] = total_blength;//this array stores cumulative branch lengths & hence is sorted.
	parent_array[currnode] = k.parent;     
      }   
      
      //int wait;
      //cin>>wait;
    } 
    delete [] parent_array;
    delete [] bl_array;
    delete [] translate_index;
    aw::tree2newick(ofs, t, t_name, t_weight);
    ofs<<endl;
    cout<<"\nThe new total branch length is: "<<total_blength;
  } 
  cout<<"\n";
  
  /* deallocate memory */
  delete [] added_leaf; 
  delete [] leafindex; 
  delete [] bl_array_initial; 
  delete [] parent_array_initial; 
  delete [] translate_index_initial; 
  
  return 0;  
}  
