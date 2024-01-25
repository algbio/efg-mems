#include <iostream>
#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include "br_index.hpp"
#include "utils.hpp"
#include "bdbwt2brindex.hpp"

using namespace bri;
using namespace std;
using namespace sdsl;

//#define EFG_MEMS_DEBUG

// Global variables assigned with optional parameters
ulint kappa = 1; // MEM threshold
bool asymmetric = false;
bool use_brindex = true; // if false, then use bdbwt
string alphabet = "";
string output_file = "";
bool fasta = false; // reading queries as fasta
bool msa_output = false; // map MEMs back to the MSA, instead of the graph

struct starting_positions {
	bit_vector nodes;
	bit_vector edges;
	rank_support_v5<> nodes_rs;
	rank_support_v5<> edges_rs;
	select_support_mcl<> nodes_ss;
	select_support_mcl<> edges_ss;
};

struct gfa_info {
	vector<string> ordered_node_ids;
	vector<int> node_lengths;
	vector<pair<int,int>> ordered_edges;
};

struct query_info {
	vector<string> ordered_query_ids;
	vector<ulint> ordered_query_lengths;
	bit_vector concat;
	rank_support_v5<> concat_rs;
	select_support_mcl<> concat_ss;
};


void help()
{
	cout << "efg-mems: locate all MEMs between queries and elastic founder graph" << endl;

	cout << "Usage: efg-mems [options] <queries> <efg>" << endl;
        cout << "   --asymmetric  use asymmetric definition for MEMs " << endl;
	cout << "   --bdbwt       use bdbwt instead of default br-index " << endl;
	cout << "   -k      MEM threshold" << endl;
	cout << "   -a      alphabet string with last three symbols regarded special, default ACGTN#0" << endl;	
	cout << "   -o      output file where lines x,i,d are outputed for MEMs Q[i..i+d-1]=T[x..x+d-1]; " << endl;
    cout << "           in symmetric mode (default) the matches are locally maximal, i.e., Q[i-1]!=T[x-1] and Q[i+d]!=T[x+d] " << endl;
	cout << "           in asymmetric mode the matches are globally maximal, i.e., Q[i..i+d] or Q[i-1..i+d-1] do not " << endl;
	cout << "           occur in T; only one occurrence T[x..x+d-1] is reported in asymmetric mode" << endl;		
	cout << "           output is formatted as fasta file with MEMs listed according to the input efg fasta" << endl;		
	cout << "   -f     file containing alphabet " << endl;
	cout << "   <queries>  index file of concatenation of queries (with extension .bri), or plain text or fasta if using --bdbwt" << endl;
	cout << "               index or plain text: concatenation format: #AGGATG#AGATGT#, where # is separator symbol" << endl;
        cout << "               fasta: alternating header >H and query lines Q; query is converted to #Q#, so MEM coordinates will be +1" << endl;
	cout << "   <efg>      elastic founder graph in GFA format (with extension .gfa)" << endl;
	exit(0);
}

bool parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;


    if (s.compare("--asymmetric") == 0)
    {

        asymmetric = true;

    }
    else if (s.compare("--bdbwt") == 0)
    {

        use_brindex = false;

    }    
    else if (s.compare("-k") == 0) 
    {

		if(ptr>=argc-1){
			cerr << "Error: missing parameter after -k option." << endl;
			help();
		}

		kappa = atoi(argv[ptr]);
		ptr++;
    
    }
    else if (s.compare("-a") == 0)
    {

		if(ptr>=argc-1){
			cerr << "Error: missing parameter after -a option." << endl;
			help();
		}
        alphabet = string(argv[ptr]);
        ptr++;
    }
    else if (s.compare("-o") == 0)
    {

		if(ptr>=argc-1){
			cerr << "Error: missing parameter after -o option." << endl;
			help();
		}
        output_file = string(argv[ptr]);  
        ptr++;
    }
    else if (s.compare("-f") == 0)
    {

		if(ptr>=argc-1){
			cerr << "Error: missing parameter after -f option." << endl;
			help();
		}
        string alphabet_file="";
        alphabet_file = string(argv[ptr]);
        ifstream as(alphabet_file);
        getline(as,alphabet);
        as.close();    
    }
    else if (s.compare("--map-to-msa") == 0)
    {
	    msa_output = true;
    }
    else
    {
        return 0;
	}
	return 1;

}

void read_gfa(ifstream& gfa, string& nodes, string& edges) {

    vector<string> node_labels;
    vector<string> edge_labels;
    vector<int> edge_node_start;
    vector<int> edge_node_end;    
    unordered_map<int, unordered_set<char>> lext; // chars to the left from nodes 
    unordered_map<int, unordered_set<char>> rext;  // chars to the right from nodes
    
    string lextnodes, rextnodes, lextedges, rextedges;
    
    string line;
    while (getline(gfa, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty or comment lines
        }

        vector<string> fields;
        string field;
        for (char c : line) {
            if (c == '\t') {
                fields.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        fields.push_back(field);

        if (fields[0] == "S") {
            node_labels.push_back(fields[2]);
        } else if (fields[0] == "L") {
            int start_node_id = stoi(fields[1]);
            bool start_node_orientation = (fields[2] == "-");
            int end_node_id = stoi(fields[3]);
            bool end_node_orientation = (fields[4] == "-");
            lext[end_node_id].insert(node_labels[start_node_id-1][node_labels[start_node_id-1].size()-1]);
            rext[start_node_id].insert(node_labels[end_node_id-1][0]);            
            string start_node_label = node_labels[start_node_id - 1];
            if (start_node_orientation) {
                start_node_label = start_node_label + "_rev";
            }
            string end_node_label = node_labels[end_node_id - 1];
            if (end_node_orientation) {
                end_node_label = end_node_label + "_rev";
            }
            edge_labels.push_back(start_node_label + ">" +end_node_label);
            edge_node_start.push_back(start_node_id);
            edge_node_end.push_back(end_node_id);
	}
    }

    nodes = alphabet[alphabet.size()-1];
    lextnodes = alphabet[alphabet.size()-1];
    rextnodes = alphabet[alphabet.size()-1];
    for (int i=0; i < node_labels.size(); i++) {
       nodes += alphabet[alphabet.size()-1]+node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // left and right chars not yet known
       if (i==0) 
          lextnodes += alphabet[alphabet.size()-3] + node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       else if (lext[i].size()==1) { 
          for (char lextchar : lext[i])
             lextnodes += lextchar + node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // unique left
       }
       else
          lextnodes += alphabet[alphabet.size()-3] + node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       if (rext[i].size()==1) {
          for (char rextchar : rext[i])
             rextnodes += alphabet[alphabet.size()-1]  + node_labels[i] + rextchar + alphabet[alphabet.size()-1]; // unique right
       }
       else
          rextnodes += alphabet[alphabet.size()-1] + node_labels[i] + alphabet[alphabet.size()-3] + alphabet[alphabet.size()-1]; //ambiguous right
    }
  
    edges = alphabet[alphabet.size()-1];
    lextedges = alphabet[alphabet.size()-1];
    rextedges = alphabet[alphabet.size()-1];
    for (int i=0; i < edge_labels.size(); i++) {
       edges += alphabet[alphabet.size()-1]+edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // left and right chars not yet known
       if (i==0) 
          lextedges += alphabet[alphabet.size()-3] + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       else if (lext[edge_node_start[i]].size()==1) { 
          for (char lextchar : lext[edge_node_start[i]])
             lextedges += lextchar + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // unique left
       }
       else
          lextedges += alphabet[alphabet.size()-3] + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       if (rext[edge_node_end[i]].size()==1) {
          for (char rextchar : rext[edge_node_end[i]])
             rextedges += alphabet[alphabet.size()-1]  + edge_labels[i] + rextchar + alphabet[alphabet.size()-1]; // unique right
       }
       else
          rextedges += alphabet[alphabet.size()-1] + edge_labels[i] + alphabet[alphabet.size()-3] + alphabet[alphabet.size()-1]; //ambiguous right
    }
   // combining left and right chars to nodes and edges
   nodes[1] = lextnodes[1];
   nodes[nodes.size()-2] = rextnodes[nodes.size()-2];
   for (int i=2; i<nodes.size()-2; i++)
      if (nodes[i+1]==alphabet[alphabet.size()-1]) {
         nodes[i] = rextnodes[i];
         nodes[i+2] = lextnodes[i+2];
      }
   edges[1] = lextedges[1];
   edges[edges.size()-2] = rextedges[edges.size()-2];
   for (int i=2; i<edges.size()-2; i++)
      if (edges[i+1]==alphabet[alphabet.size()-1]) {
         edges[i] = rextedges[i];
         edges[i+2] = lextedges[i+2];
      }
}

void read_gfa_generic(ifstream& gfa, string& nodes, string& edges, starting_positions& sp, gfa_info &efginfo) {

    vector<string> ordered_node_ids;
    vector<string> ordered_node_labels;
    std::unordered_map<std::string,int> node_indexes;
    vector<bool> is_source;
    vector<bool> is_sink;
    vector<pair<int,int>> ordered_edges;


    vector<string> edge_labels;
    vector<int> edge_node_start;
    vector<int> edge_node_end;    
    unordered_map<int, unordered_set<char>> lext; // chars to the left from nodes 
    unordered_map<int, unordered_set<char>> rext;  // chars to the right from nodes
    
    string lextnodes, rextnodes, lextedges, rextedges;
    
    string line;
    while (getline(gfa, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty or comment lines
        }

        vector<string> fields;
        string field;
        for (char c : line) {
            if (c == '\t') {
                fields.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        fields.push_back(field);

        if (fields[0] == "S") {
            ordered_node_ids.push_back(fields[1]);
            ordered_node_labels.push_back(fields[2]);
            is_source.push_back(true);
	    is_sink.push_back(true);
	    node_indexes.insert({ fields[1], ordered_node_ids.size() - 1 });
        } else if (fields[0] == "L") {
            //TODO: remove assumption that GFA file is ordered or warn the user
            assert(node_indexes.find(fields[1]) != node_indexes.end()); // maybe use at operator?
            int start_node_id = node_indexes[fields[1]];
            assert(node_indexes.find(fields[3]) != node_indexes.end()); // maybe use at operator?
            int end_node_id = node_indexes[fields[3]];
            lext[end_node_id].insert(ordered_node_labels[start_node_id][ordered_node_labels[start_node_id].size()-1]);
            rext[start_node_id].insert(ordered_node_labels[end_node_id][0]);            
	    ordered_edges.push_back(pair<int,int>(start_node_id, end_node_id));
	    is_source[end_node_id] = false;
	    is_sink[start_node_id] = false;
        }
    }

    efginfo.node_lengths.resize(ordered_node_labels.size());
    for (int i = 0; i < ordered_node_labels.size(); i++) {
	    efginfo.node_lengths[i] = ordered_node_labels[i].length();
	    if (is_source[i])
		    ordered_node_labels[i] = alphabet[alphabet.size()-3] + ordered_node_labels[i];
	    if (is_sink[i])
		    ordered_node_labels[i] = ordered_node_labels[i] + alphabet[alphabet.size()-3];
    }

    {
	    int cumul_node_length = 0;
	    for (int i = 0; i < ordered_node_labels.size(); i++) {
		    cumul_node_length += ordered_node_labels[i].length();
	    }
	    sp.nodes = bit_vector(cumul_node_length + 3*ordered_node_labels.size() + 1, 0);

	    long int d = 2;
	    for (int i = 0; i < ordered_node_labels.size(); i++) {
		    if (is_source[i]) {
			    d += 1;
			    sp.nodes[d] = 1;
			    d += ordered_node_labels[i].length() - 1;
		    } else {
			    sp.nodes[d] = 1;
			    d += ordered_node_labels[i].length();
		    }
		    d += 3;
	    }
    }
    {
	    int cumul_edge_length = 0;
	    for (int i = 0; i < ordered_edges.size(); i++) {
		    int start_node_length = ordered_node_labels[get<0>(ordered_edges[i])].length();
		    int end_node_length = ordered_node_labels[get<1>(ordered_edges[i])].length();
		    cumul_edge_length += start_node_length + end_node_length;
	    }
	    sp.edges = bit_vector(cumul_edge_length + 3*ordered_edges.size() + 1, 0);

	    long int d = 2;
	    for (int i = 0; i < ordered_edges.size(); i++) {
		    int start_node_length = ordered_node_labels[get<0>(ordered_edges[i])].length();
		    int end_node_length = ordered_node_labels[get<1>(ordered_edges[i])].length();

		    if (is_source[get<0>(ordered_edges[i])]) {
			    d += 1;
			    sp.edges[d] = 1;
			    d += start_node_length + end_node_length - 1;
		    } else {
			    sp.edges[d] = 1;
			    d += start_node_length + end_node_length;
		    }
		    d += 3;
	    }

    }

    for (int i = 0; i < ordered_edges.size(); i++) {
            int start_node_id = get<0>(ordered_edges[i]);
            string start_node_label = ordered_node_labels[start_node_id];
            int end_node_id = get<1>(ordered_edges[i]);
            string end_node_label = ordered_node_labels[end_node_id];
            edge_labels.push_back(start_node_label + ">" +end_node_label);
            edge_node_start.push_back(start_node_id);
            edge_node_end.push_back(end_node_id);
    }

    nodes = alphabet[alphabet.size()-1];
    lextnodes = alphabet[alphabet.size()-1];
    rextnodes = alphabet[alphabet.size()-1];
    for (int i=0; i < ordered_node_labels.size(); i++) {
       nodes += alphabet[alphabet.size()-1]+ordered_node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // left and right chars not yet known
       if (i==0) 
          lextnodes += alphabet[alphabet.size()-3] + ordered_node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       else if (lext[i].size()==1) { 
          for (char lextchar : lext[i])
             lextnodes += lextchar + ordered_node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // unique left
       }
       else
          lextnodes += alphabet[alphabet.size()-3] + ordered_node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       if (rext[i].size()==1) {
          for (char rextchar : rext[i])
             rextnodes += alphabet[alphabet.size()-1]  + ordered_node_labels[i] + rextchar + alphabet[alphabet.size()-1]; // unique right
       }
       else
          rextnodes += alphabet[alphabet.size()-1] + ordered_node_labels[i] + alphabet[alphabet.size()-3] + alphabet[alphabet.size()-1]; //ambiguous right
    }
  
    edges = alphabet[alphabet.size()-1];
    lextedges = alphabet[alphabet.size()-1];
    rextedges = alphabet[alphabet.size()-1];
    for (int i=0; i < edge_labels.size(); i++) {
       edges += alphabet[alphabet.size()-1]+edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // left and right chars not yet known
       if (i==0) 
          lextedges += alphabet[alphabet.size()-3] + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       else if (lext[edge_node_start[i]].size()==1) { 
          for (char lextchar : lext[edge_node_start[i]])
             lextedges += lextchar + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // unique left
       }
       else
          lextedges += alphabet[alphabet.size()-3] + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       if (rext[edge_node_end[i]].size()==1) {
          for (char rextchar : rext[edge_node_end[i]])
             rextedges += alphabet[alphabet.size()-1]  + edge_labels[i] + rextchar + alphabet[alphabet.size()-1]; // unique right
       }
       else
          rextedges += alphabet[alphabet.size()-1] + edge_labels[i] + alphabet[alphabet.size()-3] + alphabet[alphabet.size()-1]; //ambiguous right
    }
   // combining left and right chars to nodes and edges
   nodes[1] = lextnodes[1];
   nodes[nodes.size()-2] = rextnodes[nodes.size()-2];
   for (int i=2; i<nodes.size()-2; i++)
      if (nodes[i+1]==alphabet[alphabet.size()-1]) {
         nodes[i] = rextnodes[i];
         nodes[i+2] = lextnodes[i+2];
      }
   edges[1] = lextedges[1];
   edges[edges.size()-2] = rextedges[edges.size()-2];
   for (int i=2; i<edges.size()-2; i++)
      if (edges[i+1]==alphabet[alphabet.size()-1]) {
         edges[i] = rextedges[i];
         edges[i+2] = lextedges[i+2];
      }

   std::swap(ordered_node_ids, efginfo.ordered_node_ids);
   std::swap(ordered_edges, efginfo.ordered_edges);
}

void output_node_MEM_gaf(ulint ti, ulint qi, ulint length, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo)
{
	ulint node = sp.nodes_rs(ti + 1) - 1;
	ulint nodeindex = ti - sp.nodes_ss(node+1);
	ulint query = qinfo.concat_rs(qi + 1) - 1;
	ulint queryindex = qi - qinfo.concat_ss(query + 1);

	// GAF format
	output <<
		qinfo.ordered_query_ids[query] << "\t" <<     // query name
		qinfo.ordered_query_lengths[query] << "\t" << // query length
		queryindex << "\t" <<          // query start (0-based, closed)
		queryindex + length << "\t" << // query end (open)
		"+" << "\t";                  // strand

	output << ">" << efginfo.ordered_node_ids[node] << "\t";        // path
	output << efginfo.node_lengths[node] << "\t"; // path length
	output << nodeindex << "\t";                            // start position on the path
	output << nodeindex + length << "\t";                   // end position on the path
	output << 0 << "\t" << 0 << "\t" << 255 << std::endl; // FIXME
}
void output_edge_MEM_gaf(ulint ti, ulint qi, ulint length, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo)
{
	ulint edge = sp.edges_rs(ti + 1) - 1;
	ulint source = efginfo.ordered_edges[edge].first;
	ulint target = efginfo.ordered_edges[edge].second;
	ulint slength = efginfo.node_lengths[source];
	ulint tlength = efginfo.node_lengths[target];
	ulint edgeindex = ti - sp.edges_ss(edge+1);
	ulint query = qinfo.concat_rs(qi + 1) - 1;
	ulint queryindex = qi - qinfo.concat_ss(query + 1);

	// GAF format
	output <<
		qinfo.ordered_query_ids[query] << "\t" <<     // query name
		qinfo.ordered_query_lengths[query] << "\t" << // query length
		queryindex << "\t" <<          // query start (0-based, closed)
		queryindex + length << "\t" << // query end (open)
		"+" << "\t";                  // strand

	output << ">" << efginfo.ordered_node_ids[source] << ">" << efginfo.ordered_node_ids[target] << "\t";
	output << slength + tlength << "\t";
	output << edgeindex << "\t";                            // start position on the path
	output << edgeindex + length << "\t";                   // end position on the path
	output << 0 << "\t" << 0 << "\t" << 255 << std::endl; // FIXME
}

template <class T, class TS>
bool is_right_maximal(T idx, TS sample) 
{
   uchar c = idx.bwt_at(sample.rangeR.first,true);
   TS right = idx.right_extension(c,sample);
   if (sample.size() == right.size())
      return 0;
   else
      return 1;
}
    
template <class T, class TS>    
bool is_left_maximal(T idx, TS sample) 
{
   uchar c = idx.bwt_at(sample.rangeR.first,false);
   TS left = idx.left_extension(c,sample);
   if (sample.size() == left.size())
      return 0;
   else
      return 1;
}

template<class T, class TS>
ulint reportMEMs(T idx, T qidx, TS sample, TS qsample, ulint d, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo, sdsl::int_vector<> d_bwt = sdsl::int_vector<>(), sdsl::rmq_succinct_sada<> rmq = sdsl::rmq_succinct_sada<>(), ulint* saidx = NULL)
{ 
   ulint MEMcount = 0;
   // a bit naive implementation of the cross product 
   // additive factor alphabet.size()^2(t_biBWTstep+\alphabet.size()^2) slower than an optimal implementation
   std::vector<ulint>** a= new  std::vector<ulint>*[alphabet.size()];
   std::vector<ulint>** b= new  std::vector<ulint>*[alphabet.size()];
   TS** Sa= new  TS*[alphabet.size()];
   TS** Sb= new  TS*[alphabet.size()];
   bool** Ia= new  bool*[alphabet.size()];
   bool** Ib= new  bool*[alphabet.size()];

   TS left,right;
   for (ulint i=0; i<alphabet.size()-1; i++) { // not branching on graph concat separator char, as no valid match can start/end at the boundaries
      a[i] = new std::vector<ulint>[alphabet.size()];
      b[i] = new std::vector<ulint>[alphabet.size()];      
      Sa[i] = new TS[alphabet.size()];
      Sb[i] = new TS[alphabet.size()];
      Ia[i] = new bool[alphabet.size()];
      Ib[i] = new bool[alphabet.size()];      
      for (ulint j=0; j<alphabet.size()-1; j++) {
         Ia[i][j] = false;
         right = qidx.right_extension(alphabet[j],qsample);
         if (!right.is_invalid()) {
            left = qidx.left_extension(alphabet[i],right);         
            if (!left.is_invalid()) {
               Sa[i][j] = left;
               Ia[i][j] = true;
            }
         }   
         Ib[i][j] = false;
         right = idx.right_extension(alphabet[j],sample);
         if (!right.is_invalid()) {
            left = idx.left_extension(alphabet[i],right);
            if (!left.is_invalid()) {
               Sb[i][j] = left;
               Ib[i][j] = true;
            }
         }
      }
   }
   for (ulint i=0; i<alphabet.size()-1; i++) 
      for (ulint j=0; j<alphabet.size()-1; j++)
         if (Ia[i][j]) {
            for (ulint ii=0; ii<alphabet.size()-1; ii++)
               for (ulint jj=0; jj<alphabet.size()-1; jj++)
                  if (i!=ii and j!=jj and Ib[ii][jj]) {
                     if (a[i][j].size()==0)  // locate charged on the output
                        a[i][j] = qidx.locate_sample(Sa[i][j]);
                     if (b[ii][jj].size()==0)  // locate charged on the output
                        if (saidx==NULL) // no suffix array, using slower locate
                           b[ii][jj] = idx.locate_sample(Sb[ii][jj]);
                     if (d_bwt.size()==0) // no distance constraint, outputing all              
                        for (ulint iii=0; iii<a[i][j].size(); iii++)  
                           for (ulint jjj=0; jjj<b[ii][jj].size(); jjj++) {
                              MEMcount += 1;
                              if (output.is_open())
                                 output_node_MEM_gaf(b[ii][jj][jjj]+1, a[i][j][iii]+1, d, output, sp, efginfo, qinfo);
                           }
                     else { // outputing recursively using the distance constraint
                        pair<ulint,ulint> interval, leftinterval, rightinterval;
                        std::stack<pair<ulint,ulint>> intervalS;
                        interval.first = Sb[ii][jj].range.first;
                        interval.second = Sb[ii][jj].range.second;
                        intervalS.push(interval);
                        ulint argmin;
                        while (!intervalS.empty()) {
                           interval = intervalS.top();
                           intervalS.pop();
                           argmin = rmq(interval.first,interval.second);
                           if (d_bwt[argmin]<=d+1)
                              for (ulint iii=0; iii<a[i][j].size(); iii++) {
                                 MEMcount += 1;
				 if (output.is_open())
                                    output_edge_MEM_gaf(saidx[argmin]+1, a[i][j][iii]+1, d, output, sp, efginfo, qinfo);
                              }

                           leftinterval.first = interval.first;
                           leftinterval.second = argmin-1;
                           if (leftinterval.second>=leftinterval.first) 
                              intervalS.push(leftinterval);
                           rightinterval.first = argmin+1;
                           rightinterval.second = interval.second; 
                           if (rightinterval.second>=rightinterval.first) 
                              intervalS.push(rightinterval);
                        }
                     }
                  }   
         }

   for (ulint i=0; i<alphabet.size()-1; i++) {
      delete[] a[i];
      delete[] b[i];
      delete[] Sa[i];
      delete[] Sb[i];
      delete[] Ia[i];
      delete[] Ib[i];      
   }
   delete[] a;
   delete[] b;    
   delete[] Sa;
   delete[] Sb;    
   delete[] Ia;
   delete[] Ib;    

   return MEMcount;
}

template<class T, class TS>
ulint reportAMEMs(T qidx, T idx, TS qsample, TS sample, ulint d, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo)
{
   ulint MEMcount = 0;
   std::vector<ulint>** a= new  std::vector<ulint>*[alphabet.size()];
   bool** b= new  bool*[alphabet.size()];
   TS** Sa= new  TS*[alphabet.size()];
   bool** Ia= new  bool*[alphabet.size()];
   

   TS left,right;
   for (ulint i=0; i<alphabet.size()-1; i++) {
      a[i] = new std::vector<ulint>[alphabet.size()];
      b[i] = new bool[alphabet.size()];
      Sa[i] = new TS[alphabet.size()];
      Ia[i] = new bool[alphabet.size()];
      for (ulint j=0; j<alphabet.size()-1; j++) {
         Ia[i][j] = false;
         right = qidx.right_extension(alphabet[j],qsample);
         if (!right.is_invalid()) {
            left = qidx.left_extension(alphabet[i],right);
            if (!left.is_invalid()) {
               Sa[i][j] = left;
               Ia[i][j] = true;
            }
         }   
         b[i][j] = true; // maximal in text?         
         right = idx.right_extension(alphabet[j],sample);
         if (!right.is_invalid())
            b[i][j] = false; // not right-maximal  
         left = idx.left_extension(alphabet[i],sample);
         if (!left.is_invalid())
            b[i][j] = false; // not left-maximal
      }
   }
   std::vector<ulint> locations;
   for (ulint i=0; i<alphabet.size()-1; i++) 
      for (ulint j=0; j<alphabet.size()-1; j++)
         if (Ia[i][j] and b[i][j]) { 
            a[i][j] = qidx.locate_sample(Sa[i][j]);
            if (locations.size()==0)
               locations = idx.locate_sample(sample);
            for (ulint iii=0; iii<a[i][j].size(); iii++) {
               // reporting one occurrence in the text
               MEMcount += 1;
	       if (output.is_open())
                  output_edge_MEM_gaf(locations[0], a[i][j][iii]+1, d, output, sp, efginfo, qinfo);
            }
         }
   for (ulint i=0; i<alphabet.size()-1; i++) {
      delete[] a[i];
      delete[] b[i];
      delete[] Sa[i];
      delete[] Ia[i];
      
   }
   delete[] a;
   delete[] b;    
   delete[] Sa;
   delete[] Ia;

   return MEMcount;
}

template<class T, class TS>
pair<ulint,ulint> explore_mems(T tidx, T qidx, T fidx, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo, bool f = false, sdsl::int_vector<> d_bwt = sdsl::int_vector<>(), sdsl::rmq_succinct_sada<> rmq = sdsl::rmq_succinct_sada<>(), ulint* sa = NULL)
{
    TS sample(tidx.get_initial_sample(true));
    TS qsample(qidx.get_initial_sample(true));
    pair <TS,TS> node;
    std::stack<pair <TS,TS>> S; // interval pairs
    std::stack<ulint> dS; // string depths
    node.first = sample;
    node.second = qsample;
    // node is now suffix tree root
    S.push(node);
    ulint d = 0;
    dS.push(d); 
    bool MEM;
    ulint maxMEM = 0;
    ulint MEMcount = 0;
    ulint nodeCount = 0;
    
    TS fsample;
    if (f) // using fidx as filter
       fsample = fidx.get_initial_sample(true);
         
    std::stack<TS> fS; // filter text index range
    if (f)
       fS.push(fsample);
    
    while (!S.empty()) {
       nodeCount++;
       node = S.top();
       S.pop();
       d = dS.top();
       dS.pop();
       sample = node.first;
       qsample = node.second;
       
       if (f) {
          fsample = fS.top();
          fS.pop();
       }
       
       if (sample.is_invalid() or 
          qsample.is_invalid() or
          (!asymmetric and f and fsample.is_invalid())) { 
          continue; // not a valid range
       }
       if ((!is_right_maximal<T,TS>(tidx,sample) and !is_right_maximal<T,TS>(qidx,qsample)) and 
           tidx.bwt_at(sample.rangeR.first,true)==qidx.bwt_at(qsample.rangeR.first,true) ) {
          continue; // implicit node reached
       }
       MEM = 1;
       // Taking Weiner links from current node to visit nodes at string depth++
       // TODO: Push the largest interval first to limit the size of the stack
       // TODO: use as alphabet the distinct symbols appearing in the range

       // not branching with three last alphabet characters
       for (ulint j=0;j<alphabet.size()-3; j++) { 
          node.first = tidx.left_extension(alphabet[j],sample); 
          node.second = qidx.left_extension(alphabet[j],qsample);   
          S.push(node);       
          dS.push(d+1);
          if (f) {
             fS.push(fidx.left_extension(alphabet[j],fsample));
          }
          
          // range not splitting, no MEM reported
          if (sample.size()+qsample.size()==node.first.size()+node.second.size()) 
             MEM = 0;
       }
       if (MEM and d>=kappa) {
#ifdef EFG_MEMS_DEBUG
       cerr << "calling reportMEMs with depth equal to " << d;
       cerr << " on string ";
       for (ulint n = 0, i = qidx.FL(qsample.range.first); n < d; n++, i = qidx.FL(i)) {
          cerr << qidx.bwt_at(i);
       }
       cerr << std::endl;
#endif
          ulint mems;
	  if (!asymmetric) {
             mems = reportMEMs<T,TS>(tidx,qidx,sample,qsample,d,output,sp,efginfo,qinfo,d_bwt,rmq,sa);
          } else if (!f or fsample.is_invalid()) { // MEM string not found in filter index, so not a dublicate
             mems = reportAMEMs<T,TS>(tidx,qidx,sample,qsample,d,output,sp,efginfo,qinfo);
          }

          MEMcount += mems;
          if (mems > 0 && d > maxMEM)
                  maxMEM = d;
       }
    }
    cerr << "Number of recursion tree nodes visited: " << nodeCount << endl;
    return pair<ulint,ulint>(maxMEM,MEMcount);
}


template<class T, class TS>
pair<ulint,ulint> report_full_node_mems(string nodes, T qidx, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo)
{
   ulint i = nodes.size()-3;
   ulint d;
   bool report;
   ulint maxMEM = 0, MEMcount = 0;
   TS qsample;
   std::vector<ulint> locations;
   while (i>1) {
      qsample = qidx.get_initial_sample(true);
      d = 0;
      while (!qidx.left_extension(nodes[i],qsample).is_invalid() && nodes[i-1]!=alphabet[alphabet.size()-1]) {
         qsample = qidx.left_extension(nodes[i],qsample);
         report = true;
         i--;  
         d++;
      }
      if (d>1 and nodes[i-1]==alphabet[alphabet.size()-1]) { // TODO fix special case d=1
         // full node match
         locations = qidx.locate_sample(qsample);
         MEMcount += locations.size();
         for (ulint j=0; j < locations.size(); j++)
            //output << i+1 << "," << locations[j] << "," << d << endl; 
	    if (output.is_open())
               output_node_MEM_gaf(i+1, locations[j], d, output, sp, efginfo, qinfo);
         if (locations.size()>0 and d>maxMEM)
            maxMEM = d;  
      }
      // scanning to next node
      while (nodes[i]!=alphabet[alphabet.size()-1])
         i--;
      if (i>1)    
         i = i-2; // bypassing right-char
   }
   return pair<ulint,ulint>(maxMEM,MEMcount);
}


template<class T, class TS>
pair<ulint,ulint> report_suffix_mems(string edges, sdsl::int_vector<> d_edge, T qidx, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo)
{
   ulint i = edges.size()-3;
   ulint d;
   bool report;
   ulint maxMEM = 0, MEMcount = 0;
   TS qsample,q2sample;
   std::vector<ulint> locations;
   while (i>1) {
      qsample = qidx.get_initial_sample(true);
      d = 0;
      report = false;
      while (edges[i-1]!=alphabet[alphabet.size()-1]) {
         qsample = qidx.left_extension(edges[i],qsample);
         if (qsample.is_invalid())
            break;
         if (d_edge[i-1]==1) // found a node boundary
            report = true;
         i--; 
         d++;
         if (report) {
            if (d>maxMEM)
              maxMEM = d;             
            q2sample = qidx.left_extension(edges[i],qsample);
            if (q2sample.is_invalid()) { // no valid extensions, reporting the whole interval and stopping
               locations = qidx.locate_sample(qsample);
               MEMcount += locations.size();
               for (ulint jj=0; jj < locations.size(); jj++)
                     //output << i+1 << "," << locations[jj] << "," << d << endl;            
		     if (output.is_open())
                        output_edge_MEM_gaf(i+1, locations[jj], d, output, sp, efginfo, qinfo);
               break;
            }  
            // Valid extension exist, reporting the rest
            for (ulint j=0; j<alphabet.size()-1; j++)
               if (alphabet[j]!=edges[i]) {
                  q2sample = qidx.left_extension(alphabet[j],qsample);
                  if (q2sample.is_invalid())
                     continue;                 
                  locations = qidx.locate_sample(q2sample);
                  MEMcount += locations.size();
                  for (ulint jj=0; jj < locations.size(); jj++)
                     //output << i+1 << "," << locations[jj]+1 << "," << d << endl;
		     if (output.is_open())
                        output_edge_MEM_gaf(i+1, locations[jj]+1, d, output, sp, efginfo, qinfo);
               }              
         } 
      }
      // scanning to next edge
      while (edges[i]!=alphabet[alphabet.size()-1])
         i--;
      if (i>1)   
         i = i-2; // bypassing right-char
   }
   return pair<ulint,ulint>(maxMEM,MEMcount);
}

template<class T, class TS>
pair<ulint,ulint> report_prefix_mems(string edges, sdsl::int_vector<> d_edge, T qidx, ofstream& output, starting_positions &sp, gfa_info &efginfo, query_info &qinfo)
{
   ulint i = 2;
   ulint d;
   bool report;
   ulint maxMEM = 0, MEMcount = 0;
   TS qsample,q2sample;
   std::vector<ulint> locations;
   while (i<edges.size()) {
      qsample = qidx.get_initial_sample(true);
      d = 0;
      report = false;
      while (edges[i+1]!=alphabet[alphabet.size()-1]) {
         qsample = qidx.right_extension(edges[i],qsample);
         if (qsample.is_invalid())
            break;
         if (d_edge[i]==1) // found a node boundary
            report = true;
         i++; 
         d++;
         if (report) {
            if (d>maxMEM)
               maxMEM = d;
            q2sample = qidx.right_extension(edges[i],qsample);
            if (q2sample.is_invalid()) { // no valid extensions, reporting the whole interval and stopping
               locations = qidx.locate_sample(qsample);
               MEMcount += locations.size();
               for (ulint jj=0; jj < locations.size(); jj++) {
                     //output << i-d << "," << locations[jj] << "," << d << endl;            
		     if (output.is_open())
                        output_edge_MEM_gaf(i-d, locations[jj], d, output, sp, efginfo, qinfo);
               }
               break;
            }  
            // Valid extension exist, reporting the rest
            for (ulint j=0; j<alphabet.size()-1; j++)
               if (alphabet[j]!=edges[i]) {
                  q2sample = qidx.right_extension(alphabet[j],qsample);
                  if (q2sample.is_invalid())
                     continue;                 
                  locations = qidx.locate_sample(q2sample);
                  MEMcount += locations.size();
                  for (ulint jj=0; jj < locations.size(); jj++) {
                     //output << i-d << "," << locations[jj] << "," << d << endl;
		     if (output.is_open())
                        output_edge_MEM_gaf(i-d, locations[jj], d, output, sp, efginfo, qinfo);
                  }
               }              
         } 
      }
      // scanning to next edge
      while (edges[i]!=alphabet[alphabet.size()-1])
         i++;
      i = i+2; // bypassing left-char
   }
   return pair<ulint,ulint>(maxMEM,MEMcount);
}

template<class T>
void initialize_query_info(T &qidx, query_info &qinfo)
{
	// find length of string
	range_t fullrange = qidx.full_range();
	ulint length = fullrange.second - fullrange.first; // terminator char included

#ifdef EFG_MEMS_DEBUG
	// iterate over string
	cerr << "length of query index is " << length << " and query is \"" << endl;
	for (ulint i = qidx.FL(qidx.get_terminator_position()); i != qidx.get_terminator_position(); i = qidx.FL(i)) {
		cerr << qidx.bwt_at(i);
	}
	cerr << "\"" << endl;
#endif

	qinfo.concat.resize(length - 1);
	// TODO: assert
	qinfo.concat[0] = 0;
	ulint qlength = 0;
	for (ulint j = 0, i = qidx.FL(qidx.get_terminator_position()); qidx.FL(qidx.FL(i)) != qidx.get_terminator_position(); i = qidx.FL(i), j += 1) {
		char c = qidx.bwt_at(i);
		if (c == alphabet[alphabet.size()-2]) {
			qinfo.concat[j+1] = 1;
			if (qlength != 0) {
				qinfo.ordered_query_lengths.push_back(qlength);
			}
			qlength = 0;
		} else {
			qinfo.concat[j+1] = 0;
			qlength += 1;
		}
	}
	qinfo.ordered_query_lengths.push_back(qlength);

#ifdef EFG_MEMS_DEBUG
	cerr << "query lengths are ";
	for (auto l : qinfo.ordered_query_lengths) {
		cerr << l << " ";
	}
	cerr << endl;
#endif

	qinfo.concat_rs = sdsl::rank_support_v5<>(&qinfo.concat);
	qinfo.concat_ss = sdsl::select_support_mcl<>(&qinfo.concat);

	// TODO: change
	for (ulint i = 0; i < qinfo.concat_rs(length - 1); i++) {
		qinfo.ordered_query_ids.push_back("query" + std::to_string(i+1));
	}
}

template<class T, class TS>
void prepare_efg(string const &filename_efg, string &nodes_without_gt, string &edges_without_gt, T &nidx, T &eidx, sdsl::int_vector<> &d_edge, sdsl::int_vector<> &d_edge_bwt, sdsl::rmq_succinct_sada<> &rmq_edge, ulint *&sa_edge, starting_positions &sp, gfa_info &efginfo)
{
    string nodes, edges;
    ifstream efg_in(filename_efg);
    read_gfa_generic(efg_in,nodes,edges,sp,efginfo);
    efg_in.close();

    nodes_without_gt = string(nodes);    
    nodes_without_gt.erase(std::remove(nodes_without_gt.begin(), nodes_without_gt.end(), '>'), nodes_without_gt.end());
    edges_without_gt = string(edges);    
    edges_without_gt.erase(std::remove(edges_without_gt.begin(), edges_without_gt.end(), '>'), edges_without_gt.end());

    sp.nodes_rs = rank_support_v5<>(&sp.nodes);
    sp.edges_rs = rank_support_v5<>(&sp.edges);
    sp.nodes_ss = select_support_mcl<>(&sp.nodes);
    sp.edges_ss = select_support_mcl<>(&sp.edges);

    /*cout << "node select: ";
    for (int i = 1; i <= efginfo.node_lengths.size(); i++) {
	    cout << sp.nodes_ss(i) << "\t";
    }
    cout << endl;*/

    //bool* Be_suffix = new bool[edges_without_gt.size()];
    bit_vector Be_suffix(edges_without_gt.size(),0);
    //bool* Be_suffix_bwt = new bool[edges_without_gt.size()+1];    
    //bool* Be_prefix = new bool[edges_without_gt.size()];    
    bit_vector Be_prefix(edges_without_gt.size(),0);    
    //bool* Be_prefix_bwt = new bool[edges_without_gt.size()+1];        
    //bool* Bp_suffix = new bool[paths_without_gt.size()];
       
    bool suffix = true;
    bool prefix = false;

    ulint j=0;
    for (ulint i=0; i<edges.size(); i++)
       if (edges[i]==alphabet[alphabet.size()-1]) { // edge boundary
          suffix = true; 
          prefix = false;
          Be_suffix[j]=suffix;
          Be_prefix[j++]=prefix;          
       } 
       else if (edges[i]!='>') {
          Be_suffix[j]=suffix;
          Be_prefix[j++]=prefix;          
       }
       else { // node boundary
          suffix = false;
          prefix = true;
       }
    d_edge.resize(edges_without_gt.size()+1);        
    j = 0;
    for (ulint i=edges_without_gt.size()-1; i!=0; i--) {
       if (Be_prefix[i])
          j = 0;
       else 
          j++;
       if (!Be_suffix[i])
          d_edge[i]=edges_without_gt.size();
       else
          d_edge[i]=j;
    }
    d_edge[edges_without_gt.size()]=edges_without_gt.size();
  
    Be_prefix.empty();
    Be_suffix.empty();
    
    // releasing memory
    nodes.clear();
    edges.clear();
  
    /* Outputing the raw content for validation
    ofstream n_out(path_prefix+".nodes");
    n_out << nodes_without_gt;
    n_out.close();
    ofstream e_out(path_prefix+".edges");
    e_out << edges_without_gt;
    e_out.close();
    ofstream p_out(path_prefix+".paths");
    p_out << paths_without_gt;
    p_out.close();
    */

#ifdef EFG_MEMS_DEBUG
    cerr << "Node index is " << endl;
    cerr << nodes_without_gt << endl;

    cerr << "Node starting positions in the index are " << endl;
    for (ulint i=0; i < sp.nodes.size(); i++)
	    cerr << sp.nodes[i];
    cerr << endl;

    cerr << "Edge index is " << endl;
    cerr << edges_without_gt << endl;

    cerr << "Edge starting positions in the index are " << endl;
    for (ulint i=0; i < sp.edges.size(); i++)
	    cerr << sp.edges[i];
    cerr << endl;
#endif


    // Building and saving indexes if they don't exist
    if (use_brindex) {
       ifstream nidx_in(filename_efg+".nodes.bri");
       if (nidx_in.is_open()) {
          nidx.load(nidx_in);
          nidx_in.close();
       }
       else {
          nidx = T(nodes_without_gt);
          nidx.save_to_file(filename_efg + ".nodes"); 
       }
       ifstream eidx_in(filename_efg+".edges.bri");
       if (eidx_in.is_open()) {
          eidx.load(eidx_in);
          eidx_in.close();
       }
       else {
          eidx = T(edges_without_gt);
          eidx.save_to_file(filename_efg+".edges");
       } 
    }
    else { // just building bdbwt from scratch
       nidx = T(nodes_without_gt);
       eidx = T(edges_without_gt);
    }
    
    /* Converting d_edge to BWT order */
    // Computing suffix array also
    d_edge_bwt.resize(edges_without_gt.size()+1);    
    sa_edge = new ulint[edges_without_gt.size()+1];
    // assuming endmarker has lex-order 0
    j = eidx.LF(0);
    for (ulint i=0; i<edges_without_gt.size(); i++) {
       d_edge_bwt[j] = d_edge[edges_without_gt.size()-i-1];
       sa_edge[j]=edges_without_gt.size()-i-1;
       j = eidx.LF(j);    
    } 
    //d_edge.empty();
    
    rmq_edge = sdsl::rmq_succinct_sada<>(&d_edge_bwt);
        
    // releasing memory
    //nodes_without_gt.clear();
    //edges_without_gt.clear();
}

template<class T, class TS>
void prepare_query(string filename_qidx, T &qidx, std::stack<T> &qindexes, query_info &qinfo)
{
    if (use_brindex)
       qidx.load_from_file(filename_qidx);      
    else {
       // With bdbwt, creating indexes on the fly
       ifstream qidx_in(filename_qidx);
       string queries;
       getline(qidx_in,queries); // by default assuming input is one line
       //cout << queries << endl;
       if (queries[0]=='>') {// fasta
          fasta = true;
          while (getline(qidx_in,queries)) {
             //cout << queries << endl;
             if (queries.size()>1 and queries[0]!='>')
                qindexes.push(T("#"+queries+"#"));
          }
       }   
       else
          qidx = T(queries);
       qidx_in.close(); 
    }
    initialize_query_info(qidx, qinfo);
}

template<class T, class TS>
void find_mems(string filename_qidx, string filename_efg)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

    cerr << "Creating and loading indexes" << endl;
    auto tstart = high_resolution_clock::now();


    if (alphabet.size()==0) {
       alphabet = "ACGTN#$";
       // last char is a separator in graph concatenation, omitted in all algorithms
       // second last char is separator in queries, not used in MEM exploration, but used in cross-product
       // third last char marks ambiguous left- and right-extension, not used in MEM exploration, but used in cross-product
    }

    starting_positions sp;
    gfa_info efginfo;
    string nodes_without_gt, edges_without_gt;    
    T nidx, eidx;
    sdsl::int_vector<> d_edge; // distance from pos in start node to the start of end node
    sdsl::int_vector<> d_edge_bwt;
    ulint *sa_edge;
    sdsl::rmq_succinct_sada<> rmq_edge;
    prepare_efg<T,TS>(filename_efg, nodes_without_gt, edges_without_gt, nidx, eidx, d_edge, d_edge_bwt, rmq_edge, sa_edge, sp, efginfo);

    query_info qinfo;
    T qidx;
    std::stack<T> qindexes; // queries in fasta
    prepare_query<T,TS>(filename_qidx, qidx, qindexes, qinfo);

    ofstream output;
    if (output_file.size()!=0) {     
       output.open(output_file);
       if (!output) {
          cerr << "Could not open output file " << output_file << endl;
       }
    }

    auto tmem = high_resolution_clock::now();       
    cerr << "Exploring MEMs " << endl;

    ulint maxMEM;
    
    if (!fasta) // queries as concatenation
       qindexes.push(qidx);
    
    while (!qindexes.empty()) {
       qidx = qindexes.top();
       qindexes.pop();

       pair<ulint,ulint> result = explore_mems<T,TS>(nidx,qidx,nidx,output,sp,efginfo,qinfo,false);
       maxMEM = result.first;
       ulint nodeMEMs = result.second;
       cerr << "Found " << nodeMEMs << " node " << kappa << "-MEM of length at most " << maxMEM << endl;

       ulint edgeMEMs;
       if (!asymmetric) {
          result = explore_mems<T,TS>(eidx,qidx,nidx,output,sp,efginfo,qinfo,false,d_edge_bwt,rmq_edge,sa_edge);
       } else { // using node index as filter
          result = explore_mems<T,TS>(eidx,qidx,nidx,output,sp,efginfo,qinfo,true,d_edge_bwt,rmq_edge,sa_edge);
       }
       maxMEM = result.first;
       edgeMEMs = result.second;
       cerr << "Found " << edgeMEMs << " edge " << kappa << "-MEM of length at most " << maxMEM << endl;


       // full nodes MEMs are include in edge prefix and suffix MEMs
       //output << ">full node MEMs" << endl;
       //result = report_full_node_mems<T,TS>(nodes_without_gt,qidx,output,sp,efginfo,qinfo);
       //maxMEM = result.first;
       //ulint fullnodeMEMs = result.second;
       //cout << "Maximum full node MEM is of length " << maxMEM << endl;
       //cerr << "Found " << edgeMEMs << " full node MEMs of length at most " << maxMEM << endl;
    

       result = report_suffix_mems<T,TS>(edges_without_gt,d_edge,qidx,output,sp,efginfo,qinfo);
       maxMEM = result.first;
       ulint edgesufMEMs = result.second;
       cerr << "Found " << edgesufMEMs << " edge suffix MEMs of length at most " << maxMEM << endl;
    
       result = report_prefix_mems<T,TS>(edges_without_gt,d_edge,qidx,output,sp,efginfo,qinfo);
       maxMEM = result.first;
       ulint edgepreMEMs = result.second;
       cerr << "Found " << edgepreMEMs << " edge prefix MEMs of length at most " << maxMEM << endl;
    }
    output.close();
    d_edge_bwt.empty();
    delete[] sa_edge;

    auto tend = high_resolution_clock::now();
    ulint tot_time = duration_cast<microseconds>(tend-tstart).count(); 
    ulint mem_time = duration_cast<microseconds>(tend-tmem).count(); 
    cerr << "MEM exploration time     : " << mem_time << " microseconds" << endl;
    cerr << "Total time     : " << tot_time << " microseconds" << endl;
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    int ptr = 1;

    while (ptr < argc - 3) 
       parse_args(argv, argc, ptr);

    string filename_qidx(argv[ptr]);
    string filename_efg(argv[ptr+1]);
    
    if (use_brindex)      
       find_mems<br_index<>,br_sample>(filename_qidx, filename_efg);
    else
       find_mems<bdbwt_index,bd_sample>(filename_qidx, filename_efg);
}
