#include <ostream>
#include <iostream>
#include <fstream>
#include <cassert>

/*
 * The maximum weighed matching (MWM) approach to column reordering.
 */

//boost for SN and matching
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
typedef boost::property< boost::edge_weight_t, float, boost::property< boost::edge_index_t, int > >
    EdgeProperty;
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeProperty > Graph;
#ifdef RING
typedef boost::graph_traits<Graph>::edge_descriptor EdgeDescriptor ;
#endif


int main(int argc, char *argv[]) {
    if(argc!=2+1){
        std::cout << "Usage is: " << argv[0] << " <dataset> <generator ext.>" << std::endl;
        return -1;
    }

    std::string dataset(argv[1]);
    std::string gen(argv[2]);


    std::string line, pattern;
    std::ifstream tspfile (dataset + "." + gen + ".tsp");
    assert(tspfile.is_open());
    pattern = "DIMENSION : ";
    while (! line.rfind(pattern,0)==0)
    {
        assert(! tspfile.eof());
        std::getline (tspfile,line);
    }
    std::string num_vertices_str(line.begin()+pattern.size(),line.end());
    const unsigned num_vertices = std::stoi(num_vertices_str);
    //std::cout << "#vertices: " << num_vertices << std::endl;
    pattern = "EDGE_WEIGHT_SECTION";
    while (! line.compare(pattern)==0)
    {
        assert(! tspfile.eof());
        std::getline (tspfile,line);
    }
    
    Graph SN(num_vertices * 2);
    for(unsigned c1=0; c1<=num_vertices-2; ++c1){
        for(unsigned c2=c1+1; c2<=num_vertices-1; ++c2){
            std::getline(tspfile,line,' ');
            //In the following we change the sign of the weight, so that we deal with positive values.
            //We shift values rightwards as to avoid overflow errors at a later stage.
            //Finally, we increment w by one unit to avoid zero-weight edges.
            unsigned w = (-std::stoi(line) >> 8) + 1 ; 
            /*
             * Every edge represent a possible predecessor-successor relationship
             * between two columns. Predecessor columns have id in [0, num_vertices),
             * whereas successors take values in [num_vertices, 2*num_vertices). 
             */
            boost::add_edge(c1,c2+num_vertices,EdgeProperty(w),SN); 
            //std::cout << "paro: " << c1 << " " << c2 << " " << w << std::endl; 
#ifdef RING
            boost::add_edge(c2,c1+num_vertices,EdgeProperty(w+1),SN); //unitary increment to force matching
#endif
        }
        std::getline(tspfile,line); //useless empty string
    }

    tspfile.close();

    std::vector< boost::graph_traits< Graph >::vertex_descriptor > mate1(
        num_vertices * 2);
    boost::maximum_weighted_matching(SN, &mate1[0]);
    /*std::cout << "Found a weighted matching:" << std::endl;
    std::cout << "Matching size is " << matching_size(SN, &mate1[0])
              << ", total weight is " << matching_weight_sum(SN, &mate1[0])
              << std::endl;
    std::cout << std::endl;*/

    std::ofstream pairfile;
#ifdef RING
    pairfile.open(dataset+ "." + gen + ".ring_pairs");
#else
    pairfile.open(dataset+ "." + gen + ".pairs");
#endif

    //writing pair file
    for(size_t fst=0; fst<num_vertices; ++fst){
        std::string fst_str = std::to_string(fst);
        bool is_matched = mate1[fst]!=-1;
        size_t snd = mate1[fst]-num_vertices;
#ifdef RING
        assert(is_matched);
        std::string snd_str = std::to_string(snd);
        std::pair<EdgeDescriptor, bool> e_descr = boost::edge(fst, mate1[fst], SN);
        assert(e_descr.second);
        int weight = boost::get(boost::edge_weight_t(), SN, e_descr.first);
        std::string w_str = std::to_string(weight);
        pairfile << fst_str << " " << snd_str << " " << w_str << std::endl;
#else
        std::string snd_str = is_matched ? std::to_string(snd) : "-1";
        pairfile << fst_str << " " << snd_str << std::endl;
#endif

    }

    pairfile.flush();
    pairfile.close();
   
}
