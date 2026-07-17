#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include <cassert>

/*
Computing the column similarity matrix via local pruning.
*/

typedef std::vector<uint32_t> seq_t;

#define MAXUI32 0x7fffffff 
#define MAXU32 ~0x0 

//a triple, representing a column pair and the similarity score.
struct col_col_wth {
    uint32_t c1=MAXU32,c2=MAXU32;
    int32_t w=MAXUI32;

    friend bool operator<(const col_col_wth& l, const col_col_wth& r)
    {
        return l.w<r.w;
    }
};
typedef struct col_col_wth ccw_t;

typedef std::priority_queue<ccw_t> heap;

//reading the data set
void inline read_matrix(std::string filepath, size_t nrows, size_t ncols, seq_t &ss, seq_t &vs, seq_t &rs){
    //we use a compressed representation.
    ss.clear(); //pointers to the starting position of each column
    vs.clear(); //values
    rs.clear(); //rows indexes

    ss.reserve(ncols+1);
    std::string line;

    for(unsigned c=0;c<ncols;++c){
        ss.push_back(vs.size());

        std::ifstream istrm(filepath+"_cols/"+std::to_string(c)+".txt");
        assert(!istrm.fail());
        while (std::getline(istrm,line,'@')){
            unsigned val, row;
            val = (unsigned)stoul(line);//stoi(line);
            vs.push_back(val);
            std::getline(istrm,line);
            row = (unsigned)stoul(line);//stoi(line);
            rs.push_back(row);
        }
        istrm.close();
    }
    ss.push_back(vs.size());
}

//writing the similarity matrix to a file
void inline write_simmat(std::string filename, size_t nrows, size_t ncols, seq_t &ss, seq_t &vs, seq_t &rs, const unsigned k_par){
    
    std::vector<std::vector<ccw_t>> max_heaps;
    uint32_t heap_size = std::min<uint32_t>(k_par,ncols-1);
    max_heaps.reserve(ncols);
    for(auto i=0;i<ncols;++i)
        max_heaps.emplace_back(heap_size);

    for(uint32_t c1=0; c1<ncols-1; ++c1){
        for(uint32_t c2=c1+1; c2<ncols; ++c2){
            std::vector<std::pair<unsigned,unsigned>> vpairs;
            //merge-sort-like intersection
            unsigned i1=ss[c1],i2=ss[c2];
            while (i1<ss[c1+1] && i2<ss[c2+1]) {
                if(rs[i1]<rs[i2]){
                    ++ i1;
                }else if(rs[i1]>rs[i2]){
                    ++ i2;
                }else{
                    vpairs.emplace_back(vs[i1],vs[i2]);
                    ++ i1;
                    ++ i2;
                }
            }
            std::sort(vpairs.begin(),vpairs.end());
            auto iter = vpairs.begin()+1; //skip first pair
            while(
                iter < vpairs.end() &&
                iter->first == 0 &&
                iter->second == 0
            ) ++iter;
            int32_t distance = 0;//cs;
            for(; iter < vpairs.end(); ++iter){
                if (
                    iter->first==(iter-1)->first &&
                    iter->second==(iter-1)->second &&
                    iter->first != 0 && iter->second!= 0
                ){
                    -- distance;
                }
            }
            for (auto c : {c1,c2}) {
                auto top = max_heaps[c].front();
                if (top.w > distance) {
                    //pop operations
                    std::pop_heap(max_heaps[c].begin(), max_heaps[c].end());
                    max_heaps[c].pop_back();

                    ccw_t newv;
                    newv.c1 = c1;
                    newv.c2 = c2;
                    newv.w = distance;

                    //push operations
                    max_heaps[c].push_back(newv);
                    std::push_heap(max_heaps[c].begin(),max_heaps[c].end());
                }
            }
        }
    }

    //extract ccw from all heaps
    std::vector<ccw_t> edges;
    for(auto max_heap : max_heaps){
        for(auto e : max_heap){
            edges.push_back(e);
        }
    }

    auto row_wise_cmp = [](ccw_t e, ccw_t other){
        return  !(
          e.c1 < other.c1 || 
          (e.c1  == other.c1) && e.c2 < other.c2
        );
    };
    std::sort<>(edges.begin(),edges.end(),row_wise_cmp);

    std::ofstream outfile;
    
    //TSP
    outfile.open (filename+".pruned_local_"+std::to_string(k_par)+".tsp");
    assert(outfile.is_open());

    outfile << "NAME: " << filename << std::endl;
    outfile << "COMMENT: Column reordering problem, pruned local version" << std::endl;
    outfile << "TYPE : TSP" << std::endl;
    outfile << "DIMENSION : " << std::to_string(ncols) << std::endl;
    outfile << "EDGE_WEIGHT_TYPE: EXPLICIT" << std::endl;
    outfile << "EDGE_WEIGHT_FORMAT : UPPER_ROW" << std::endl;
    outfile << "EDGE_WEIGHT_SECTION" << std::endl;

    auto pop_next = [&edges](){
        if (edges.empty()) return ccw_t{};
        auto back = edges.back();
        edges.pop_back();
        return back;
    };
    ccw_t next = pop_next();
    for(auto c1=0;c1<ncols-1;++c1){
        for(auto c2=c1+1;c2<ncols;++c2){
            while (next.c1 < c1 || (next.c1 == c1 && next.c2<c2)){
                next = pop_next();
            }
            assert(!(next.c1 < c1 || (next.c1 == c1 && next.c2<c2)));
            if (next.c1 == c1 && next.c2 == c2) {
                outfile << std::to_string(next.w);
                next = pop_next();
            }else{
                outfile << "0";
            }
            outfile << " ";
        }
        outfile << "\n";
        outfile.flush();

    }
    outfile << "EOF" << std::endl;
    outfile.close();

    // PARAMETER FILE
    outfile.open (filename+".pruned_local_"+std::to_string(k_par)+".par");
    assert(outfile.is_open());

    outfile << "PROBLEM_FILE = " << filename << ".pruned_local_" << std::to_string(k_par) << ".tsp" << std::endl;
    outfile << "MOVE_TYPE = 5" << std::endl;
    outfile << "PATCHING_C = 3" << std::endl;
    outfile << "PATCHING_A = 2" << std::endl;
    outfile << "RUNS = 10" << std::endl;
    outfile << "TOUR_FILE = " << filename << ".pruned_local_" << std::to_string(k_par) << ".tsp.solution" << std::endl;
    outfile.flush();
    outfile.close();
}

int main(int argc, char** argv) {
    if (argc!=4+1){
        std::cout << "Usage is: " << argv[0] << " <filepath> <nrows> <ncols> <k par>\n";
        return -1;
    }

    std::string filepath(argv[1]);
    const unsigned nrows = atoi(argv[2]), ncols = atoi(argv[3]), k_par = atoi(argv[4]);

    //matrix ms;
    seq_t ss,vs,rs;

    read_matrix(filepath,nrows,ncols,ss,vs,rs);
    write_simmat(filepath,nrows,ncols,ss,vs,rs,k_par);


}
