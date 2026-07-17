#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include <cassert>
#include <thread>

//boost
#include <boost/program_options.hpp>
namespace po = boost::program_options;

/*
Computing the column similarity matrix via local pruning.
*/

#define MAXUI32 0x7fffffff
#define MAXU32 ~0x0

constexpr size_t MAX_PAIRS_PER_MEGABYTE = (1 << 20)/(8+4);

//a triple, representing a column pair and the similarity score.
struct col_col_wth {
    uint32_t c1=MAXU32,c2=MAXU32;
    int32_t w=MAXUI32;

    bool has_value() {
        return c1!=MAXU32 or c2!=MAXU32 or w!=MAXUI32;
    }
    friend bool operator<(const col_col_wth& l, const col_col_wth& r)
    {
        return l.w<r.w;
    }
};
typedef struct col_col_wth ccw_t;
typedef std::vector<ccw_t> max_heap_t;

//reading the data set
void inline read_csv(
    std::string &csv_path,
    size_t nrows,
    size_t ncols,
    std::vector<std::vector<double>> &vs,
    std::vector<std::vector<uint32_t>> &rs,
    size_t mega_thrs
    )
{
    std::ifstream ifs;
    ifs.open (csv_path, std::ifstream::in);
    assert(ifs.is_open());

    //struct init
    vs.clear(); vs.reserve(ncols);
    rs.clear(); rs.reserve(ncols);
    for(size_t c=0; c<ncols; ++c) {
        vs.emplace_back();
        rs.emplace_back();
    }

    //params
    size_t cur_pairs=0;
    //if mega_thrs is null, there's no limit.
    const size_t npairs = mega_thrs ? mega_thrs * MAX_PAIRS_PER_MEGABYTE : std::numeric_limits<size_t>::max();
    std::string token;

    //business-logic code
    for(size_t r=0; r<nrows and cur_pairs<npairs; ++r){
        for(size_t c=0; c<ncols and cur_pairs<npairs; ++c){
            const char delim = (c+1 == ncols) ? '\n' : ',';
            std::getline(ifs, token, delim);
            double v = stod(token);
            vs[c].push_back(v);
            rs[c].push_back(r);
            ++ cur_pairs;
        }
    }

    //closing
    ifs.close();
    for(size_t c=0; c<ncols; ++c) {
        vs[c].shrink_to_fit();
        rs[c].reserve(rs[c].size());
        assert(vs[c].size() == rs[c].size());
    }
    return;
}

void compute_local_simmat(
    std::string &filename,
    const size_t ncols,
    std::vector<std::vector<double>> &vs,
    std::vector<std::vector<uint32_t>> &rs,
    const unsigned k_par,
    const unsigned tid,
    const unsigned ts,
    std::vector<max_heap_t> &max_heaps
) {
    max_heaps.clear();

    uint32_t heap_size = std::min<uint32_t>(k_par, ncols - 1);
    max_heaps.reserve(ncols);
    for (auto i = 0; i < ncols; ++i)
        max_heaps.emplace_back(heap_size);

    for (uint32_t c1 = tid; c1 < ncols - 1; c1+=ts)
        //std::for_each(std::execution::par, uint32_t(0), uint32_t(ncols), [&](uint32_t c1)
    {
        for (uint32_t c2 = c1 + 1; c2 < ncols; ++c2) {
            std::vector<std::pair<double, double>> vpairs;
            //merge-sort-like intersection
            unsigned i1 = 0, i2 = 0;
            while (i1 < rs[c1].size() and i2 < rs[c2].size()) {
                const uint32_t r1 = rs[c1][i1], r2 = rs[c2][i2];

                if (r1 < r2) {
                    ++i1;
                } else if (r1 > r2) {
                    ++i2;
                } else {
                    const double v1 = vs[c1][i1], v2 = vs[c2][i2];
                    vpairs.emplace_back(v1, v2);
                    ++i1;
                    ++i2;
                }
            }
            std::sort(vpairs.begin(), vpairs.end());
            auto iter = vpairs.begin() + 1; //skip first pair
            while (
                    iter < vpairs.end() and
                    iter->first == 0 and
                    iter->second == 0
                    )
                ++iter;
            int32_t distance = 0;//cs;
            for (; iter < vpairs.end(); ++iter) {
                if (
                        iter->first == (iter - 1)->first and
                        iter->second == (iter - 1)->second and
                        iter->first != 0 and iter->second != 0
                        ) {
                    --distance;
                }
            }
            for (auto c: {c1, c2}) {
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
                    std::push_heap(max_heaps[c].begin(), max_heaps[c].end());
                }
            }
        }
    } //end parallel_for
    //});
};

void compute_global_simmat(
        const size_t ncols
        , const unsigned k_par
        , const unsigned ts
        , std::vector<std::vector<max_heap_t>> &max_heaps
)
{
    //single thread hereafter
    //for(size_t c=tid; c<cs; c+=ts){
    for(size_t c=0; c<ncols; ++c){
        //aggregating everything to max_heaps[0][c]
        max_heaps[0][c].reserve(k_par * ts);
        for(size_t t=1; t<ts; ++t) {
            max_heaps[0][c].insert(
                    max_heaps[0][c].end(),
                    max_heaps[t][c].begin(),
                    max_heaps[t][c].end()
                    );
        }
        //find top elements
        std::sort(max_heaps[0][c].begin(), max_heaps[0][c].end());
        max_heaps[0][c].resize(k_par);
    }
}


void write_simmat(
        std::vector<max_heap_t> &max_heaps,
        std::string &filepath,
        unsigned k_par,
        unsigned cs
        ){
    //extract ccw from all heaps
    std::vector<ccw_t> edges;
    for(auto max_heap : max_heaps){
        for(auto e : max_heap){
            if (e.has_value())
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
    outfile.open (filepath + ".pruned_local_" + std::to_string(k_par) + ".tsp");
    assert(outfile.is_open());

    outfile << "NAME: " << filepath << std::endl;
    outfile << "COMMENT: Column reordering problem, pruned local version" << std::endl;
    outfile << "TYPE : TSP" << std::endl;
    outfile << "DIMENSION : " << std::to_string(cs) << std::endl;
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
    for(auto c1=0;c1<cs-1;++c1){
        for(auto c2=c1+1;c2<cs;++c2){
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
    outfile.open (filepath + ".pruned_local_" + std::to_string(k_par) + ".par");
    assert(outfile.is_open());

    outfile << "PROBLEM_FILE = " << filepath << ".pruned_local_" << std::to_string(k_par) << ".tsp" << std::endl;
    outfile << "MOVE_TYPE = 5" << std::endl;
    outfile << "PATCHING_C = 3" << std::endl;
    outfile << "PATCHING_A = 2" << std::endl;
    outfile << "RUNS = 10" << std::endl;
    outfile << "TOUR_FILE = " << filepath << ".pruned_local_" << std::to_string(k_par) << ".tsp.solution" << std::endl;
    outfile.flush();
    outfile.close();
}

int main(int argc, char** argv) {
    if (argc!=6+1){
        std::cout << "Usage is: " << argv[0] << " <filepath> <nrows> <ncols> <k par> <ts> <mega thrs>" << std::endl;
        return 1;
    }

    //args
    std::string filepath(argv[1]);
    const unsigned
        nrows = std::stoi(argv[2]),
        ncols = std::stoi(argv[3]),
        k_par = std::stoi(argv[4]),
        ts = std::stoi(argv[5]) ,
        mega_thrs = std::stoi(argv[6]); //mega_thrs==0 means no memory limit

    //data
    std::vector<std::vector<double>> vs;
    std::vector<std::vector<uint32_t>> rs;
    std::vector<std::vector<max_heap_t>> hss(ts);

    //synchro
    std::vector<std::thread> threads(ts);

    read_csv(filepath,nrows,ncols,vs,rs, mega_thrs);
    for(unsigned tid=0; tid<ts; ++tid) {
        threads[tid] = std::thread(
                compute_local_simmat,
                std::ref(filepath), ncols, std::ref(vs), std::ref(rs), k_par, tid, ts, std::ref(hss[tid])
                );
    }
    for(unsigned tid=0; tid<ts; ++tid) {
        threads[tid].join();
    }
    compute_global_simmat(ncols, k_par, ts, hss);
    write_simmat(hss[0], filepath, k_par, ncols);

    return 0;
}
