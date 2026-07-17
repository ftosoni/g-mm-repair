#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <cstring>

//fetch data from disc
void fetch_block(size_t tid, size_t rows, size_t cols, size_t blocks, std::vector<double> &ds, std::string &matrix_path){
    std::string infilepath;
    infilepath += matrix_path;
    infilepath += ".";
    infilepath += std::to_string(blocks);
    infilepath += ".";
    infilepath += std::to_string(tid);

    const size_t row_block_size = (rows + blocks - 1) / blocks;
    const size_t row_bgn = tid * row_block_size;
    const size_t row_end = std::min<size_t>((tid + 1) * row_block_size, rows);
    assert(row_end > row_bgn);
    const size_t to_fetch = cols * (row_end-row_bgn);

    ds.resize(to_fetch);
    std::ifstream infile(infilepath, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "Could not open the file - '"
                  << infilepath << "'" << std::endl;
        exit(EXIT_FAILURE);
    }

    char *ptr = (char *)(&ds[0]);
    while (infile.get(*ptr)) {
        ++ ptr;
    }
    assert(ptr == (char *)(&ds[to_fetch]));
    infile.close();
    return;
};

void fetch_vector(size_t cols, std::vector<double> &cs, std::string &y_path){
    cs.resize(cols);
    std::ifstream infile(y_path, std::ios::in);
    const size_t to_fetch = cols;

    if (!infile.is_open()) {
        std::cerr << "Could not open the file - '"
                  << y_path << "'" << std::endl;
        exit(EXIT_FAILURE);
    }

    char *ptr = (char *)(&cs[0]);
    while (infile.get(*ptr)) {
        ++ ptr;
    }
    assert(ptr == (char *)(&cs[to_fetch]));
    infile.close();
    return;
};

//business-logic functions
inline void normalise_lambda(size_t tid, size_t cols, size_t blocks, std::vector<double> &cs, double lambda){
    const size_t col_block_size = (cols + blocks - 1) / blocks;
    const size_t col_bgn = tid * col_block_size;
    const size_t col_end = std::min<size_t>((tid + 1) * col_block_size, cols);

    for(size_t c=col_bgn; c<col_end; ++c ){
        assert(c < cs.size());
        cs[c] /= lambda;
    }
};

inline void right_mult(size_t tid, size_t rows, size_t cols, size_t blocks, std::vector<double> &ds, std::vector<double> &rs, std::vector<double> &cs){
    const size_t row_block_size = (rows + blocks - 1) / blocks;
    const size_t row_bgn = tid * row_block_size;
    const size_t row_end = std::min<size_t>((tid + 1) * row_block_size, rows);

    for(size_t row=row_bgn; row<row_end; ++row){
        rs[row] = 0;
    }
    for(size_t r=row_bgn; r<row_end; ++r){
        for(size_t c=0; c<cols; ++c) {
            size_t pos = (r-row_bgn)*cols + c;
            assert(pos < ds.size());
            rs[r] += ds[pos] * cs[c];
        }
    }
    return;
};

inline void left_mult(size_t tid, size_t rows, size_t cols, size_t blocks, std::vector<double> &ds, std::vector<double> &rs, std::vector<std::vector<double>> &css){
    const size_t row_block_size = (rows + blocks - 1) / blocks;
    const size_t row_bgn = tid * row_block_size;
    const size_t row_end = std::min<size_t>((tid + 1) * row_block_size, rows);

    for(size_t r=row_bgn; r<row_end; ++r){
        for(size_t c=0; c<cols; ++c) {
            size_t pos = (r-row_bgn)*cols + c;
            assert(pos < ds.size());
            css[tid][c] += ds[pos] * rs[r];
        }
    }
    return;
};

inline void local_reduce(size_t tid, size_t cols, size_t blocks, std::vector<double> &cs, std::vector<std::vector<double>> &css, std::vector<double> &local_lambdas) {
    const size_t col_block_size = (cols + blocks - 1) / blocks;
    const size_t col_bgn = tid * col_block_size;
    const size_t col_end = std::min<size_t>((tid + 1) * col_block_size, cols);

    double lambda = std::numeric_limits<double>::min();
    for (size_t c = col_bgn; c < col_end; ++c) {
        assert(c < cs.size());
        cs[c] = 0;
        for (size_t tid = 0; tid < blocks; ++tid) {
            assert(tid < css.size());
            assert(c < css[tid].size());
            cs[c] += css[tid][c];
            css[tid][c] = 0;
        }
        if (cs[c] > lambda) lambda = cs[c];
    }
    local_lambdas[tid] = lambda;
    return;
};

inline double global_reduce(std::vector<double> &local_lambdas){
    double lambda = std::numeric_limits<double>::min();
    for (double local_lambda : local_lambdas)
        if (local_lambda > lambda)
            lambda = local_lambda;
    return lambda;
};

//synchronisation functions
inline void join_all (std::vector<std::thread> &ts) {
    for(auto &t : ts) t.join();
};

int main(int argc, char** argv) {
    if (argc != 5 + 1) {
        std::cerr << "Usage: " << argv[0] << " <matrix base path> <y path> <rows> <cols> <par. degree>" << std::endl;
        exit(1);
    }

    //args
    std::string matrix_path(argv[1]), y_path(argv[2]);
    const size_t rows = strtol(argv[3], NULL, 10), cols = strtol(argv[4], NULL, 10), blocks = strtol(argv[5], NULL, 10);

    //business-logic structures
    std::vector<double> rs(rows, 0), cs;
    std::vector<std::vector<double>> dss(blocks);
    double lambda = 1;

    //synchronisation structures
    std::vector<std::thread> ts(blocks);
    std::vector<std::vector<double>> css; //map-reduce aux vects
    css.reserve(blocks);
    for(size_t tid=0; tid<blocks; ++tid){
        css.emplace_back(cols, 0x0);
    }
    std::vector<double> local_lambdas(blocks);

    fetch_vector(cols, cs, y_path);
    for(size_t tid=0; tid<blocks; ++tid)
        ts[tid] = std::thread(
            fetch_block
            ,tid
            ,rows
            ,cols
            ,blocks
            ,std::ref(dss[tid])
            ,std::ref(matrix_path)
            );
    join_all(ts);

    const size_t ITERATIONS = 500;
    for(size_t i=0; i<ITERATIONS; ++i) {
        for(size_t tid=0; tid<blocks; ++tid)
            ts[tid] = std::thread(
                normalise_lambda
                ,tid
                //,rows
                ,cols
                ,blocks
                ,std::ref(cs)
                ,lambda
                );
        join_all(ts);
        for(size_t tid=0; tid<blocks; ++tid)
            ts[tid] = std::thread(
                right_mult
                ,tid
                ,rows
                ,cols
                ,blocks
                ,std::ref(dss[tid])
                ,std::ref(rs)
                ,std::ref(cs)
                );
        join_all(ts);
        for(size_t tid=0; tid<blocks; ++tid)
            ts[tid] = std::thread(
                left_mult
                ,tid
                ,rows
                ,cols
                ,blocks
                ,std::ref(dss[tid])
                ,std::ref(rs)
                ,std::ref(css)
                );
        join_all(ts);
        for(size_t tid=0; tid<blocks; ++tid)
            ts[tid] = std::thread(
                local_reduce
                ,tid
                //,rows
                ,cols
                ,blocks
                ,std::ref(cs)
                ,std::ref(css)
                ,std::ref(local_lambdas)
                );
        join_all(ts);
        lambda = global_reduce(local_lambdas);
    }
    std::cout << "iterations: " << ITERATIONS << std::endl;
    std::cout << "lambda: " << lambda << std::endl;
}