
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

int main(int argc, char** argv) {
    /******************
     * parameters are:
     *  - dataset: requires that <dataset>.vco already exists
     *  - nrows: number of rows in the dataset
     *  - ncols: number of columns in the dataset
     *  - solution file: any .sol file prepresenting a column reordering
     * output:
     *   <dataset>.vc contains the reordered file
     ******************/
    if (argc!=4+1){
        std::cout << "Usage is: " << argv[0] << " <dataset> <nrows> <ncols> <solution_file>\n";
        return -1;
    }
    std::string dataset = argv[1];
    unsigned nrows = atoi(argv[2]);
    unsigned ncols = atoi(argv[3]);
    std::string solution_path = argv[4];

    /********************************
     * read permutation (sol file)
     *******************************/

    std::ifstream solfile ;
    solfile.open(solution_path);
    assert(solfile.is_open());
    std::vector<uint32_t> perm;
    perm.reserve(ncols);
    std::string pattern = "TOUR_SECTION", line;
    std::getline(solfile, line) ;
    do std::getline(solfile, line); while(line != pattern);
    for(auto col = 0; col<ncols ; ++col){
        std::getline(solfile, line);
        uint32_t p = std::stoul(line);
        p -= 1; //0-based instead of 1-based
        perm.push_back(p);
    }
    std::getline(solfile,line);
    assert(line == "-1");
    assert(perm.size() == ncols);
    solfile.close();

#if DEBUG
    for(auto p : perm)
        std::cout << p << " ";
    std::cout << std::endl;
#endif

    /********************************
     * write permuted vc file
     *******************************/
    std::ifstream infile;
    std::ofstream outfile;
    infile.open(dataset + ".vco");
    outfile.open(dataset + ".vc",std::ofstream::binary);
    // outfile.open(solution_path + ".vc",std::ofstream::binary);

    assert(infile.is_open());
    assert(outfile.is_open());

    std::unordered_map<uint32_t,uint32_t> row_map;
#ifdef DEBUG
    for(auto row=0; row<1; ++row){
#else
    for(auto row=0; row<nrows; ++row){
#endif
        row_map.clear();
        row_map.reserve(ncols);

        uint32_t tmp;
        infile.read(reinterpret_cast<char *>(&tmp), sizeof(tmp)); //fst int

        //read
#ifdef DEBUG
        std::cout << "reading ..." << std::endl;
#endif
        while (tmp){
            tmp -= 1 ;
            uint32_t pos = tmp / ncols ;
            uint32_t col = tmp % ncols ;
            assert(col<ncols);
            assert(row_map.find(col) == row_map.end());
            row_map[col] = pos ;

            infile.read(reinterpret_cast<char *>(&tmp), sizeof(tmp));
#if DEBUG
            std::cout << tmp << ", pos: " << pos << ", col: " << col << std::endl;
#endif
        }
#if DEBUG
        std::cout << "EOS" << std::endl;
#endif
        
        //write
#ifdef DEBUG
        std::cout << "writing ..." << std::endl;
#endif
        for(auto col : perm) {
            auto search = row_map.find(col) ;
            if (search != row_map.end()){
                unsigned pos = search->second;
                tmp = pos*ncols + col + 1 ;
                outfile.write(reinterpret_cast<char *>(&tmp), sizeof(tmp));
#if DEBUG
                std::cout << tmp << ", pos: " << pos << ", col: " << col << std::endl;
#endif
            }
        }
#if DEBUG
        std::cout << "EOS" << std::endl;
#endif
        tmp  = 0x0;
        outfile.write(reinterpret_cast<char *>(&tmp), sizeof(tmp));
        outfile.flush();
    }
    infile.close();
    outfile.close();
}//end main
