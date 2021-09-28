#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <sstream>

using namespace std;
std::random_device rd;
std::mt19937 generator(rd());

int weighted_depth_choice(std::vector<std::pair<int,float>> depm) {
    std::uniform_real_distribution<float> rng(0, 1);
    float rn = rng(generator);
    float cum = 0;
    for (auto kv: depm) {
        cum += kv.second;
        if (rn < cum) {
            return kv.first;
        }
    }
    return depm.back().first;
}

std::vector<std::pair<int,float>> parse_depths (std::string filename) {
    std::vector<std::pair<int,float>> depm;
    int cov;
    float prop;
    bool error;
    std::string input;
    std::ifstream fin;
    fin.open(filename);
    error=false;

    while (true) {
        std::getline(fin,input);
        if (!fin) break; //check for eof

        std::istringstream buffer(input);
        buffer >> cov >> prop;

        if (!buffer || !buffer.eof()) {
            error=true;
            break;
        }
        std::pair<int,float> pp(cov,prop);
        depm.push_back(pp);
    }
    if (error) {
        cout << "parsing issue\n" << std::endl;
    }
    return depm;
}

std::pair<int,int> generate_freq(size_t gen, float s, std::vector<std::pair<int,float>> depm) {
    float current_count = 1;
    float cur_popsize;
    for (float g = gen; g <= 25; g++) {
        float pre_popsize = pow(2,g-1);
        cur_popsize = pow(2,g);
        float sel;
        int nv;
        if (g >= 14) {
            //realize s as a penalty to 1, where 0 is neutral and small positive values = small negative selection. 1 is lethal. 
            sel = 1-s;
            float prob = (current_count * sel) / ((pre_popsize - current_count) + (current_count * sel));
            std::binomial_distribution<int> bd(cur_popsize, prob);
            nv = bd(generator); //* 2;
        } else {
            sel = 1;
            //before generation 14, we allow neither drift nor selection
            //as it is expected in development that there are 14 rounds of synchronized cell doubling with no lineage replacement
            nv = current_count * 2;
        }
        
        // std::cerr << current_count << '-';
        // std::cerr << g << "-";
        // std::cerr << sel << ",";
        current_count = nv;
        if (current_count == 0) {
            //if its gone down to 0, time to stop.
            //because its just gonna be 0 forever.
            break;
        }
    }
    //std::cerr << "\n";

    // std::cerr << "f:" << current_count/cur_popsize << ",";
    //return current_count/cur_popsize;
    //do some additional calculation. get a circle draw value conditioned on my real depth distribution I parsed in.
    int d = weighted_depth_choice(depm);
    if ((current_count == 0) || (current_count == cur_popsize)) {
        std::pair<int,int> fout(0,d);
        return fout;
    }
    // std::cerr << "d:" << d << ",";
    //simulate a circle draw
    int cc = 0;
    for (int i = 0; i < d; i++) {
        std::uniform_real_distribution<float> pg(0, 1);
        float p = pg(generator);
        //if we detect the alt, update the remaining distribution (sampling without replacement)
        if (p <= (current_count/cur_popsize)) {
            cc++;
            current_count -= 1;
            cur_popsize -= 1;
        }
        //if we run out of real individuals to draw (this will basically never happen, but...), break
        if ((cc >= current_count) || (cc == d)) {
            break;
        }
    }
    std::pair<int,int> o(cc,d);
    return o;
}

int main() {
    size_t permutations = 10000;

    std::vector<int> generations = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    std::vector<float> selv;
    selv.emplace_back(0.0);
    for (int pi = 1; pi <= 5; pi++) {
        for (float s = 1; s < 10; s++) {
            float sp = (s * pow(10,-pi));
            selv.emplace_back(sp);
        }
    }

    std::vector<std::pair<int,float>> depm = parse_depths("restricted_depth_values.tsv");

    std::ofstream out;
    out.open("smallsel_pis.tsv");
    out << "generation\tselection\tdepth\tcircle\n";
    std::cerr << "Beginning simulations.\n";
    size_t completed = 0;
    for (auto g: generations) {
        for (auto s: selv) {
            //out << g << "\t" << s << "\t";
            for (auto p = 0; p < permutations; p++) {
                std::pair<int,int> fd = generate_freq(g,s,depm);
                out << g << '\t' << s << '\t' << fd.second << '\t' << fd.first << '\n';
            }
            //out << '\n';
            completed++;
            std::cerr << "Simulations for gen:" << g << " selection:" << s << " complete. Total count:" << completed << "\n";
        }
    }
    out.close();
    return 0;
}
