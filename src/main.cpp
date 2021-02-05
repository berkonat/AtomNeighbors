#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

#include <chrono>

#include "utils.hpp"
#include "neigh.hpp"

int main(int argc, const char* argv[])
{
    int i, j;//, k, d, id;
    std::vector<double> pos;
    std::vector<double> box;
    std::vector<double> box_lengths;
    std::vector<double> rads;
    std::vector<std::vector<int> > neigh_list;
    pos = read_xyz("a1t.dat", &box_lengths);
    rads = read_radius("r1t.dat"); 

    //for (auto i = rads.begin(); i != rads.end(); ++i)
    //    std::cout << *i << std::endl;

    for (i = 0; i<int(box_lengths.size()); i++) {
        std::vector<double> u(3, 0.0);
        u[i] = box_lengths[i];
        box.insert(box.end(), u.begin(), u.end());
    }

    for (i = 0; i<3; i++) {
        for (j = 0; j<3; j++)
            std::cout << box[i*3+j] << " ";
        std::cout << std::endl;
    }
   
    std::cout << " Number of points: " << pos.size()/3 << std::endl;
    std::cout << " Number of raduises: " << rads.size() << std::endl;

    // Measure time to calculate neighbours
    auto start_time = std::chrono::high_resolution_clock::now();
    //neigh_list = get_neighs(pos, pos, 6.0, &box);
    neigh_list = get_neighs(pos, 6.0, &box);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    //std::cout << " This took " << elapsed_seconds.count() << " s." << std::endl;
    //std::cout << " This took " << (end_time - start_time)/std::chrono::milliseconds(1) << " s." << std::endl;
    
    //std::cout << " Number of neigh_list: " << neigh_list.size() << std::endl;
    int nneighs = neigh_list.size();
    //#pragma omp parallel for reduction(+:nneighs)
    //for (int i = 0; i < neigh_list.size(); ++i) 
    //    nneighs += neigh_list[i].size();

    std::cout << " Number of total neighs: " << nneighs << std::endl;

    /*
    for (i = 0; i < neigh_list.size(); ++i){
        //std::cout << " Number of neigh for " << i << " : " << neigh_list[i].size() << std::endl;
        for (j = 0; j < neigh_list[i].size(); ++j){
            id = int(neigh_list[i][j][0]);
            std::cout << " atom " << id << " pos: ";
            for (d = 0; d < 3;++d) std::cout << pos[id][d] << " ";
            std::cout << " neighs : ";
            for (k = 0; k < neigh_list[i][j].size(); ++k)
                std::cout << neigh_list[i][j][k] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    */
    //for (auto i = neigh_list.begin(); i != neigh_list.end(); ++i)
    //    for (auto j = (*i).begin(); j != (*i).end(); ++j)
    //        std::cout << *j << std::endl;
    std::cout << " This took " << elapsed_seconds.count() << " s." << std::endl;
    return 0;
}

