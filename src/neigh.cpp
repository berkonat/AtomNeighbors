#include "neigh.hpp"
#include "utils.hpp"
#include <iostream>

std::vector<double> center_of_cell(std::vector<double> &box,
                                   double* tol)
{
    if(tol == NULL) *tol = 1.0E-6;
    std::vector<double> center = {0.5,0.5,0.5};
    return gemv(box, center);
}

std::vector<std::vector<int> > get_neighs(std::vector<double> &points, 
                                          double cutoff,
                                          std::vector<double>* box, 
                                          std::vector<std::vector<int> >* S, 
                                          std::vector<std::vector<double> >* D, 
                                          std::vector<std::vector<double> >* R, 
                                          std::vector<double>* ps,
                                          int* iterlimit)
{
    std::vector<double> vzero(3, 0.0);
    bool rtn_vecs = (R != NULL) ? true : false;
    bool rtn_dists = (D != NULL) ? true : false;
    bool rtn_shifts = (S != NULL) ? true : false;
    //int itlimit = 10;
    bool nops = false;
    bool havebox = false;
    std::vector<int> img_boxes;
    std::vector<int> img_box;
    std::vector<double> box_len(3, 0.0);
    std::vector<double> vbox(3, 0.0);
    std::vector<int> shift(3, 0);

    int d, i, j, k;

    if (box != NULL){
        havebox = true;
        for(d=0;d<3;d++) {
            for(i=0;i<3;i++) vbox[i] = (*box)[d*3+i];
            box_len[d] = norm((*box));
            shift[d] = int(cutoff / box_len[d]);
        }
    }
    std::cout << " shift: [" << shift[0] << "," << shift[1] << "," << shift[2]  << "]" << std::endl;
    
    //if(iterlimit != NULL) itlimit = *iterlimit;
    
    if (ps == NULL){
        ps = &points;
        nops = true;
    } 
    int NP, MP;
    NP = int(ps->size()/3);
    MP = int(points.size()/3);

    std::vector<std::vector<int> > neigh_list(NP);
    if (rtn_dists) { 
        (*D).resize(NP);
        //for(i=0;i<N;i++)
        //    (*D)[i].resize(M);
    }
    if (rtn_shifts) {
        (*S).resize(NP);
        //for(i=0;i<N;i++)
        //    (*S)[i].resize(M);
    }
    if (rtn_vecs) {
        (*R).resize(NP);
        //for(i=0;i<N;i++)
        //    (*R)[i].resize(M);
    }
        
    for(i=-shift[0];i<=shift[0];i++)
        for(j=-shift[1];j<=shift[1];j++)
            for(k=-shift[2];k<=shift[2];k++)
                img_boxes.insert(img_boxes.end(), {i,j,k});
                //std::cout << " max imgs: " << i << "," << j << "," << k  << std::endl;
    
    std::vector<std::vector<double> > box_pp(MP);
    if (havebox){
        std::vector<double> box_vecs;
        box_vecs = gemv(*box, img_boxes, 3);
        int NB = box_vecs.size()/3;

        #pragma omp parallel for schedule(guided)
        for(int j=0; j<MP; j++){
            std::vector<double> vec;
            std::vector<double> pv(points.cbegin() + j, points.cbegin() + j + 3);
            for(int b=0;b<NB;b++){
                std::vector<double> bv(box_vecs.cbegin() + b, box_vecs.cbegin() + b + 3);
                vec = vadd(pv, bv);
                box_pp[j].insert(box_pp[j].end(), vec.begin(), vec.end()); 
            }
        }
    }

    #pragma omp parallel for schedule(guided)
    for(int p=0; p<NP; p++){
        std::vector<double> dist;
        std::vector<double> vec;
        std::vector<int> local_J;
        std::vector<double> local_D;
        std::vector<int> local_S;
        std::vector<double> local_R;
        std::vector<double> pp((*ps).cbegin() + p, (*ps).cbegin() + p + 3);

        int i, b;//, nD;
        if (havebox){
            int k = nops ? p : 0;
            for(int j=k; j<MP; j++){
                vec = vsub(box_pp[j], pp, 3);
                dist = norm(vec, 3);
                for(b=0; b<int(dist.size()); b++){
                    if (dist[b] < cutoff && dist[b] > 0.0) { 
                        local_J.push_back(j);
                        if (rtn_dists) local_D.push_back(dist[b]);
                        if (rtn_shifts) { 
                            for(i=0;i<3;i++) img_box[i] = img_boxes[b*3 + i];
                            local_S.insert(local_S.end(), img_box.begin(), img_box.end());
                        }
                        if (rtn_vecs) local_R.insert(local_R.end(), vec.begin(), vec.end());
                    }
                }
            }
        } else {
            vec = vsub(points, pp, 3);
            dist = norm(vec, 3);
            for(b=0; b<int(dist.size()); b++){
                if (dist[b] < cutoff && dist[b] > 0.0) { 
                    local_J.push_back(j);
                    if (rtn_dists) local_D.push_back(dist[b]);
                    if (rtn_shifts) { 
                        for(i=0;i<3;i++) img_box[i] = img_boxes[b*3 + i];
                        local_S.insert(local_S.end(), img_box.begin(), img_box.end());
                    }
                    if (rtn_vecs) local_R.insert(local_R.end(), vec.begin(), vec.end());
                }
            }
        }
        neigh_list[p] = local_J;
        if (rtn_shifts) (*S)[p] = local_S; 
        if (rtn_dists) (*D)[p] = local_D;
        if (rtn_vecs) (*R)[p] = local_R;
    }
    return neigh_list;
}

