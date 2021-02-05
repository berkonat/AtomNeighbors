#include <vector>

std::vector<double> center_of_cell(std::vector<double> &box,
                                   double* tol = NULL);

std::vector<std::vector<int> > get_neighs(std::vector<double> &points, 
                                          double cutoff,
                                          std::vector<double>* box = NULL, 
                                          std::vector<std::vector<int> >* S = NULL, 
                                          std::vector<std::vector<double> >* D = NULL, 
                                          std::vector<std::vector<double> >* R = NULL, 
                                          std::vector<double>* ps = NULL,
                                          int* iterlimit = NULL);

