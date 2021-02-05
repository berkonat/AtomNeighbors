#include <vector>
#include <string>

std::vector<double> read_radius(std::string fname);

std::vector<double> read_xyz(std::string fname,
                             std::vector<double>* box_lengths);

template<typename T>
std::vector<T> vecslice(std::vector<T> const &A, int start, int step);

template<typename T>
std::vector<T> vecslice(std::vector<T> const &A, int start, int step, int stride);

template<typename T1, typename T2>
std::vector<T1> gemv(std::vector<T1> const &A, std::vector<T2> const &x);

template<typename T1, typename T2>
std::vector<T1> gemv(std::vector<T1> const &A, std::vector<T2> const &B,
                     int ncolA, int* axis = NULL);

template<typename T>
double norm(std::vector<T> const &u);

template<typename T>
std::vector<double> norm(std::vector<T> const &A, int ncolA, int* axis = NULL);

template<typename T1, typename T2>
std::vector<T1> vsub(std::vector<T1> const &u, T2 const &v);

template<typename T1, typename T2>
std::vector<T1> vsub(std::vector<T1> const &u, std::vector<T2> const &v);

template<typename T1, typename T2>
std::vector<T1> vsub(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA);

template<typename T1, typename T2>
std::vector<T1> vadd(std::vector<T1> const &u, T2 const &v);

template<typename T1, typename T2>
std::vector<T1> vadd(std::vector<T1> const &u, std::vector<T2> const &v);

template<typename T1, typename T2>
std::vector<T1> vadd(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA);

template<typename T1, typename T2>
std::vector<T1> vmul(std::vector<T1> const &u, std::vector<T2> const &v);

template<typename T1, typename T2>
std::vector<T1> vmul(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA);

template<typename T1, typename T2>
std::vector<T1> vmul(std::vector<T1> const &u, T2 const &v);

template<typename T1, typename T2>
std::vector<T1> vdiv(std::vector<T1> const &u, std::vector<T2> const &v);

template<typename T1, typename T2>
std::vector<T1> vdiv(std::vector<T1> const &u, T2 const v);

template<typename T1, typename T2>
std::vector<T1> vdiv(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA);

