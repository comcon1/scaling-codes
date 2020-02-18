#ifndef PDESOLVERC1
#define PDESOLVERC1

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>

#include <boost/format.hpp>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/bindings/detail/property_map.hpp>
#include <boost/numeric/bindings/views.hpp>
#include <boost/numeric/bindings/ublas/vector_expression.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack.hpp>

using std::cout;
using std::endl;
using std::flush;
namespace ubs = boost::numeric::ublas;
using boost::format;

typedef std::pair<double,double> dpair;

template <typename TParams> class CNSolver {
  private:
    // parameters
    std::ofstream &outFile;
    std::vector<ubs::vector<double> > vars;
    ubs::matrix<double,ubs::column_major> lpM,
                                          unM;
    std::vector<ubs::matrix<double> > auxLPMP; 
    std::vector<ubs::vector<double> > auxDL, auxDU, auxD;

    std::vector<ubs::vector<double> > vars_DUMP;
    std::vector<ubs::vector<double> > auxDL_DUMP, auxDU_DUMP, auxD_DUMP;

    std::vector<ubs::vector<double> > RK_k1, RK_k2, RK_k3, RK_k4, RK_k0;

    // do not touch it
    void defaults();
    void performRKstep(std::vector<ubs::vector<double> > &);
    void step(std::vector<ubs::vector<double> > &, double &); 

  public:
    CNSolver(std::ofstream &);
    void runSolver();
    void initialize();

  protected:
    double stopAfter;
    double timeStep;
    double logEvery;
    double reacLength;
    int    reacN;
    double reacStep;
    int    nVars;

    virtual void setInitials(std::vector<ubs::vector<double> >&) = 0;

    TParams pp;
    virtual void setParams(TParams &) = 0;

    ubs::vector<double> TMPR_vec_un, TMPR_vec_t1, TMPR_vec_t2, TMPR_vec_t3;
    virtual void reaction(std::vector<ubs::vector<double> > &) = 0;
};

#include "solver_impl.hpp"

#endif
