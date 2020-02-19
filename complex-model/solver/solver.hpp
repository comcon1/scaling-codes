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
typedef ubs::vector<double> dvector_type;
typedef ubs::vector<std::complex<double> > cvector_type;

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
    void performRKstep(std::vector<ubs::vector<double> > &, double = 0);
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


    dvector_type spec_var;
    dvector_type TMPR_vec_un, 
      TMPR_vec_t1, TMPR_vec_t2, TMPR_vec_t3, TMPR_vec_t4;
    cvector_type TMPR_cpx_un, 
      TMPR_cpx_t1, TMPR_cpx_t2, TMPR_cpx_t3, TMPR_cpx_t4,
      TMPR_cpx_t5, TMPR_cpx_t6, TMPR_cpx_t7, TMPR_cpx_t8,
      TMPR_cpx_t9, TMPR_cpx_t10, TMPR_cpx_t11, TMPR_cpx_t12,
      TMPR_cpx_t13, TMPR_cpx_t14, TMPR_cpx_t15, TMPR_cpx_t16;
    virtual void reaction(std::vector<ubs::vector<double> > &, double = 0) = 0;

};

    
dvector_type vector_realize(const cvector_type&);
cvector_type vector_complexize(const dvector_type&);

#include "solver_impl.hpp"

#endif
