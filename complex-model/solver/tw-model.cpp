#include <iostream>
#include <vector>
#include <cctype>
#include <cmath>
#include <complex>

#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/program_options/value_semantic.hpp"

#include "boost/lexical_cast/try_lexical_convert.hpp"

using std::vector;
namespace po = boost::program_options;
using std::string;
typedef std::complex<double> dcomplex;


#include "solver.hpp"
#include "boost/numeric/ublasx/operation/sqrt.hpp"
#include "boost/numeric/ublasx/operation/tanh.hpp"
#include "boost/numeric/ublasx/operation/exp.hpp"
#include "boost/numeric/ublasx/operation/sqr.hpp"
#include "boost/numeric/ublasx/operation/element_pow.hpp"
#include "boost/numeric/ublasx/operation/sign.hpp"
namespace ubx = boost::numeric::ublasx;

namespace boost {  namespace program_options {

vector<po::option> ignore_numbers(std::vector<std::string>& args) {
        vector<po::option> result;
        int pos = 0;
        while(!args.empty()) {
            const auto& arg = args[0];
            double num;
            if(boost::conversion::try_lexical_convert(arg, num)) {
                result.push_back(po::option());
                po::option& opt = result.back();

                opt.position_key = pos++;
                opt.value.push_back(arg);
                opt.original_tokens.push_back(arg);

                args.erase(args.begin());
            } else {
                break;
            }
        }

        return result;
    }
} }

struct MyParams {
  double bmp, Kd, sb, sa, na, nx, nm, mmp3, nu, h, smad_threshold, 
         ca, cx, tb, tld, de_init, caf_a, caf_x, caf_threshold,
        smad_tau, caf_tau, smad_mean ;
  std::vector<double> D;
  std::vector<dpair> bFlux;
};

class MySolver: public CNSolver<MyParams> {
    vector<double> cmdParams;
  public:
    MySolver(std::ofstream &of_, vector<double> cmd_): CNSolver<MyParams>(of_) { cmdParams = cmd_; }

  protected:

    virtual void reaction(std::vector<ubs::vector<double> > &vars_, double t_) {
      auto _c = vars_[0], _n = vars_[1], caf = vars_[2], de = vars_[3] ;
      auto &_1 = TMPR_vec_un,
           &_bmpa = TMPR_vec_t2,
           &_smad = TMPR_vec_t3,
           &_tmp = TMPR_vec_t4;
           /* temporaries for bmpa calculations */
      auto &u1 = TMPR_cpx_un, 
           &c = TMPR_cpx_t1,
           &n = TMPR_cpx_t2,
           &_t0 = TMPR_cpx_t3,
           &c3 = TMPR_cpx_t4, 
           &c2n = TMPR_cpx_t5,
           &cn2 = TMPR_cpx_t6,
           &c2 = TMPR_cpx_t7,
           &n2 = TMPR_cpx_t8,
           &cn = TMPR_cpx_t9,
           &n3 = TMPR_cpx_t10,
           &res_y = TMPR_cpx_t11;
      c = vector_complexize(_c); n = vector_complexize(_n);
      c3  = ubs::element_prod(c, ubx::sqr(c) );
      c2n = ubs::element_prod( ubx::sqr(c), n);
      cn2 = ubs::element_prod( ubx::sqr(n), c);
      c2  = ubx::sqr(c);
      n2  = ubx::sqr(n);
      cn  = ubs::element_prod(c,n);
      n3  = ubs::element_prod(n, ubx::sqr(n) );

    _t0 =  ubx::element_pow( 2.6968322100701824e16*u1 + 5.6595440196711e15*c + 3.95902006695e14*c2  
                            + 9.23150925e12*c3 - 7.98154665879636e16*n - 2.4131419699725e15*cn 
                            + 2.21764138875e14*c2n + 8.050040609832e16*n2 - 2.304054763875e15*cn2 
                            - 2.792943963425e16*n3 + 
                            ubx::sqrt(
                                4.*ubx::element_pow(
                                  - 5.6652092289e10*u1 - 7.9259661e9*c - 2.772225e8*c2 + 1.117783836e11*n 
                                  - 4.4397225e9*cn - 5.86147225e10*n2, 3)
                                + ubx::sqr(
                                    2.6968322100701824e16*u1 + 5.6595440196711e15*c + 3.95902006695e14*c2 
                                  + 9.23150925e12*c3 - 7.98154665879636e16*n - 2.4131419699725e15*cn 
                                  + 2.21764138875e14*c2n + 8.050040609832e16*n2 - 2.304054763875e15*cn2 
                                  - 2.792943963425e16*n3 ) 
                                      ), 
                            (1./3) );
 
      res_y = -1.4284693950432112e-6*(-238017.*u1 - 16650.*c - 483350.*n) 
               + dcomplex(8.998793299727685e-7,1.5586367201938738e-6)* 
                  ubs::element_div(-5.6652092289e10*u1 - 7.9259661e9*c - 2.772225e8*c2 
                          + 1.117783836e11*n - 4.4397225e9*cn - 5.86147225e10*n2, _t0)
               - dcomplex(5.668884550989927e-7,-9.818796064556836e-7)*_t0;

      _bmpa = 
        
    ubs::element_prod( _1 - ubx::sign(_n - _1 * 1e-6),  0.5 * 
        // solution for zero noggin
        (ubx::sqrt(_1 * 0.3 + ubx::sqr(_1 * 0.7 - _c)/4.) + 0.5*(_1 * 0.7 - _c))
                     ) +
    ubs::element_prod( _1 + ubx::sign(_n - _1 * 1e-6) , 0.5 * 
        // solution for non-zero noggin
        vector_realize(
          u1 - ubs::element_div(3.33*ubs::element_prod(c,res_y), (50.*n-46.67*res_y) ) - res_y) 
                     );


      spec_var = _bmpa;

      _smad = ubs::element_div ( _1, 
                                ( ubs::element_div( _1, ubx::sqr( _bmpa) ) 
                                   + _1 / (pp.sb*pp.sb) )
                               ) * (pp.sa/pp.sb/pp.sb);

      // chordin
      vars_[0] = ( _1 - ubx::tanh(de*5) ) * (pp.ca/2)
                                -
                            _c * pp.cx * pp.tld;

      // noggin2
      vars_[1] = ( _1 - ubx::tanh(de*5) ) * (pp.na/2) 
                                - 
                  _n * ( (pp.nm*pp.mmp3 + 1)*pp.nx )  ;

      // caf
      vars_[2] = ( _1 - ubx::tanh(de*5) ) * (pp.caf_a/2)
                                -
                            caf * pp.caf_x;

      // delta
      vars_[3] = 
        ( de * pp.h - ubs::element_prod(de, ubs::element_prod( de, de) ) 
          + 
         (_1 * pp.caf_threshold - caf) *
              std::exp( -std::pow(t_/pp.caf_tau, 2) )
          +
         (_smad - _1 * pp.smad_threshold) *
              std::exp( -0.5 * std::pow( (t_-pp.smad_mean)/pp.smad_tau, 2) ) 
        ) * pp.nu;

//      cout << " --> " << vars_[1] << endl;

    } // end CNSolver::reaction

    virtual void setParams(MyParams &pp_) {
      // Parameters for setting up
      timeStep = 8.0;  // not defined in paper
      stopAfter = cmdParams[0] * 3600;   // T
      logEvery = 1000; // my choice
      reacLength = cmdParams[1];  // L
      reacN = int(reacLength / 3.0); // not defined in paper
      nVars = 4;

      pp_.na = cmdParams[2];        // 1e-4
      pp_.nx = cmdParams[3];        // 5e-5
      pp_.mmp3 = cmdParams[4];      // 1
      pp_.sa = cmdParams[5];        // 8
      pp_.sb = cmdParams[6];        // 1
      pp_.smad_threshold            
             = cmdParams[7];        // 2.0
      pp_.nu = cmdParams[8];        // 1.0e-4
      pp_.h  = cmdParams[9];        // 1.5
      pp_.tb = cmdParams[10];       // ?
      pp_.ca = cmdParams[11];       // 1.0e-3
      pp_.cx = cmdParams[12];       // 1.0e-3
      pp_.de_init = cmdParams[13];  // 400 vs 200
      // ---- ---- ----
      pp_.caf_a = cmdParams[14];      // 1e-3 
      pp_.caf_x = cmdParams[15];      // 1e-3
      pp_.caf_threshold = cmdParams[16];      // 0.6
      // timing params. Strictly set up
      pp_.smad_tau  = 1 * 3600; //
      pp_.caf_tau   = 1 * 3600;  
      pp_.smad_mean = 2 * 3600;

      pp_.nm = 3;
      pp_.D = { 15.0, 0.5, 0, 0};
      pp_.bFlux = { dpair( 0 , 0 ), dpair( 0 , 0 ), dpair( 0 , 0 ), dpair( 0 , 0 ) };
      
      pp_.tld = 1/(1 + pp_.mmp3*pp_.mmp3/pp_.tb/pp_.tb); // derivative parameter

      cout << "\n ** P A R A M E T E R S ** " << endl;
      cout << format("BMP4: %5.2f  Kc: %5.2f  Kn: %5.2f") 
        %  1.0 % 0.3 % 0.1;
      cout << endl;
      cout << format("na: %6.1e  nx: %6.1e  ca: %6.1e  cx: %6.1e") 
        % pp_.na % pp_.nx % pp_.ca % pp_.cx;
      cout << endl;
      cout << format("nm: %5.2f mmp3: %5.2f t_b: %5.2f") 
        % pp_.nm % pp_.mmp3 % pp_.tb;
      cout << endl;
      cout << format("nu: %6.1e h: %5.2f Sa: %5.2f Sb: %5.2f SMAD_thr: %5.2f") 
        % pp_.nu % pp_.h % pp_.sa % pp_.sb % pp_.smad_threshold;
      cout << endl;
      cout << format("ncoeff: %5.2f tau_1: %5.2f mu_2: %5.2f tau_2: %5.2f")
       % 0.84 % 1 % 2 % 1;
      cout << endl;
      cout << format("De Init: %7.1f Diff. coeff Chd/N2: %5.2f / %5.2f") 
        % pp_.de_init % pp_.D[0] % pp_.D[1];
      cout << endl;
      cout << format("CAF_a: %6.1e CAF_x: %6.1e CAF_thr: %5.2f") 
        % pp_.caf_a % pp_.caf_x % pp_.caf_threshold;
      cout << endl;
      cout << " ============================= " << endl;
      // TODO: output
      
    }

    virtual void setInitials(std::vector<ubs::vector<double> > &vars_) {
//        vars_[2] = TMPR_vec_un * 0; // <1...1>*0
        double aux = 0;

        for(int j=0; j<reacN; j++) {
          aux = ((double)reacLength*j/reacN)/pp.de_init;
          vars_[0](j) = (pp.ca/pp.cx)*exp(-aux*aux/2.);         // Chd
          vars_[1](j) = (pp.na/pp.nx)*exp(-aux*aux/2.);         // N2
          vars_[2](j) = (pp.caf_a/pp.caf_x)*exp(-aux*aux/2.);   // Caf
          vars_[3](j) = -1.0 + ( (aux > 2.0) ? 2.0 : aux );     // De
        }
        // initialize debug vector
        spec_var = dvector_type(reacN, 0.0);
    }

};

/// DO NOT CALL VIRTUAL FROM CONSTRUCTOR!!


int main(int ac, char **av) {

  // defining parameters of command line call
  //
  po::options_description desc("Allowed options");
  
  desc.add_options()
      ("help,h", "produce help message")
      ("output-file,O", po::value<string>(), "set output file")
      ("input-params", po::value< vector<double> >(), 
          "input parameters: D, not defined")
  ;

  po::positional_options_description pos;
  pos.add("input-params", -1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(ac, av)
              .extra_style_parser(&po::ignore_numbers)
              .options(desc)
              .positional(pos)
              .run(), vm);
    po::notify(vm);
  } catch(...) {
    cout << "Illegal parameter string... Use -h for help." << endl;
    return 1;
  }

  if (vm.count("help")) {
    cout << desc << "\n";
    cout << "TODO: parameter description" << endl;
    return 1;
  }

  if (! vm.count("output-file") ) {
    cout << "Output file missing. Run with -h for help." << endl;
    return 1;
  }

  if (! vm.count("input-params") ) {
    cout << "Missing model parameters. Run with -h for help." << endl;
    return 1;
  }

  string ofName = vm["output-file"].as<string>();

  vector<double> rParams = vm["input-params"].as<vector<double> >();

  if (rParams.size() != 17) {
    cout << "Incorrect number of parameters. See help carefully." << endl;
    return 1;
  } else {
      for (auto i: rParams)
          cout << i << " ";
      cout << endl;
  }

  std::ofstream of( ofName, std::ios::out);

  of << "# write by main.x" << endl;
  MySolver s(of, rParams);
  s.initialize();
  s.runSolver();

  of.close();

  return 0;
}
