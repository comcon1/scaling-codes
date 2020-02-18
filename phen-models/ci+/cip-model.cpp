// чувак, попробуй, чтобы морфоген хуярил контрактор в зоне своей экспрессии
// тогда контрактор будет развиваться комплементарно морфогену
// что в целом сделает паттерн более заебись


#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include <iostream>
#include <vector>

using std::vector;
namespace po = boost::program_options;
using std::string;


#include "solver.hpp"
#include <cctype>
#include <cmath>

struct MyParams {
  double am1, am2, am3, ae, ae2, be, t0;
  int h;
  std::vector<double> D;
  std::vector<dpair> bFlux;
};

class MySolver: public CNSolver<MyParams> {
    vector<double> cmdParams;
  public:
    MySolver(std::ofstream &of_, vector<double> cmd_): CNSolver<MyParams>(of_) { cmdParams = cmd_; }

  protected:

    virtual void reaction(std::vector<ubs::vector<double> > &vars_) {
      auto m = vars_[0], e = vars_[1];
      auto &_un = TMPR_vec_un,
           &_mh = TMPR_vec_t1;

      // mh = (m/t0)^h
      _mh = m / pp.t0;
      for (int i=1; i<pp.h; ++i)
           _mh = ubs::element_prod( _mh, m) / pp.t0;
      
//      cout << _mh << endl;

      vars_[0] =  - ubs::element_prod( 
             ubs::element_div(m * pp.am2, _un+ubs::element_prod(m,m)*pp.am3), 
             ubs::element_div( ubs::element_prod(e,e) , _un*pp.am1*pp.am1+ubs::element_prod(e,e)) 
                                     );
      
            //   ubs::element_div ( m * pp.am1, _un + e ) * -1.0
            // - ubs::element_div ( ubs::element_prod(m, m) * pp.am2, _un + e ) ;

      vars_[1] = 
              ubs::element_div ( _mh * pp.be, _un  + _mh ) 
            - ubs::element_div ( e * pp.ae , _un*pp.ae2 + _mh ); 
      
//      cout << " --> " << vars_[1] << endl;

    } // end CNSolver::reaction

    virtual void setParams(MyParams &pp_) {
      // Parameters for setting up
      timeStep = 0.1;  // not defined in paper
      stopAfter = cmdParams[12]; // T
      logEvery = 1000; // my choice
      reacLength = cmdParams[11];  // 
      reacN = int(reacLength / 1.5); // not defined in paper
      nVars = 2;

      pp_.am1 = cmdParams[2];
      pp_.am2 = cmdParams[3];
      pp_.am3 = cmdParams[4];
      pp_.ae = cmdParams[5];
      pp_.ae2 = cmdParams[6];
      pp_.be  = cmdParams[8];
      pp_.t0   = cmdParams[9];
      pp_.h  = (int) (cmdParams[10]);

      
      pp_.D = { cmdParams[0], cmdParams[1]};
      pp_.bFlux = { dpair( -cmdParams[7], 0), dpair( 0 , 0 ) };

      cout << format("Diff. coeff 1: %10.5f") % pp_.D[0];
      cout << endl;
      cout << format("Diff. coeff 2: %10.5f") % pp_.D[1];
      cout << endl;
      cout << format("am1 = %10.5f") % pp_.am1;
      cout << endl;
      cout << format("am2 = %10.5f") % pp_.am2;
      cout << endl;
      cout << format("am3 = %10.5f") % pp_.am3;
      cout << endl;      
      cout << format("ae = %10.5f") % pp_.ae;
      cout << endl;
      cout << format("be = %10.5f") % pp_.be;
      cout << endl;
      cout << format("t0 = %10.5f") % pp_.t0;
      cout << endl;
      
    }

    virtual void setInitials(std::vector<ubs::vector<double> > &vars_) {
        for(int j=0; j<reacN; j++) {
          vars_[0](j) = 0.0;  // M
          vars_[1](j) = 0.0;  // E
        }
    }

};

/// DO NOT CALL VIRTUAL FROM CONSTRUCTOR!!


int main(int ac, char **av) {

  // defining parameters of command line call
  //
  po::options_description desc("Allowed options");
  
  desc.add_options()
      ("help", "produce help message")
      ("output-file,O", po::value<string>(), "set output file")
      ("input-params", po::value< vector<double> >(), 
          "input parameters: Dm, De, am1, am2, am3, ae, mFlux, be, t0, h, L, T")
  ;

  po::positional_options_description pos;
  pos.add("input-params", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).
            options(desc).positional(pos).run(), vm);
  po::notify(vm);

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

  if (rParams.size() != 13) {
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
