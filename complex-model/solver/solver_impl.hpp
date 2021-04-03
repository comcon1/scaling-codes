#ifndef PDESOLVERC1
#error You should not include solver_impl.hpp as a module.\
  It is just an implementation of the template class methods.
#endif


template <typename TParams> 
CNSolver<TParams>::CNSolver(std::ofstream &ofs): outFile(ofs) {
  defaults();
}

template <typename TParams> 
void CNSolver<TParams>::defaults() {
  timeStep = 0.01;  // dt
  stopAfter = 100; // T
  logEvery = 1; // log step
  reacLength = 1;  // l
  reacN = 500; // n
  nVars = 3; //TODO: move to constructor param
  pp.D = { 1, 1 };
  pp.bFlux = { dpair(0,0), 
               dpair(0,0) };

}

template <typename TParams> 
void CNSolver<TParams>::initialize() {

  setParams(pp);
  // calculated params
  reacStep = (double)reacLength/reacN; // h

  for (int i=0; i < nVars; ++i) {
    vars.push_back( dvector_type(reacN, 0.0) );
    vars_DUMP.push_back( dvector_type(reacN, 0.0) );
  }

  // lapplacian lpM
  lpM =  ubs::matrix<double,ubs::column_major>(reacN, reacN);
  for (int j = 0; j < reacN; ++j) {
    lpM(j,j) = -2;
    if (j < reacN-1)
      lpM(j,j+1) = 1;
    if (j > 0)
      lpM(j,j-1) = 1;
  }
  lpM(0,0) =-1;
  lpM(reacN-1,reacN-1) =-1;
  lpM *= 1/(reacStep*reacStep);

  // unit matrix unM
  unM = ubs::matrix<double,ubs::column_major> (reacN, reacN);
  for (int j = 0; j < reacN; ++j) { unM(j,j) = 1;}

  //creating matrix Lapp (plus) and Lapm (minus) for
  //the Crank-Nicolson's method 
  //  [E/dt-D*/2*Lap]*v2=[E/dt+D/2*Lap]*v1
  //Lapm*v2 = Lapp*v1
  //v2 =  v1 + J*dt
  for (int i = 0; i < nVars; ++i) {
    ubs::matrix<double> Lapp =  lpM*pp.D[i]/2 + (1/timeStep)*unM;
    ubs::matrix<double,ubs::column_major> Lapm = 
          -lpM*pp.D[i]/2 + (1/timeStep)*unM;

    // apply constant fluxes
    if (pp.bFlux[i].first != 0) {
        cout << 
          format("Flux is non-zero for variable %d for Left boundary") 
          % i << endl;
        Lapp(0,0) = +1;
        Lapp(0,1) = -1;
        Lapm(0,0) = -1;
        Lapm(0,1) = +1;
    }
    if (pp.bFlux[i].second != 0) {
        cout << 
          format("Flux is non-zero for variable %d for Right boundary") 
          % i << endl;
        Lapp(reacN-1,reacN-1) = -1; // check this signs!!
        Lapp(reacN-1,reacN-2) = +1;
        Lapm(reacN-1,reacN-1) = +1;
        Lapm(reacN-1,reacN-2) = -1;
    }
    // end apply constant fluxes

    auxLPMP.push_back(Lapp);


    dvector_type dl(reacN-1,0), du(reacN-1,0), d(reacN,0);
    for (int j=0;j<reacN;j++) {
        if (j>0) dl(j-1) = Lapm(j,j-1);
        d(j) = Lapm(j,j);
        if (j<reacN-1) du(j) = Lapm(j,j+1);
    }
    auxDL.push_back(dl); auxDL_DUMP.push_back(dl);
    auxDU.push_back(du); auxDU_DUMP.push_back(du);
    auxD.push_back(d);   auxD_DUMP.push_back(d);

  }

  // RK initialization

  for (int i=0; i<nVars; ++i) {
    RK_k1.push_back( dvector_type(reacN, 0.0) );
    RK_k2.push_back( dvector_type(reacN, 0.0) );
    RK_k3.push_back( dvector_type(reacN, 0.0) );
    RK_k4.push_back( dvector_type(reacN, 0.0) );
    RK_k0.push_back( dvector_type(reacN, 0.0) );
  }

  // temporaries reaction

  TMPR_vec_un = dvector_type(reacN, 1.0);
  TMPR_vec_t1 = dvector_type(reacN);
  TMPR_vec_t2 = dvector_type(reacN);
  TMPR_vec_t3 = dvector_type(reacN);
  TMPR_vec_t4 = dvector_type(reacN);
  TMPR_vec_t5 = dvector_type(reacN);
  TMPR_vec_t6 = dvector_type(reacN);

  TMPR_cpx_un  = cvector_type(reacN, dcomplex(1.0,0.0));
  TMPR_cpx_t1  = cvector_type(reacN);
  TMPR_cpx_t2  = cvector_type(reacN);
  TMPR_cpx_t3  = cvector_type(reacN);
  TMPR_cpx_t4  = cvector_type(reacN);
  TMPR_cpx_t5  = cvector_type(reacN);
  TMPR_cpx_t6  = cvector_type(reacN);
  TMPR_cpx_t7  = cvector_type(reacN);
  TMPR_cpx_t8  = cvector_type(reacN);
  TMPR_cpx_t9  = cvector_type(reacN);
  TMPR_cpx_t10 = cvector_type(reacN);
  TMPR_cpx_t11 = cvector_type(reacN);
  TMPR_cpx_t12 = cvector_type(reacN);

} // end CNSolver::initialize

/*
 * Make Runge-Kutta step.
 * vars_ is vector [nVars] of space vectors [reacN] 
 * t_ is time value (for non-autonomous systems)
 * Funcion
 */
template <typename TParams> 
void CNSolver<TParams>::performRKstep
      (std::vector<dvector_type> &vars_, double t_) {
  for (int i=0; i<nVars; ++i) 
    ubs::noalias(RK_k1[i]) = vars_[i];
  reaction(RK_k1, t_);
  for (int i=0; i<nVars; ++i) 
    ubs::noalias(RK_k2[i]) = vars_[i]+timeStep*RK_k1[i]/2;
  reaction(RK_k2, t_);
  for (int i=0; i<nVars; ++i) 
    ubs::noalias(RK_k3[i]) = vars_[i]+timeStep*RK_k2[i]/2;
  reaction(RK_k3, t_);
  for (int i=0; i<nVars; ++i) 
    ubs::noalias(RK_k4[i]) = vars_[i]+timeStep*RK_k3[i];
  reaction(RK_k4, t_);
  for (int i=0; i<nVars; ++i) 
    ubs::noalias(RK_k0[i]) = (RK_k1[i]+RK_k2[i]*2+RK_k3[i]*2+RK_k4[i])/6;

  for (int i=0; i<nVars; ++i) 
    vars_[i] += RK_k0[i]*timeStep;

} // end performRKstep


template <typename TParams> 
void CNSolver<TParams>::step
      (std::vector<ubs::vector<double> > &vars_, double &t_) {
    performRKstep(vars_, t_);

    for (int i=0; i<nVars; ++i) {
      ubs::noalias(vars_DUMP[i]) = ubs::prod(auxLPMP[i], vars_[i]);

      // apply non-zero flux
      if (pp.bFlux[i].first != 0) {
        vars_DUMP[i](0) += pp.bFlux[i].first/pp.D[i] * reacStep;
      }
      if (pp.bFlux[i].second != 0) {
        vars_DUMP[i](reacN-1) += pp.bFlux[i].second/pp.D[i] * reacStep;
      }
      // end apply non-zero flux

      auxD_DUMP[i] = auxD[i]; auxDL_DUMP[i] = auxDL[i]; auxDU_DUMP[i] = auxDU[i];
      boost::numeric::bindings::lapack::gtsv(reacN,
                auxDL_DUMP[i], auxD_DUMP[i], auxDU_DUMP[i], vars_DUMP[i]);
      vars_[i] = vars_DUMP[i];
    }

} // end CNSolver::step

template <typename TParams> 
void CNSolver<TParams>::runSolver() {
    setInitials(vars);
    //TODO: check if outFile is opened
    double t = 0;
    int co = 0;
    int logCo = 0;
    std::ostringstream os;
    while (t < stopAfter) {

        // logging
        os.clear();
        os.str("");
        if ( (abs(t) < 1e-10) || ( (int) (t / logEvery) > logCo) ) {

            os << format("%10.5f") % t;
            // do not print the last one
            for (int i=0; i<nVars; ++i) {
                for (auto x: vars[i])
                    os << format("%15.8f") % x;
            }
            // print spec_var as well
            for (auto x: spec_var)
                os << format("%15.8f") % x;

            outFile << os.str() << endl;

            cout << format("\n%10.5f") % t << flush;
        }

        step(vars, t);
        if (co % 500 == 0)
            cout << "." << flush;
        // step counter
        logCo = (int) (t / logEvery);
        t += timeStep;
        co++;
    }
} // end runSolver

dvector_type vector_realize (const cvector_type &c) {
  int sz = c.size();
  dvector_type result(c.size());
  for(int i=0; i<sz; ++i) 
    result(i) = c(i).real();
  return result;
}

cvector_type vector_complexize (const dvector_type &d) {
  int sz = d.size();
  cvector_type result(d.size());
  for(int i=0; i<sz; ++i) 
    result(i) = dcomplex(d(i), 0);
  return result;
}
