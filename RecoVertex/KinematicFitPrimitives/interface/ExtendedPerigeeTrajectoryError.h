#ifndef ExtendedPerigeeTrajectoryError_H
#define ExtendedPerigeeTrajectoryError_H

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class ExtendedPerigeeTrajectoryError
{ 
public:
 ExtendedPerigeeTrajectoryError(): weightAvailable(false),vl(false)
 {
    std::cerr << "called without weights!\n";
}

 ExtendedPerigeeTrajectoryError(const AlgebraicSymMatrix66& covariance):
                               cov(covariance),weightAvailable(false),
			       vl(true)
 {

    std::cerr << "called with weights!\n";
}


/**
 * Access methods
 */

 bool isValid() const
 {return vl;}

 bool weightIsAvailable() const
 {return weightAvailable;}

 const AlgebraicSymMatrix66 & covarianceMatrix()const
 {return cov;}
 
 const AlgebraicSymMatrix66 & weightMatrix(int & error)const
 {
  error = 0;
  if(! weightIsAvailable()) {
    weight = cov.Inverse(error);
   if(error != 0) {
        LogDebug("RecoVertex/ExtendedPerigeeTrajectoryError") << "unable to invert covariance matrix\n";
        fprintf(stderr, "unable to invert covariance matrix!\n");
    }
   weightAvailable = true;
  }
    std::cerr << "cov = " << cov << '\n';
    std::cerr << "weight = " << weight << '\n';
    if (cov(1,1) < 0) {
        std::cerr << "cov(1,1) < 0\n";
        exit(1);
    }
    if (cov(2,2) < 0) {
        std::cerr << "cov(2,2) < 0\n";
        exit(1);
    }
    if (weight(1,1) < 0) {
        std::cerr << "weight(1,1) < 0\n";
        exit(1);
    }
    if (weight(2,2) < 0) {
        std::cerr << "weight(2,2) < 0\n";
        exit(1);
    }
  return weight;
 }
 
private:
 AlgebraicSymMatrix66 cov;
 mutable AlgebraicSymMatrix66 weight;
 mutable bool weightAvailable;
 mutable bool vl;
};
#endif
