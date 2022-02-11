#ifndef KinematicParametersError_H
#define KinematicParametersError_H
#include "signal.h"

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"

#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToCurvilinear.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h"

/**
 * Class to store the error matrix
 * for (x,y,z,p_x,p_y,p_z,m)
 * particle parametrization
 *
 * Kirill Prokofiev January 2003
 */


class KinematicParametersError{

public:
 KinematicParametersError()
 {vl = false;}

 KinematicParametersError(const AlgebraicSymMatrix77& er):
                             theCovMatrix(er)
 {vl = true;

    std::cerr << "KinematicParametersError = " << er << '\n';
    raise(SIGINT);
    if (fabs(er(0,0) - 0.00041556) < 1e-8) raise(SIGINT);
    for (int i = 0; i < 7; i++) {
        if (er(i,i) < 0) {
            std::cerr << "er(" << i << "," << i << ") < 0!\n";
            exit(1);
        }
    }
}
 
  KinematicParametersError(const CartesianTrajectoryError& err, float merr) {
    theCovMatrix.Place_at(err.matrix(),0,0);
    theCovMatrix(6,6) = merr * merr;
    vl = true;
    if (fabs(theCovMatrix(0,0) - 0.00041556) < 1e-8) raise(SIGINT);
    std::cerr << "KinematicParametersError = " << theCovMatrix << '\n';
    raise(SIGINT);
  }
  
/**
 * access methods
 */ 
 
 AlgebraicSymMatrix77 const & matrix() const {
    std::cerr << "matrix() = " << theCovMatrix << '\n';
    return theCovMatrix;
}

  //AlgebraicSymMatrix77  & matrix() {
  AlgebraicSymMatrix77  matrix() {
    std::cerr << "matrix() = " << theCovMatrix << '\n';
    return theCovMatrix;
}
 
 
 bool isValid() const
 {return vl;}

private:
 AlgebraicSymMatrix77 theCovMatrix;
 bool vl;
};
#endif

