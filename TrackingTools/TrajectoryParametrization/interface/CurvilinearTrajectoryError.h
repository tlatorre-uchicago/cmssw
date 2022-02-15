#ifndef _TRACKER_CURVILINEARTRAJECTORYERROR_H_
#define _TRACKER_CURVILINEARTRAJECTORYERROR_H_

#include "DataFormats/GeometrySurface/interface/TrivialError.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "DataFormats/Math/interface/Error.h"
#include "TVectorD.h"

/* Compute the V = LDL^T factorization of the Matrix V.
 *
 * Taken from https://github.com/nojhan/cholesky/blob/master/cholesky.h. */
int ldl(TMatrixD &V, TMatrixD &D)
{
    int i, j, k;
    int N = V.GetNrows();
    // example of an invertible matrix whose decomposition is undefined
    assert(V(0,0) != 0); 
    TMatrixD L(N,N);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            D(i,j) = 0;
            V(i,j) = 0;
        }
    }

    D(0,0) = V(0,0);

    for (j = 0; j < N; j++) {
        L(j, j) = 1;

        D(j,j) = V(j,j);
        for (k = 0; k <= j-1; k++) {
            D(j,j) -= L(j,k) * L(j,k) * D(k,k);
        }

        for (i = j+1; i < N; i++) {
            for (k = 0; k <= j-1; k++) {
                L(i,j) -= L(i,k)*L(j,k) * D(k,k);
            }
            L(i,j) /= D(j,j);
        }
    }
    //this->_L = root(L, D);
    return 0;
}


/* Compute the final symetric matrix: _L = L D^1/2
 * remember that V = ( L D^1/2) ( L D^1/2)^T
 * the factorization is thus L*D^1/2 */
TMatrixD root(TMatrixD& L, TMatrixD &D)
{
    int i;

    // fortunately, the square root of a diagonal matrix is the square 
    // root of all its elements
    TMatrixD sqrt_D = D;
    for (i = 0; i < D.GetNrows(); i++) {
        sqrt_D(i,i) = sqrt(D(i,i));
    }

    return L*sqrt_D;
}

bool isHermitian(TMatrixD A)
{
    int i, j;

    if (A.GetNrows() != A.GetNcols()) {
        fprintf(stderr, "error: A is not square!\n");
        return false;
    }

    for (i = 0; i < A.GetNrows(); i++) {
        for (j = 0; j < A.GetNcols(); j++) {
            if (A(i,j) != A(j,i))
                return false;
        }
    }

    return true;
}

TMatrixD eye(int N)
{
    int i, j;

    TMatrixD D(N,N);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i == j)
                D(i,j) = 1;
            else
                D(i,j) = 0;
        }
    }

    return D;
}

TMatrixD eig(TMatrixD &V, TMatrixD &U)
{
    int i;
    TMatrixD D = eye(V.GetNrows());
    TVectorD e(V.GetNrows());

    U = V.EigenVectors(e);
    for (i = 0; i < V.GetNrows(); i++)
        D(i,i) = e(i);

    return D;
}

/* modchol_ldlt  Modified Cholesky algorithm based on LDL' factorization.
 * [L D,P,D0] = modchol_ldlt(A,delta) computes a modified Cholesky
 * factorization P*(A + E)*P' = L*D*L', where P is a permutation matrix, L is
 * unit lower triangular, and D is block diagonal and positive definite with
 * 1-by-1 and 2-by-2 diagonal blocks.  Thus A+E is symmetric positive definite,
 * but E is not explicitly computed.  Also returned is a block diagonal D0 such
 * that P*A*P' = L*D0*L'.  If A is sufficiently positive definite then E = 0
 * and D = D0.  The algorithm sets the smallest eigenvalue of D to the
 * tolerance delta, which defaults to sqrt(eps)*norm(A,'fro').  The LDL'
 * factorization is compute using a symmetric form of rook pivoting proposed by
 * Ashcraft, Grimes and Lewis.
 * 
 * Reference:
 * S. H. Cheng and N. J. Higham. A modified Cholesky algorithm based on a
 * symmetric indefinite factorization. SIAM J. Matrix Anal. Appl.,
 * 19(4):1097-1110, 1998. doi:10.1137/S0895479896302898,
 * 
 * Authors: Bobby Cheng and Nick Higham, 1996; revised 2015. */
int modchol_ldlt(TMatrixD A, double delta)
{
    int i, j, k, ii, N;

    if (!isHermitian(A)) {
        fprintf(stderr, "Must supply symmetric matrix.");
        return 1;
    }

    N = A.GetNrows();

    TMatrixD D(N,N);

    ldl(A,D); 
    TMatrixD DMC = eye(N);

    // Modified Cholesky perturbations.
    k = 1;
    while (k <= N) {
        if ((k == N) || D(k,k+1) == 0) { // 1-by-1 block
            if (D(k,k) <= delta)
                DMC(k,k) = delta;
            else
                DMC(k,k) = D(k,k);
            k = k+1;
        } else { // 2-by-2 block
            //E = D(k:k+1,k:k+1);
            TMatrixD E = D.GetSub(k,k+1,k,k+1);
            TMatrixD U(2,2);
            TMatrixD T = eig(E,U);
            for (ii = 1; ii <= 2; ii++) {
                if (T(ii,ii) <= delta)
                    T(ii,ii) = delta;
            }
            TMatrixD temp = U*T*U.T();
            //DMC(k:k+1,k:k+1) = (temp + temp.Transpose())/2;  // Ensure symmetric.
            for (i = k; i <= k+1; i++)
                for (j = k; j <= k+1; j++)
                    DMC(i,j) = (temp(i,j) + temp(j,i))/2;  // Ensure symmetric.
            k = k + 2;
        }
    }

    return 0;
}

/** Parametrization of the error matrix in the curvilinear frame.
 *  This frame is tangent to the track at the point of definition,
 *  with Z_T parallel to the track. X_T is in the global xy plane 
 *  and points to the left when looking into the direction of the track,
 *  and Y_T forms a right-handed frame with X_T and Z_T.
 * 
 *  The error along Z_T is therefore zero.
 *  The parameters are <BR>
 *    sigma^2( charge / abs_momentum) <BR>
 *    sigma^2( lambda) <BR>
 *    sigma^2( phi) <BR>
 *    sigma^2( x_transverse)) <BR>
 *    sigma^2( y_transverse)) <BR> <BR>
 *
 *  Please note that lambda and phi are defined in the global frame. Lambda is the helix
 *  dip angle (pi/2 minus theta (polar angle)), while phi is the angle of 
 *  inclination with the global x-axis in the transverse (global xy) plane.
 */

class CurvilinearTrajectoryError {
public:

  /// parameter dimension
  enum { dimension = 5 };
  /// 5 parameter covariance matrix
  typedef math::Error<dimension>::type MathCovarianceMatrix;

// construct
  CurvilinearTrajectoryError() {}

  CurvilinearTrajectoryError(InvalidError) : theCovarianceMatrix(ROOT::Math::SMatrixNoInit()) {theCovarianceMatrix(0,0)=-99999.e10;}

  /** Constructing class from a full covariance matrix. The sequence of the parameters is
   *  the same as the one described above.
   */
  CurvilinearTrajectoryError(const AlgebraicSymMatrix55& aCovarianceMatrix) :
    theCovarianceMatrix(aCovarianceMatrix) {
    int i;
    TVectorD eigValues(5);
    // need to copy the matrix since diagonalization is not
    // available for SMatrix but there must be a better way
    // to copy the matrix ...
    TMatrixDSym diag(5);
    for (i=0; i<5; ++i ) {
        for ( unsigned int j=0; j<5; ++j ) {
            diag(i,j) = theCovarianceMatrix(i,j);
        }
    }
  
    diag.EigenVectors(eigValues);
    for (i = 0; i < 5; i++) {
        if (eigValues(i) < 0) {
            fprintf(stderr, "negative eigenvalue in curvilinear trajectory error!\n");
        }
    }
 }
  template<typename M55>
  CurvilinearTrajectoryError(const M55& aCovarianceMatrix) :
    theCovarianceMatrix(aCovarianceMatrix) {
    int i;
    TVectorD eigValues(5);
    // need to copy the matrix since diagonalization is not
    // available for SMatrix but there must be a better way
    // to copy the matrix ...
    TMatrixDSym diag(5);
    for (i=0; i<5; ++i ) {
        for ( unsigned int j=0; j<5; ++j ) {
            diag(i,j) = theCovarianceMatrix(i,j);
        }
    }
  
    diag.EigenVectors(eigValues);
    for (i = 0; i < 5; i++) {
        if (eigValues(i) < 0) {
            fprintf(stderr, "negative eigenvalue in curvilinear trajectory error!\n");
        }
    }
 }


  bool invalid() const { return theCovarianceMatrix(0,0)<-1.e10;}
  bool valid() const { return !invalid();}

  // not really full check of posdef
  bool posDef() const { 
    return (theCovarianceMatrix(0,0)>=0) && (theCovarianceMatrix(1,1)>=0) && 
      (theCovarianceMatrix(2,2)>=0) && (theCovarianceMatrix(3,3)>=0) && (theCovarianceMatrix(4,4)>=0);
  }


// access

  /** Returning the covariance matrix.
   */
  const AlgebraicSymMatrix55 &matrix() const {
    return theCovarianceMatrix;
  }

 AlgebraicSymMatrix55 &matrix() {
    return theCovarianceMatrix;
  }


  /** Enables the multiplication of the covariance matrix with the scalar "factor".
   */

  void operator *= (double factor) {
    theCovarianceMatrix *= factor;
  }

  void zeroFieldScaling(double factor){
    double root_of_factor = sqrt(factor);
    //scale the 0 indexed covariance by the factor
    for (unsigned int i=1;i!=5;++i)      theCovarianceMatrix(i,0)*=root_of_factor;

    //scale all others by the scared factor
    for (unsigned int i=1;i!=5;++i)  for (unsigned int j=i;j!=5;++j) theCovarianceMatrix(i,j)*=factor;
    //term 0,0 is not scaled at all
  }

  operator MathCovarianceMatrix & () { return theCovarianceMatrix; }
  operator const MathCovarianceMatrix &() const { return theCovarianceMatrix; }

private:
  AlgebraicSymMatrix55 theCovarianceMatrix;
};

#endif
