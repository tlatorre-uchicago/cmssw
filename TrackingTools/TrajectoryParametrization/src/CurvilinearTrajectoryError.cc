#include "DataFormats/GeometrySurface/interface/TrivialError.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "DataFormats/Math/interface/Error.h"
#include "TVectorD.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TDecompBK.h"

/* Compute the V = LDL^T factorization of the Matrix V.
 *
 * Taken from https://github.com/nojhan/cholesky/blob/master/cholesky.h. */
int ldl(const TMatrixD &V, TMatrixD &D)
{
    int i, j, k;
    int N = V.GetNrows();
    // example of an invertible matrix whose decomposition is undefined
    assert(V(0,0) != 0); 
    TMatrixD L(N,N);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            L(i,j) = 0;
            D(i,j) = 0;
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
int modchol_ldlt(TMatrixDSym A, TMatrixD &DMC, double delta)
{
    int i, j, k, ii, N;

    if (!isHermitian(A)) {
        fprintf(stderr, "Must supply symmetric matrix.");
        return 1;
    }

    N = A.GetNrows();

    TMatrixD D(N,N);

    fprintf(stderr, "blah = \n");
    TDecompBK D0 = TDecompBK(A); 
    D = TMatrixD(D0.GetU()).T();
    fprintf(stderr, "D0 = \n");
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fprintf(stderr, "%.2e ", D(i,j));
        }
        fprintf(stderr, "\n");
    }
    DMC = eye(N);

    // Modified Cholesky perturbations.
    k = 0;
    while (k < N) {
        if ((k == N-1) || D(k,k+1) == 0) { // 1-by-1 block
            fprintf(stderr, "1-by-1 block\n");
            if (D(k,k) <= delta)
                DMC(k,k) = delta;
            else
                DMC(k,k) = D(k,k);
            k = k+1;
        } else { // 2-by-2 block
            fprintf(stderr, "2-by-2 block\n");
            //E = D(k:k+1,k:k+1);
            TMatrixD E = D.GetSub(k,k+1,k,k+1);
            TMatrixD U(2,2);
            TMatrixD T = eig(E,U);
            for (ii = 0; ii < 2; ii++) {
                if (T(ii,ii) <= delta)
                    T(ii,ii) = delta;
            }
            TMatrixD temp = U*T*U.T();
            //DMC(k:k+1,k:k+1) = (temp + temp.Transpose())/2;  // Ensure symmetric.
            for (i = k; i <= k+1; i++)
                for (j = k; j <= k+1; j++)
                    DMC(i,j) = (temp(i-k,j-k) + temp(j-k,i-k))/2;  // Ensure symmetric.
            k = k + 2;
        }
    }

    return 0;
}

CurvilinearTrajectoryError::CurvilinearTrajectoryError(const AlgebraicSymMatrix55& aCovarianceMatrix) :
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
  CurvilinearTrajectoryError::CurvilinearTrajectoryError(const M55& aCovarianceMatrix) :
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
