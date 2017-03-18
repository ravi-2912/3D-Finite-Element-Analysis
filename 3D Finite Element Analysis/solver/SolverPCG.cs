using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.fea;
using FEA3D.model;
using FEA3D.elem;
using FEA3D.util;

namespace FEA3D.solver
{
    
    // Preconditioned conjugate gradient (PCG) solver.
    // Matrix storage: sparse-row format
    public class SolverPCG : Solver 
    {

        // Pointers to matrix rows
        private int[] prow;
        // Column numbers for non-zero values in A
        private int[] coln;
        // Nonzero values of the global stiffness matrix by rows
        private double[] A;
        //mAdcOW.DataStructures.Array<double> A;
        // Working arrays
        private double[] b;
        private double[] r;
        private double[] w;
        private double[] p;
        private double[] md;

        // Constructor for PCG solver.
        public SolverPCG() {

            // Set structure of matrix A
            setSparseRowStructure();

            //A = new double[prow[neq]];
            string path = AppDomain.CurrentDomain.BaseDirectory;
            A = new double[neq];

            b = new double[neq];
            r = new double[neq];
            w = new double[neq];
            p = new double[neq];
            md = new double[neq];
        }

        // Set sparse row structure for storage of nonzero
        //  coefficients of the global stiffness matrix.
        private void setSparseRowStructure() {
            prow = new int[neq+1];
            int lrow = fem.nDf *((fem.nDf == 2) ? FE.maxRow2D : FE.maxRow3D);
            coln = new int[neq*lrow];
            for (int i = 0; i <= fem.nNod; i++)
                prow[i] = i*lrow*fem.nDf;

            // Create nodal sparse-row matrix structure
            for (int i = 0; i < prow[fem.nNod]; i++) coln[i] = -1;
            // Diagonal entry - first in row
            for (int i = 0; i < fem.nNod; i++)  coln[prow[i]] = i;
            for (int iel = 0; iel < fem.nEl; iel++) {
                foreach (int anInd in fem.elems[iel].ind) 
                {
                    if (anInd == 0) continue;
                    int ii = anInd - 1;  // Hyperrow
                    foreach (int anInd1 in fem.elems[iel].ind) 
                    {
                        if (anInd1 == 0) continue;
                        int jj = anInd1 - 1;  // Hypercolumn
                        int k;
                        for (k=prow[ii]; k<prow[ii+1]; k++) 
                        {
                            // If column already exists
                            if (coln[k] == jj) break;
                            if (coln[k] == -1) 
                            {
                                coln[k] = jj;
                                break;
                            }
                        }
                        if (k == prow[ii + 1])
                            UTIL.errorMsg("PCG sparse-row structure: not enough space for node " + ii);
                    }
                }
            }

            // Compress
            int p = 0;
            for (int i = 0; i < fem.nNod; i++) {
                int k = prow[i];
                prow[i] = p;
                for (int j = k; j < prow[i + 1]; j++) {
                    if (coln[j] == -1) break;
                    coln[p++] = coln[j];
                }
            }
            prow[fem.nNod] = p;

            // Transform to degrees of freedom
            int pdof = p*fem.nDf *fem.nDf;
            for (int i = fem.nNod - 1; i >= 0; i--) {
                int deln = (prow[i+1] - prow[i]);
                p -= deln;
                for (int k = fem.nDf; k > 0; k--) {
                    prow[i*fem.nDf+k] = pdof;
                    pdof -= deln*fem.nDf;
                    for (int j = prow[i+1]-prow[i]-1; j>=0; j--) {
                        for (int m = fem.nDf-1; m >= 0; m--) {
                            coln[pdof+j*fem.nDf+m] =
                                    coln[p+j]*fem.nDf + m;
                        }
                    }
                }
            }
            lengthOfGSM = (int) (prow[neq]*1.5);
        }

        // Assemble element matrix to the global stiffness matrix
        public override void assembleESM() {
            for (int j = 0; j < nindf; j++) {
                int jj = indf[j] - 1;
                if (jj >= 0) {
                    for (int i = 0; i < nindf; i++) {
                        int ii = indf[i] - 1;
                        if (ii >= 0) {
                            // Sparse row format (full matrix)
                            int k;
                            for (k=prow[jj]; k<prow[jj+1]; k++)
                                if (coln[k] == ii) break;
                            if (i <= j) A[k] += Element.kmat[i][j];
                            else        A[k] += Element.kmat[j][i];
                        }
                    }
                }
            }
        }

        // Solve equation system by PCG method.
        // x - right-hand side/solution (in/out)
        public override int solve(double[] x) {

            if (newMatrix) {
                displacementBC();
                newMatrix = false;
            }
            return pcg(x);
        }

        // Apply displacement boundary conditions:
        private void displacementBC() 
        {
            IEnumerator<Dof> it = fem.defDs.GetEnumerator();
            //ListIterator it = fem.defDs.listIterator(0);
            Dof d;
            while (it.MoveNext()) {
                d = (Dof) it.Current;
                int j = d.dofNum-1;
                for (int k = prow[j]; k < prow[j+1]; k++)
                    if (coln[k] == j) A[k] = FE.bigValue;
                    else              A[k] = 0.0;
            }
        }

        // PCG solution method.
        // x - right-hand side/solution (in/out)
        public int pcg(double[] x) {

            diagonalPreconditioner();

            // Save x[] in b[] and calculate initial x
            for (int i = 0; i < neq; i++) {
                b[i] = x[i];
                x[i] = x[i]*md[i];
            }

            // r = b - A*x and initinal error
            matrixVectorProduct(x, r);
            for (int i = 0; i < neq; i++) r[i] = b[i] - r[i];
            double gamma0 = 1;
            double gammai = 1;
            int iter;
            for (iter = 0; iter < FE.maxIterPcg; iter++) {
                //  w = (M-1)*r
                for (int i = 0; i < neq; i++)
                    w[i] = md[i]*r[i];
                // gam = (r,w)
                double gammai1 = gammai;
                gammai = 0;
                for (int i = 0; i < neq; i++)
                    gammai += r[i]*w[i];
                if (iter == 0) {
                    gamma0 = gammai;
                    System.Array.ConstrainedCopy(w, 0, p, 0, neq);
                }
                else {
                    double rg = gammai /gammai1;
                    for (int i = 0; i < neq; i++)
                        p[i] = w[i] + rg*p[i];
                }
                // w = A*p
                matrixVectorProduct(p, w);
                double beta = 0;
                for (int i = 0;  i < neq; i++) beta += p[i]*w[i];
                double alpha = gammai/beta;
                // Update x and r, calculate error
                for (int i = 0; i < neq; i++) {
                    x[i] += alpha*p[i];
                    r[i] -= alpha*w[i];
                }
                double err = Math.Sqrt(gammai /gamma0);
                if (err < FE.epsPCG) return (iter + 1);
            }
            return (iter);
        }

        // Diagonal preconditioner md = (Diag(A))-1.
        private void diagonalPreconditioner() {
            for (int j = 0; j < neq; j++) {
                int i;
                for (i = prow[j]; i < prow[j+1]; i++)
                        if (coln[i] == j) break;
                md[j] = 1.0/A[i];
            }
        }

        //  Sparse matrix-vector product y = A*x.
        private void matrixVectorProduct(double[] x, double[] y) {
            if (FE.tunedSolver) {
                if (fem.nDf == 2) {  // tuned for nDf = 2
                    for (int j = 0; j < neq; j++) {
                        double s = 0;
                        for (int i = prow[j]; i < prow[j+1]; i+=2)
                            s += A[i  ]*x[coln[i  ]]
                               + A[i+1]*x[coln[i+1]];
                        y[j] = s;
                    }
                }
                else {    // tuned for nDf = 3
                    for (int j = 0; j < neq; j++) {
                        double s = 0;
                        for (int i = prow[j]; i < prow[j+1]; i+=3)
                            s += A[i  ]*x[coln[i  ]]
                               + A[i+1]*x[coln[i+1]]
                               + A[i+2]*x[coln[i+2]];
                        y[j] = s;
                    }
                }
            }
            else {        // not tuned
                for (int j = 0; j < neq; j++) {
                    double s = 0;
                    for (int i = prow[j]; i < prow[j+1]; i++)
                        s += A[i]*x[coln[i]];
                    y[j] = s;
                }
            }
        }

    }

}
