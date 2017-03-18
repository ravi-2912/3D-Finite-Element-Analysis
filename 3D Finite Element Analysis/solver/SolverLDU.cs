using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.elem;
using FEA3D.model;
using FEA3D.fea;


namespace FEA3D.solver
{
    // Profile LDU symmetric solver.
    // Upper symmetric part of the global stiffness matrix is
    // stored by columns of variable height (profile storage).
    public class SolverLDU : Solver 
    {

        // Pointers to matrix columns
        private int[] pcol;
        // Global stiffness matrix
        private double[] A;
        //mAdcOW.DataStructures.Array<double> A;
        
        // Constructor for LDU symmetric solver.
        public SolverLDU() 
        {
            // Create profile of the global stiffness matrix
            setProfile();
            string path = AppDomain.CurrentDomain.BaseDirectory;
            //A = new mAdcOW.DataStructures.Array<double>(pcol[neq], path); 
            A = new double[pcol[neq]];         
            
        }

        // Calculate profile of GSM: set column pointers pcol[]
        private void setProfile() 
        {
            pcol = new int[neq + 1];
            for (int i = 0; i < neq+1; i++) pcol[i] = 0;

            for (int iel = 0; iel < fem.nEl; iel++) {
                foreach (int jind in fem.elems[iel].ind) {
                    // Calculate profile: nodal hypercolumn length
                    // is in the first entry for a hypercolumn
                    if (jind > 0) {
                        int ncol = (jind - 1) * fem.nDf + 1;
                        int icolh = pcol[ncol];
                        foreach (int kind in fem.elems[iel].ind) {
                            if (kind > 0) icolh =
                                    Math.Max(icolh, jind - kind);
                        }
                        pcol[ncol] = icolh;
                    }
                }
            }
            //  Transform hypercolumn lengths to column addresses
            int ic = 0;
            for (int i = 0; i < fem.nNod; i++) {
                int icolh = pcol[i*fem.nDf + 1]*fem.nDf;
                for (int j = 0; j < fem.nDf; j++) {
                    ic++;
                    icolh++;
                    pcol[ic] = pcol[ic - 1] + icolh;
                }
            }
            lengthOfGSM = pcol[neq];
        }

        // Assemble element matrix to the global stiffness matrix
        public override void assembleESM() 
        {
            // Assemble all contributions to the upper part of GSM
            for (int j = 0; j < nindf; j++)
            {
                int jj = indf[j] - 1;
                if (jj >= 0)
                {
                    for (int i = 0; i < nindf; i++)
                    {
                        int ii = indf[i] - 1;
                        if (ii >= 0)
                        {
                            // Profile format (columns top->bottom)
                            if (ii <= jj)
                            {
                                int iad = pcol[jj + 1] - jj + ii - 1;
                                if (i <= j)
                                    A[iad] += Element.kmat[i][j];
                                else
                                    A[iad] += Element.kmat[j][i];
                            }
                        }
                    }
                }
            }
        }

        // Solve equation system by direct LDU method.
        // x - right-hand side/solution (in/out)
        public override int solve(double[] x) 
        {
            if (newMatrix)
                {
                    displacementBC();
                    if (FE.tunedSolver)
                        lduDecompositionTuned();
                    else
                        lduDecomposition();
                    newMatrix = false;
                }
            lduFrwdBksb(x);
            
            return 0;
        }

        // Apply displacement boundary conditions
        private void displacementBC() 
        {
            IEnumerator<Dof> it = fem.defDs.GetEnumerator();//.listIterator(0);
            Dof d;
            while (it.MoveNext()) 
            {
                d = it.Current;
                A[pcol[d.dofNum] - 1] = FE.bigValue;
                
            }
        }

        // LDU decomposition for symmetric matrix
        private int lduDecomposition() 
        {
            // Working array
            double[] w = new double[neq];

            // UtDU decomposition
            for (int j = 1; j < neq; j++)
            {
                int jfirst = j - (pcol[j + 1] - pcol[j]) + 1;
                int jj = pcol[j + 1] - j - 1;
                for (int i = jfirst; i < j; i++)
                    w[i] = A[i + jj] / A[pcol[i + 1] - 1];
                for (int i = j; i < neq; i++)
                {
                    int ifirst = i - (pcol[i + 1] - pcol[i]) + 1;
                    int ii = pcol[i + 1] - i - 1;
                    double s = 0.0;
                    for (int m = Math.Max(jfirst, ifirst);
                             m < j; m++) s += A[m + ii] * w[m];
                    A[j + ii] -= s;
                }
            }
            for (int j = 0; j < neq; j++)
            {
                int jfirst = j - (pcol[j + 1] - pcol[j]) + 1;
                int jj = pcol[j + 1] - j - 1;
                for (int i = jfirst; i < j; i++)
                    A[i + jj] /= A[pcol[i + 1] - 1];
            }
            return 0;
        }

        // Forward reduction and backsubstitution
        private void lduFrwdBksb(double[] x)
        {
            //  b =(U)-T*b
            for (int j = 1; j < neq; j++)
            {
                int jfirst = j - (pcol[j + 1] - pcol[j]) + 1;
                int jj = pcol[j + 1] - j - 1;
                for (int i = jfirst; i < j; i++)
                    x[j] -= A[i + jj] * x[i];
            }
            //  b = (U)-1*(D)-1*b
            for (int i = neq - 1; i >= 0; i--)
            {
                x[i] /= A[pcol[i + 1] - 1];
                for (int j = i + 1; j < neq; j++)
                {
                    int jfirst = j - (pcol[j + 1] - pcol[j]) + 1;
                    if (i >= jfirst)
                    {
                        int jj = pcol[j + 1] - j - 1;
                        x[i] -= A[i + jj] * x[j];
                    }
                }
            }
        }

        // Tuned LDU decomposition for symmetric matrix
        // (block-block tuning, block size = 2 - 2D; 3 - 3D)
        private int lduDecompositionTuned()
        {
            

            double s0, s1, s2, t0, t1, t2, u0, u1, u2;
            int ndf = fem.nDf;
            double[] w = new double[ndf*neq];
            double[] z = new double[neq];

            int maxcol = 0;
            for (int i = 0; i < neq; i++)
                maxcol = Math.Max(maxcol, pcol[i+1] - pcol[i]);

            // UtDU decomposition
            for (int i = 0; i < neq; i++)
                z[i] = 1.0/A[pcol[i+1]-1];

            for (int jc = 1; jc <= neq; jc += ndf) {
                int jfirst = jc - (pcol[jc+1-1] - pcol[jc-1]) + 1;
                int n1 = pcol[jc+ndf-1] - pcol[jc+ndf-2];
                int jw = n1 - (jc + ndf - 1);
                for (int j = jc; j < jc + ndf; j++) {
                    int jj = pcol[j] - j;
                    for (int i = jfirst; i < j; i++)
                        w[i+jw+n1*(j-jc)-1] = A[i+jj-1]*z[i-1];
                    for (int i = j; i < jc + ndf; i++) {
                        int ifirst = i - (pcol[i] - pcol[i-1]) + 1;
                        int ii = pcol[i] - i;
                        s0 = 0.0;
                        for (int m = Math.Max(jfirst, ifirst);
                                 m < j; m++)
                            s0 += A[m+ii-1]*w[m+jw+n1*(j-jc)-1];
                        A[j+ii-1] -= s0;
                    }
                    z[j-1] = 1.0/A[pcol[j]-1];
                    w[(j-jc+1)*n1-1] = z[j-1];
                }
                if (ndf == 2) {
                    for (int i = jc+ndf;
                         i<Math.Min(neq, jc+ndf+maxcol); i+=ndf) {
                        int ifirst = i - (pcol[i] - pcol[i-1]) + 1;
                        int ii = pcol[i]-i;
                        int ii1 = pcol[i+1]-(i+1);
                        s0 = s1 = t0 = t1 = 0;
                        for (int m = Math.Max(jfirst, ifirst)-1;
                                                  m < jc-1; m++) {
                            s0 += A[m+ii]*w[m+jw];
                            s1 += A[m+ii]*w[m+jw+n1];
                            t0 += A[m+ii1]*w[m+jw];
                            t1 += A[m+ii1]*w[m+jw+n1];
                        }
                        if (jc >= ifirst) {
                            A[jc+ii-1] -= s0;
                            A[jc+ii1-1] -= t0;
                            s1 += A[jc+ii-1]*w[jc+jw+n1-1];
                            A[jc+1+ii-1] -= s1;
                            t1 += A[jc+ii1-1]*w[jc+jw+n1-1];
                            A[jc+1+ii1-1] -= t1;
                        }
                    }
                }
                else {   // tuned for nDf = 3
                    for (int i = jc+ndf;
                         i<Math.Min(neq, jc+ndf+maxcol); i+=ndf) {
                        int ifirst = i - (pcol[i] - pcol[i-1]) + 1;
                        int ii = pcol[i]-i;
                        int ii1 = pcol[i+1]-(i+1);
                        s0 = s1 = t0 = t1 = 0;
                        int ii2 = pcol[i+2]-(i+2);
                        s2 = t2 = u0 = u1 = u2 = 0;
                        for (int m = Math.Max(jfirst, ifirst)-1;
                                                  m < jc-1; m++) {
                            s0 += A[m+ii]*w[m+jw];
                            s1 += A[m+ii]*w[m+jw+n1];
                            s2 += A[m+ii]*w[m+jw+2*n1];
                            t0 += A[m+ii1]*w[m+jw];
                            t1 += A[m+ii1]*w[m+jw+n1];
                            t2 += A[m+ii1]*w[m+jw+2*n1];
                            u0 += A[m+ii2]*w[m+jw];
                            u1 += A[m+ii2]*w[m+jw+n1];
                            u2 += A[m+ii2]*w[m+jw+2*n1];
                        }
                        if (jc >= ifirst) {
                            A[jc+ii-1] -= s0;
                            A[jc+ii1-1] -= t0;
                            s1 += A[jc+ii-1]*w[jc+jw+n1-1];
                            A[jc+1+ii-1] -= s1;
                            t1 += A[jc+ii1-1]*w[jc+jw+n1-1];
                            A[jc+1+ii1-1] -= t1;
                            A[jc+ii2-1] -= u0;
                            u1 += A[jc+ii2-1]*w[jc+jw+n1-1];
                            A[jc+1+ii2-1] -= u1;
                            s2 += A[jc+ii-1]*w[jc+jw+2*n1-1]
                                +A[jc+1+ii-1]*w[jc+jw+1+2*n1-1];
                            A[jc+2+ii-1] -= s2;
                            t2 += A[jc+ii1-1]*w[jc+jw+2*n1-1]
                                +A[jc+1+ii1-1]*w[jc+jw+1+2*n1-1];
                            A[jc+2+ii1-1] -= t2;
                            u2 += A[jc+ii2-1]*w[jc+jw+2*n1-1]
                                +A[jc+1+ii2-1]*w[jc+jw+1+2*n1-1];
                            A[jc+2+ii2-1] -= u2;
                        }
                    }
                }
            }
            for (int j = 1; j <= neq; j++) {
                int jfirst = j - (pcol[j+1-1] - pcol[j-1]) + 1;
                int jj = pcol[j+1-1] - j;
                for (int i = jfirst; i < j; i++)
                    A[i+jj-1] *= z[i-1];
            }
            

            return 0;

        }

    }
}
