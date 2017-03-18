using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.model;
using FEA3D.elem;
using FEA3D.fea;

namespace FEA3D.solver
{
    // Solution of the global equation system
    public abstract class Solver 
    {
        // Default solver
        public static solvers solver = solvers.ldu;

        protected static FeModel fem;
        // number of equations
        protected static int neq;
        // length of global stiffness matrix
        public static int lengthOfGSM;
        // elem connectivities - degrees of freedom
        protected int[] indf;
        // number of degrees of freedom for element
        protected int nindf;
        // Indicator of new global matrix
        protected bool newMatrix;

        /*public static enum Solvers {
            ldu {Solver create()
                    {return new SolverLDU();} },
            pcg {Solver create()
                    {return new SolverPCG();} };
            abstract Solver create();
        }*/

        public enum solvers
        {
            ldu, pcg
        };

        private static Solver create()
        {
            switch (solver)
            {
                case solvers.ldu:
                    return new SolverLDU();

                case solvers.pcg:
                    return new SolverPCG();
                    
                default:
                    return new SolverLDU();
            }
        }

        

        public static Solver newSolver(FeModel fem) 
        {
            Solver.fem = fem;
            neq = fem.nEq;
            
            return create();
        }

        // Assemble global stiffnes matrix
        public void assembleGSM() 
        {
            Element elm;
            indf = new int[FE.maxNodesPerElem * fem.nDf];

            for (int iel = 0; iel < fem.nEl; iel++)
            {
                for (int i = 0; i < fem.elems[iel].ind.Length; i++)
                {
                    for (int k = 0; k < fem.nDf; k++)
                        indf[fem.nDf * i + k] =
                            (fem.elems[iel].ind[i] - 1) * fem.nDf + k + 1;
                }
                nindf = fem.elems[iel].ind.Length * fem.nDf;
                elm = fem.elems[iel];
                elm.setElemXy();
                elm.stiffnessMatrix();
                assembleESM();
            }
            // Indicate that new global matrix appeared
            newMatrix = true;
        }

        // Add element stiffness matrix to GSM
        public virtual void assembleESM() {}

        // Solve global equation system
        // x - right-hand side/solution (in/out)
        public virtual int solve(double[] x) 
        {
            return 0;
        }

    }
}
