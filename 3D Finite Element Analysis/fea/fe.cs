using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEA3D.fea
{
    public class FE
    {
        // Main method of application: JFEM/JMGEN/JVIS
        public static int main;
        public const int FEM = 0, MGEN = 1, VIS = 2;

        public static readonly int maxNodesPerElem = 20;

        // Big value for displacements boundary conditions
        public static double bigValue = 1.0e64;
        // Solution tuning
        public static bool tunedSolver = true;

        // Error tolerance for PCG solution of FE equation system
        public static readonly double epsPCG = 1.0e-10;
        // Constants for PCG solution method
        public static int maxRow2D = 21,
                          maxRow3D = 117,
                          maxIterPcg = 10000;

        // Integration scheme for elastic-plastic problems
        public static bool epIntegrationTANGENT = false;
    }
}
