using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.util;

namespace FEA3D.model
{
    // Load data
    public class FeLoadData 
    {

        protected FeScanner RD;
        public static String loadStepName;
        // Load scale multiplier
        protected double scaleLoad;
        // Relative residual norm tolerance
        public static double residTolerance = 0.01;
        // Maximum number of iterations (elastic-plastic problem)
        public static int maxIterNumber = 100;
        // Degrees of freedom with node forces
        protected List<Dof> nodForces;
        // Element face surface loads
        protected List<ElemFaceLoad> surForces;
        // Temperature increment
        public static double[] dtemp;

        // Increment of force load
        public static double[] dpLoad;
        // Total force load
        public static double[] spLoad;
        // Increment of fictitious thermal loading
        public static double[] dhLoad;
        // Displacement increment
        public static double[] dDispl;
        // Total displacements
        public static double[] sDispl;
        // Right-hand side of global equation system
        public static double[] RHS;

        // Working arrays
        protected static int[] iw = new int[8];
        protected static double[] dw = new double[8];
        protected static double[][] box = new double[2][]{new double[3], new double[3]};

        protected enum vars {
            loadstep, scaleload, residtolerance, maxiternumber,
            nodforce, surforce, boxsurforce, nodtemp, anglesurforce, surforcelist, elemanglesurforce,
            includefile, end, NONE
        };

    }
}
