using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEA3D.elem
{
    // Stresses and equivalent strains at integration point
    public class StressContainer
    {
        // Accumulated stress
        public double [] sStress;

        // StressContainer increment
        public double[] dStress;

        // Accumulated equivalent plastic strain
        public double sEpi;

        // Equivalent plastic strain increment
        public double dEpi;

        public StressContainer(int nDim) 
        {
            sStress = new double[2*nDim];
            dStress = new double[2*nDim];
        }
    }
}
