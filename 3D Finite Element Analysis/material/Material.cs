using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.elem;

namespace FEA3D.material
{
    public class Material
    {

        // StressContainer state (plstrain/plstress/axisym/threed)
        protected String stressState;
        // Elasticity modulus
        public double e;
        // Poisson's ratio
        public double nu;
        // Thermal expansion
        protected double alpha;
        // Yield stress
        protected double sY;
        // Hardening coefficient
        protected double km;
        // Hardening power
        protected double mm;

        public static Material newMaterial(String matPhysLaw, String stressState)
        {
            if (matPhysLaw.Equals("elastic"))
                return new ElasticMaterial(stressState);
            else return new ElasticPlasticMaterial(stressState);
        }

        // Given strain increment at integration point ip
        // element elm, compute stress dsig increment
        public virtual void strainToStress(Element elm, int ip) { }

        // Set elastic properties
        public void setElasticProp(double e, double nu, double alpha)
        {
            this.e = e;
            this.nu = nu;
            this.alpha = alpha;

        }

        // Set plastic properties
        public void setPlasticProp(double sY, double km, double mm)
        {
            this.sY = sY;
            this.km = km;
            this.mm = mm;
        }

        // Returns Lame constant lambda
        public double getLambda()
        {
            return (stressState.Equals("plstress")) ? e * nu / ((1.0 + nu) * (1.0 - nu)) : e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        }

        // Returns shear modulus
        public double getMu()
        { return 0.5 * e / (1.0 + nu); }

        // Returns Poisson's ratio
        public double getNu() 
        { return nu; }

        // Returns thermal expansion coefficient
        public double getAlpha() 
        { return alpha; }

        // Returns Youngs Modulus
        public double getElasticModulus()
        { return e; }

        // Compute elasticity matrix emat
        public virtual void elasticityMatrix(double[][] emat) { }

    }
}
