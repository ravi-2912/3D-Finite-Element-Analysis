using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.elem;
using FEA3D.util;
using FEA3D.fea;

namespace FEA3D.material
{
    //  Constitutive relations for elastic-plastic material
    public class ElasticPlasticMaterial : ElasticMaterial
    {

        // Stress at the beginning of increment
        private static double[] sig0 = new double[6];
        // Stress at the end of increment
        private static double[] sig = new double[6];
        // Derivatives of yield function
        private static double[] a = new double[6];
        // Elasticity matrix
        private static double[][] Emat = new double[6][]{new double[6], new double[6], new double[6], 
                                                          new double[6], new double[6], new double[6]};
        // Ea = Emat*a
        private static double[] Ea = new double[6];
        // Plastic strain increment
        private static double[] depsp = new double[6];
        // Equivalent plastic strain
        private static double epi;
        // Shear modulus
        private static double G;
        private static readonly double beta = 0.1;
        private Midpoint midpoint;

        
        public ElasticPlasticMaterial(String stressState) : base(stressState)
        {
            if (stressState.Equals("plstress"))
                UTIL.errorMsg("Elastic-plastic material is not implemented for plane stress");
            if (!FE.epIntegrationTANGENT)
                midpoint = new Midpoint(this);
        }

        // Elastic-plastic stress increment.
        // elm - element,
        // ip - integration point within element
        public override void strainToStress(Element elm, int ip)
        {

            elasticityMatrix(Emat);
            G = getMu();

            // Elastic stress increment dsig due to deps
            
            base.strainToStress(elm, ip);

            for (int i = 0; i < lv; i++) {
                sig0[i] = elm.str[ip].sStress[i];
                sig[i] = sig0[i] + dsig[i];
            }
            // Equivalent plastic strain
            epi = elm.str[ip].sEpi;
            double epi0 = epi;
            double f1 = yieldFunction(sig, epi);

            // Elastic point
            if (f1 < 0.0) return;

            // Elastic-plastic point
            double f0 = yieldFunction(sig0, epi);
            double r;
            if (f0 < 0.0) {
                r = -f0/(f1 - f0);
                for (int i = 0; i < lv; i++)
                    sig[i] = sig0[i] + dsig[i]*r;
                double f = yieldFunction(sig, epi);
                derivYieldFunc(sig, a);
                double c1 = 0.0;
                for (int i = 0; i < lv; i++) c1 += a[i]*dsig[i];
                r = r - f/c1;
            }
            else r = 0.0;

            // Number of subincrements ( = 1 for midpoint method)
            int nsub = (FE.epIntegrationTANGENT) ?
                    (int) (f1/(beta*sY))+1 : 1;

            for (int i = 0; i < lv; i++) {
                sig[i] = sig0[i] + dsig[i]*r;
                dsig[i] = (1.0 - r)*dsig[i]/nsub;
                deps[i] = (1.0 - r)*deps[i]/nsub;
            }
            // Subincrement loop: tangent or midpoint method
            for (int isub = 0; isub < nsub; isub++) {
                if (FE.epIntegrationTANGENT)
                     tangentStressIncrement();
                else midpoint.stressIncrement();
            }
            for (int i = 0; i < lv; i++)
                elm.str[ip].dStress[i] = sig[i] - sig0[i];
            elm.str[ip].dEpi = epi - epi0;
        }

        // Compute elastic-plastic increment by tangent method.
        // Update stresses sig and equivalent plastic strain epi
        private void tangentStressIncrement() {

            double H = slopeH(epi);
            derivYieldFunc(sig, a);
            double dlambda = 0.0;
            for (int i = 0; i < lv; i++) {
                double s = 0.0;
                for (int j = 0; j < lv; j++) s += Emat[i][j]*a[j];
                Ea[i] = s;
                dlambda += a[i]*dsig[i];
            }
            dlambda /= (H + 3*G);
            if (dlambda < 0.0) dlambda = 0.0;
            for (int i = 0; i < lv; i++) {
                sig[i] += dsig[i] - dlambda*Ea[i];
                depsp[i] = dlambda*a[i];
            }
            epi += eqPlastStrain(depsp);

            // Stress correction
            double f = yieldFunction(sig, epi);
            double c1 = 0.0;
            for (int i = 0; i < lv; i++) c1 += a[i]*a[i];
            for (int i = 0; i < lv; i++) sig[i] -= a[i]*f/c1;
        }

        // Yield function.
        // s - stresses,
        // epi - equivalent plastic strain,
        // returns  yield function value
        private double yieldFunction(double[] s, double epi) {
            double sm, seq;
            if (stressState.Equals("threed")) {
                sm = (s[0] + s[1] + s[2])/3;
                seq = Math.Sqrt(3.0*(0.5*((s[0] - sm)*(s[0] - sm)
                        + (s[1] - sm)*(s[1] - sm)
                        + (s[2] - sm)*(s[2] - sm))
                        + s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
            }
            else {
                sm = (s[0] + s[1] + s[3])/3;
                seq = Math.Sqrt(3.0*(0.5*((s[0] - sm)*(s[0] - sm)
                        + (s[1] - sm)*(s[1] - sm)
                        + (s[3] - sm)*(s[3] - sm)) + s[2]*s[2]));
            }
            return seq - yieldRadius(epi);
        }

        // Radius of yield surface Y = sY + k*ep^m
        private double yieldRadius(double ep) {
            if (ep <= 0.0) return sY;
            else return (km*Math.Pow(ep, mm) + sY);
        }

        // Derivatives of yield function.
        // s - stresses (in),
        // a - derivatives of yield function (out)
        private void derivYieldFunc(double[] s, double[] a) {
            double sm, seq;

            if (stressState.Equals("threed")) {
                sm = (s[0] + s[1] + s[2])/3;
                a[0] = s[0] - sm;
                a[1] = s[1] - sm;
                a[2] = s[2] - sm;
                a[3] = 2*s[3];
                a[4] = 2*s[4];
                a[5] = 2*s[5];
                seq = Math.Sqrt(
                        0.5*(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
                        + s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
            }
            else {
                sm = (s[0] + s[1] + s[3])/3;
                a[0] = s[0] - sm;
                a[1] = s[1] - sm;
                a[2] = 2*s[2];
                a[3] = s[3] - sm;
                seq = Math.Sqrt(0.5*(a[0]*a[0] + a[1]*a[1]
                        + a[3]*a[3]) + s[2]*s[2]);
            }

            for (int i = 0; i < lv; i++)
                a[i] = 0.5*Math.Sqrt(3.0)/seq*a[i];
        }

        // Returns slope of deformation curve
        // epi - equivalent plastic strain
        protected double slopeH(double epi) {
            if (km == 0.0) return 0.0;
            else if (mm == 1.0) return km;
            else return (epi == 0.0) ? 0.0 : km*mm*Math.Pow(epi, mm - 1.0);
        }

        // Returns equivalent plastic strain
        // dp - pastic strains (in)
        private double eqPlastStrain(double[] dp) {
            if (stressState.Equals("threed"))
                return Math.Sqrt(
                   (2*(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2])
                   + dp[3]*dp[3] + dp[4]*dp[4] + dp[5]*dp[5])/3.0);
            else {
                if (stressState.Equals("plstress"))
                    dp[3] = -(dp[0] + dp[1]);
                return Math.Sqrt((2*(dp[0]*dp[0] + dp[1]*dp[1]
                        + dp[3]*dp[3]) + dp[2]*dp[2])/3.0);
            }
        }

        // Midpoint method for integration of constitutive
        // relations for the von Mises hardening material
        class Midpoint 
        {
            ElasticPlasticMaterial epm;
            // Strain deviator
            double[] ed = new double[6];
            // Stress deviator
            double[] sd = new double[6];
            // Trial stress deviator
            double[] sdtr = new double[6];
            // Yield function derivatives
            double[] a0 = new double[6];

            double[] sal = new double[6];
            double[] salbar = new double[6];
            double[] b = new double[6];
            static double alpha = 0.5;

            public Midpoint(ElasticPlasticMaterial e)
            {
                epm = e;
            }

            
            // Elastic-plastic stress increment by midpoint method.
            // Update stresses sig and
            //    equivalent plastic strain epi
            public void stressIncrement() 
            {
                double SQ32 = Math.Sqrt(1.5);
                double SQ23 = Math.Sqrt(2.0/3.0);
                double tolerance = 1.0e-5;
                // Transform strains to tensor components
                if (lv == 4) deps[2] *= 0.5;
                else for (int i = 3; i < 6; i++) deps[i] *= 0.5;
                double depsm = deviator(deps, ed);
                double sigm = deviator(sig, sd);
                double sigeq0 = SQ32*norm(sd);
                for (int i = 0; i < lv; i++) {
                    sdtr[i] = sd[i] + 2*G*ed[i];
                    a0[i] = 1.5*sd[i]/sigeq0;
                }
                double lambda = SQ23*norm(ed);
                // Find lambda by Newton-Raphson iteration
                double epi1, sigeq;
                for (; ;) {
                    for (int i = 0; i < lv; i++)
                        sal[i] = sdtr[i] -
                                2 * G * lambda * (1 - alpha) * a0[i];
                    double salmod = norm(sal);
                    for (int i = 0; i < lv; i++) {
                        salbar[i] = sal[i]/salmod;
                        b[i] = (1 - alpha) * a0[i] +
                                alpha * SQ32 * salbar[i];
                    }
                    double bmod = norm(b);
                    epi1 = epi + lambda*SQ23*bmod;
                    sigeq = epm.yieldRadius(epi1);
                    double phi = SQ32*salmod -
                            3.0 * G * alpha * lambda - sigeq;
                    if (Math.Abs(phi) < tolerance * epm.sY) break;
                    double phiPrime = -2.0*SQ32*G
                            * (1 - alpha) * dyadicProduct(a0, salbar)
                            - 3.0 * G * alpha - epm.slopeH(epi1) * SQ23 * bmod;
                    double lambda1 = lambda - phi/phiPrime;
                    lambda = (lambda1 <= 0.0) ? 0.5*lambda:lambda1;
                }
                epi = epi1;
                for (int i = 0; i < lv; i++)
                    sd[i] = sal[i] / (1 + 3 * G * alpha * lambda / sigeq);
                double dsigm = depsm * epm.e / (1.0 - 2.0 * epm.nu);
                stressFromDeviator(sd, sigm + dsigm, sig);
            }

            // Compute deviator.
            // s - stress,
            // d - deviator (out),
            // returns  mean value
            double deviator(double[] s, double[] d) {
                double sm;
                if (lv == 4) {
                    sm = (s[0] + s[1] + s[3])/3;
                    d[0] = s[0] - sm;
                    d[1] = s[1] - sm;
                    d[3] = s[3] - sm;
                    d[2] = s[2];
                }
                else {
                    sm = (s[0] + s[1] + s[2])/3;
                    d[0] = s[0] - sm;
                    d[1] = s[1] - sm;
                    d[2] = s[2] - sm;
                    d[3] = s[3];
                    d[4] = s[4];
                    d[5] = s[5];
                }
                return (sm);
            }

            // Compute stress s = d + sm.
            // d - deviator,
            // sm - mean stress,
            // s - stress (out)
            void stressFromDeviator(double[] d, double sm,
                                    double[] s) {
                if (lv == 4) {
                    s[0] = d[0] + sm;
                    s[1] = d[1] + sm;
                    s[3] = d[3] + sm;
                    s[2] = d[2];
                }
                else {
                    s[0] = d[0] + sm;
                    s[1] = d[1] + sm;
                    s[2] = d[2] + sm;
                    s[3] = d[3];
                    s[4] = d[4];
                    s[5] = d[5];
                }
            }

            // Returns norm = sqrt(Sij*Sij)
            double norm(double[] s) {
                if (lv == 4)
                    return (Math.Sqrt(s[0]*s[0] + s[1]*s[1] +
                            s[3]*s[3] + 2*s[2]*s[2]));
                else
                    return (Math.Sqrt(s[0]*s[0] + s[1]*s[1] +
                            s[2]*s[2] + 2*(s[3]*s[3] + s[4]*s[4] +
                            s[5]*s[5])));
            }

            // Returns dyadic product = aij*bij
            double dyadicProduct(double[] a, double[] b) {
                if (lv == 4)
                    return (a[0]*b[0] + a[1]*b[1] + a[3]*b[3] +
                            2*a[2]*b[2]);
                else
                    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] +
                            2*(a[3]*b[3] + a[4]*b[4] + a[5]*b[5]));
            }

        }

    }
}
