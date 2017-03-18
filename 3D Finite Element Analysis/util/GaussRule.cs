using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEA3D.util
{
    public class GaussRule
    {
        // Abscissas of the Gauss rule
        public double[] xii, eti, zei;
        
        // Integration weights
        public double[] wi;
        
        // Total namber of integration poins
        public int nIntPoints;

        // Abscissas and weights for 1, 2 and 3-point rules
        private static readonly double[][] X = new double[][]{ new double[]{0.0}, 
                                                               new double[]{-1.0/Math.Sqrt(3.0), 1.0/Math.Sqrt(3.0)}, 
                                                               new double[]{-Math.Sqrt(0.6), 0.0, Math.Sqrt(0.6)} 
                                                           };
        private static readonly double[][] W = new double[][]{ new double[]{2.0}, 
                                                               new double[]{1.0, 1.0}, 
                                                               new double[]{5.0/9.0, 8.0/9.0, 5.0/9.0}
                                                           };
        
        // Abscissas and weights for 14-point rule (3D)
        private static readonly double a = 0.7587869106393281, b = 0.7958224257542215;
        private static readonly double[] X14 = new double[] { -a, a, -a, -a, a, -a, a, a, -b, b, 0.0, 0.0, 0.0, 0.0 };
        private static readonly double[] Y14 = new double[] { -a, -a, a, -a, a, a, -a, a, 0.0, 0.0, -b, b, 0.0, 0.0 };
        private static readonly double[] Z14 = new double[] { -a, -a, -a, a, -a, a, a, a, 0.0, 0.0, 0.0, 0.0, -b, b };
        private static readonly double Wa = 0.3351800554016621, Wb = 0.8864265927977839;

        // Construct Gauss integration rule.
        // nGauss - number of Gauss points in each direction
        // (excluding 14-point rule),
        // nDim - number of dimensions
        public GaussRule(int nGauss, int nDim) 
        {
            if (!((nGauss >= 1 && nGauss <= 3) || nGauss == 14))
                UTIL.errorMsg("nGauss has forbidden value: " + nGauss);
            if (!(nDim >= 1 && nDim <= 3))
                UTIL.errorMsg("GaussRule: nDim has forbidden value: "+ nDim);

            if (nGauss == 14) nIntPoints = 14;
            else
            {
                nIntPoints = 1;
                for (int i = 0; i < nDim; i++)
                    nIntPoints *= nGauss;
            }

            xii = new double[nIntPoints];
            wi = new double[nIntPoints];
            if (nDim > 1) eti = new double[nIntPoints];
            if (nDim > 2) zei = new double[nIntPoints];

            if (nGauss == 14)
            {
                for (int i = 0; i < nGauss; i++)
                {
                    xii[i] = X14[i];
                    eti[i] = Y14[i];
                    zei[i] = Z14[i];
                    wi[i] = (i < 8) ? Wa : Wb;
                }
            }
            else
            {
                int ip = 0;
                int n = nGauss - 1;
                switch (nDim)
                {
                    case 1:
                        for (int i = 0; i < nGauss; i++)
                        {
                            xii[ip] = X[n][i];
                            wi[ip++] = W[n][i];
                        }
                        break;

                    case 2:
                        for (int i = 0; i < nGauss; i++)
                        {
                            for (int j = 0; j < nGauss; j++)
                            {
                                xii[ip] = X[n][i];
                                eti[ip] = X[n][j];
                                wi[ip++] = W[n][i] * W[n][j];
                            }
                        }
                        break;

                    case 3:
                        for (int i = 0; i < nGauss; i++)
                        {
                            for (int j = 0; j < nGauss; j++)
                            {
                                for (int k = 0; k < nGauss; k++)
                                {
                                    xii[ip] = X[n][i];
                                    eti[ip] = X[n][j];
                                    zei[ip] = X[n][k];
                                    wi[ip++] = W[n][i] * W[n][j] * W[n][k];
                                }
                            }
                        }
                        break;
                }
            }
        }
    }
}
