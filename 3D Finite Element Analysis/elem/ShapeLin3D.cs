using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.util;

namespace FEA3D.elem
{
    // Quadratic 3D shape functions and their derivatives

    public class ShapeLin3D
    {
        // Degeneration check.
        // The only degeneration is: 0=3 & 4=7
        public static int degeneration(int[] ind)
        {
            // Per face two degeneration cases giving total of 12
            // Only one is considered here.

            // Element should be linear
            if ((ind[0] == ind[3] && ind[4] == ind[7]))
                return 1;
            else return 0;
        }

        // Shape functions.
        // xi, et, ze - local coordinates;
        // ind - element connectivities;
        // n - shape functions (out)
        public static void shape(double xi, double et, double ze, int[] ind, double[] n)
        {
            double s0 = 1 + xi;
            double t0 = 1 + et;
            double d0 = 1 + ze;
            double s1 = 1 - xi;
            double t1 = 1 - et;
            double d1 = 1 - ze;


            // Vertex nodes
            n[0] = 0.125 * s1 * t1 * d1;
            n[1] = 0.125 * s0 * t1 * d1;
            n[2] = 0.125 * s0 * t0 * d1;
            n[3] = 0.125 * s1 * t0 * d1;

            n[4] = 0.125 * s1 * t1 * d0;
            n[5] = 0.125 * s0 * t1 * d0;
            n[6] = 0.125 * s0 * t0 * d0;
            n[7] = 0.125 * s1 * t0 * d0;

            // Modification of functions due to degeneration
            if (degeneration(ind) == 1)
            {
                double i1 = n[0] + n[3];
                double i2 = n[4] + n[7];
                n[0] = n[3] = i1;
                n[4] = n[7] = i2;

            }
        }

        // Derivatives of shape functions
        //    with respect to global coordinates xy.
        // xi, et, ze - local coordinates;
        // ind - element connectivities;
        // xy - nodal coordinates;
        // dnxy - derivatives of shape functions (out);
        // returns  determinant of the Jacobian matrrix
        public static double deriv(double xi, double et, double ze, int[] ind, double[][] xy, double[][] dnxy)
        {
            // Derivatives with respect to local coordinates d
            double[][] d = new double[8][];
            for (int i = 0; i < 8; i++)
                d[i] = new double[3];

            double s0 = 1 + xi;
            double t0 = 1 + et;
            double d0 = 1 + ze;
            double s1 = 1 - xi;
            double t1 = 1 - et;
            double d1 = 1 - ze;

            // Vertex nodes
            d[0][0] = -0.125 * t1 * d1;
            d[0][1] = -0.125 * s1 * d1;
            d[0][2] = -0.125 * s1 * t1;

            d[1][0] =  0.125 * t1 * d1;
            d[1][1] = -0.125 * s0 * d1;
            d[1][2] = -0.125 * s0 * t1;

            d[2][0] = 0.125 * t0 * d1;
            d[2][1] = 0.125 * s0 * d1;
            d[2][2] = -0.125 * s0 * t0;

            d[3][0] = -0.125 * t0 * d1;
            d[3][1] = 0.125 * s1 * d1;
            d[3][2] = -0.125 * s1 * t0;

            d[4][0] = -0.125 * t1 * d0;
            d[4][1] = -0.125 * s1 * d0;
            d[4][2] = 0.125 * s1 * t1;

            d[5][0] = 0.125 * t1 * d0;
            d[5][1] = -0.125 * s0 * d0;
            d[5][2] = 0.125 * s0 * t1;

            d[6][0] = 0.125 * t0 * d0;
            d[6][1] = 0.125 * s0 * d0;
            d[6][2] = 0.125 * s0 * t0;

            d[7][0] = -0.125 * t0 * d0;
            d[7][1] = 0.125 * s1 * d0;
            d[7][2] = 0.125 * s1 * t0;

            // Modification of derivatives due to degeneration
            if (degeneration(ind) == 1)
            {
                for(int i = 0; i < 3; i++)
                {
                    double i1 = d[0][i] + d[3][i];
                    d[0][i] = d[3][i] = i1;

                    double i2 = d[4][i] + d[7][i];
                    d[4][i] = d[7][i] = i2;
                }
            }
            // Jacobian matrix ja
            double[][] ja = new double[3][] { new double[3], new double[3], new double[3] };
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    ja[j][i] = 0.0;
                    for (int k = 0; k < 8; k++)
                        ja[j][i] += d[k][i] * xy[k][j];
                }
            // Determinant of Jacobian matrix det
            double det =
                ja[0][0] * (ja[1][1] * ja[2][2] - ja[2][1] * ja[1][2])
              - ja[1][0] * (ja[0][1] * ja[2][2] - ja[0][2] * ja[2][1])
              + ja[2][0] * (ja[0][1] * ja[1][2] - ja[0][2] * ja[1][1]);
            if (det <= 0)
                UTIL.errorMsg("Negative/zero Jacobian determinant for 8N element " + (float)det);
            // Jacobian inverse ja1
            double[][] ja1 = new double[3][] { new double[3], new double[3], new double[3] };
            double v = 1.0 / det;
            ja1[0][0] = (ja[1][1] * ja[2][2] - ja[2][1] * ja[1][2]) * v;
            ja1[1][0] = (ja[0][2] * ja[2][1] - ja[0][1] * ja[2][2]) * v;
            ja1[2][0] = (ja[0][1] * ja[1][2] - ja[0][2] * ja[1][1]) * v;
            ja1[0][1] = (ja[1][2] * ja[2][0] - ja[1][0] * ja[2][2]) * v;
            ja1[1][1] = (ja[0][0] * ja[2][2] - ja[2][0] * ja[0][2]) * v;
            ja1[2][1] = (ja[0][2] * ja[1][0] - ja[0][0] * ja[1][2]) * v;
            ja1[0][2] = (ja[2][1] * ja[1][0] - ja[2][0] * ja[1][1]) * v;
            ja1[1][2] = (ja[0][1] * ja[2][0] - ja[0][0] * ja[2][1]) * v;
            ja1[2][2] = (ja[0][0] * ja[1][1] - ja[1][0] * ja[0][1]) * v;

            for (int k = 0; k < 8; k++)
                for (int i = 0; i < 3; i++)
                {
                    dnxy[k][i] = 0.0;
                    for (int j = 0; j < 3; j++)
                        dnxy[k][i] += ja1[i][j] * d[k][j];
                }
            return det;
        }

        // Two-dimensional shape functions and derivatives
        //    for a face of 3d 8-20n element.
        // xi, et - local coordinates;
        // ind - element connectivities;
        // an - shape functions (out);
        // dn  derivatives of shape functions
        //      with respect to xi and et (out)
        public static void shapeDerivFace(double xi, double et, int[] ind, double[] an, double[][] dn)
        {
            double x1 = 1 - xi;
            double x2 = 1 + xi;
            double e1 = 1 - et;
            double e2 = 1 + et;

            // Shape functions for corner nodes
            an[0] = 0.25 * (1.0 - xi) * (1.0 - et);
            an[1] = 0.25 * (1.0 + xi) * (1.0 - et);
            an[2] = 0.25 * (1.0 + xi) * (1.0 + et);
            an[3] = 0.25 * (1.0 - xi) * (1.0 + et);

            // Derivatives for corner nodes
            dn[0][0] = -0.25 * (1.0 - et);
            dn[0][1] = -0.25 * (1.0 - xi);
            dn[1][0] = 0.25 * (1.0 - et);
            dn[1][1] = -0.25 * (1.0 + xi);
            dn[2][0] = 0.25 * (1.0 + et);
            dn[2][1] = 0.25 * (1.0 + xi);
            dn[3][0] = -0.25 * (1.0 + et);
            dn[3][1] = 0.25 * (1.0 - xi);

            // Modification due to degeneration
            int deg = 0;
            for (int i = 0; i < 3; i++)
                if (ind[i] == ind[i + 1])
                    deg = i + 1;
            
            if (deg > 0)
            {
                double sum = an[deg - 1] + an[deg];
                an[deg - 1] = an[deg] = sum;

                double sum0 = dn[deg - 1][0] + dn[deg][0];
                double sum1 = dn[deg - 1][1] + dn[deg][1];
                dn[deg - 1][0] = dn[deg][0] = sum0;
                dn[deg - 1][1] = dn[deg][1] = sum1;
            }

            
        }
    }
}
