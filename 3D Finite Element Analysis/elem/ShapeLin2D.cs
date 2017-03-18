using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.util;

namespace FEA3D.elem
{
    // Quadratic 2D shape functions and their derivatives
    public class ShapeLin2D
    {
        // Degeneration check.
        // If element is triangular then the method return
        // a  number (starting from 1) of the degenerated node.
        // 1 is the node 0, 2 is node 1.
        // ind - connectivity numbers
        static int degeneration(int[] ind)
        {
            int deg = 0;
            for (int i = 0; i < 3; i++)
            {
                if (ind[i] == ind[i + 1])
                {
                    deg = i + 1;//(i + 5) % 8///////;
                    //return deg;
                }
            }
            return deg;
        }

        // Shape functions.
        // xi, et - local coordinates;
        // ind[4] - element connectivities;
        // an[??] - shape functions (out)
        public static void shape(double xi, double et, int[] ind, double[] an)
        {
            // Shape functions of corner nodes

            an[0] = 0.25 * (1.0 - xi) * (1.0 - et);
            an[1] = 0.25 * (1.0 + xi) * (1.0 - et);
            an[2] = 0.25 * (1.0 + xi) * (1.0 + et);
            an[3] = 0.25 * (1.0 - xi) * (1.0 + et);

            // Modification of functions due to degeneration
            int deg = degeneration(ind);
            if (deg > 0)
            {
                double sum = an[deg - 1] + an[deg];
                an[deg - 1] = an[deg] = sum;
            }
        }

        // Derivatives of shape functions
        // with respect to global coordinates x and y.
        // xi, et - local coordinates;
        // ind[4] - element connectivities;
        // xy[4][2] - nodal coordinates;
        // dnxy[4][2] - derivatives of shape functions (out);
        // returns  determinant of the Jacobian matrrix
        public static double deriv(double xi, double et, int[] ind, double[][] xy, double[][] dnxy)
        {
            // Derivatives in local coords dN/dXi, dN/dEta
            // Midside nodes
            double[][] dnxe = new double[4][] { new double[2], new double[2], new double[2], new double[2] };
            
            // Corner nodes
            dnxe[0][0] = -0.25 * (1.0 - et);
            dnxe[0][1] = -0.25 * (1.0 - xi);
            dnxe[1][0] = 0.25 * (1.0 - et);
            dnxe[1][1] = -0.25 * (1.0 + xi);
            dnxe[2][0] = 0.25 * (1.0 + et);
            dnxe[2][1] = 0.25 * (1.0 + xi);
            dnxe[3][0] = -0.25 * (1.0 + et);
            dnxe[3][1] = 0.25 * (1.0 - xi);

            // Modification of derivatives due to degeneration
            int deg = degeneration(ind);
            if (deg > 0)
            {
                double sum0 = dnxe[deg - 1][0] + dnxe[deg][0];
                double sum1 = dnxe[deg - 1][1] + dnxe[deg][1];
                dnxe[deg - 1][0] = dnxe[deg][0] = sum0;
                dnxe[deg - 1][1] = dnxe[deg][1] = sum1;
            }

            // Jacobian matrix
            double[][] aj = new double[2][] { new double[2], new double[2] };
            for (int j = 0; j < 2; j++)
                for (int i = 0; i < 2; i++)
                {
                    aj[i][j] = 0.0;
                    for (int k = 0; k < 4; k++)
                        aj[i][j] += dnxe[k][j] * xy[k][i];
                }
            double det = aj[0][0] * aj[1][1] - aj[0][1] * aj[1][0];
            // Zero or negative determinant
            if (det <= 0)
                UTIL.errorMsg("Negative/zero Jacobian determinant for 8N element " + (float)det);
            // Jacobian inverse
            double aj00 = aj[1][1] / det;
            aj[1][1] = aj[0][0] / det;
            aj[0][0] = aj00;
            aj[1][0] = -aj[1][0] / det;
            aj[0][1] = -aj[0][1] / det;

            // Derivatives in global coordinates dN/dx, dN/dy
            for (int k = 0; k < 4; k++)
                for (int i = 0; i < 2; i++)
                    dnxy[k][i] = aj[0][i] * dnxe[k][0] + aj[1][i] * dnxe[k][1];

            return det;
        }

        // One-dimensional linear shape functions and
        //    their derivatives in local coordinates
        // xi - local coordinate;
        // an[2] - shape functions (out);
        // dndxi[3] - derivatives of shape functions (out)
        public static void shapeDerivFace(double xi, double[] an, double[] dndxi)
        {
            double x1 = 0.5 * (1.0 - xi);
            double x2 = 0.5 * (1.0 + xi);

            an[0] = x1;
            an[1] = x2;
            dndxi[0] = -0.5;
            dndxi[1] = 0.5;
        }
    }
}
