using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.util;

namespace FEA3D.elem
{
    // Quadratic 2D shape functions and their derivatives
    public class ShapeQuad2D
    {
        // Degeneration check.
        // If element is triangular then the method return
        // a local number (starting from 0) of the midside node
        // opposite to degenerated side.
        // ind - connectivity numbers
        static int degeneration(int[] ind)
        {
            int deg = 0;
            for (int i = 0; i < 7; i += 2) {
                if (ind[i] == ind[i + 1]) {
                    deg = (i + 5) % 8;
                    break;
                }
            }
            return deg;
        }

        // Shape functions.
        // xi, et - local coordinates;
        // ind[8] - element connectivities;
        // an[8] - shape functions (out)
        public static void shape(double xi, double et, int[] ind, double[] an)
        {
            // Shape functions of midside nodes
            an[1] = an[3] = an[5] = an[7] = 0.0;
            if (ind[1] > 0) an[1] = 0.5 * (1.0 - xi * xi) * (1.0 - et);
            if (ind[3] > 0) an[3] = 0.5 * (1.0 - et * et) * (1.0 + xi);
            if (ind[5] > 0) an[5] = 0.5 * (1.0 - xi * xi) * (1.0 + et);
            if (ind[7] > 0) an[7] = 0.5 * (1.0 - et * et) * (1.0 - xi);

            // Shape functions of corner nodes
            an[0] = 0.25 * (1.0 - xi) * (1.0 - et) - 0.5 * (an[7] + an[1]);
            an[2] = 0.25 * (1.0 + xi) * (1.0 - et) - 0.5 * (an[1] + an[3]);
            an[4] = 0.25 * (1.0 + xi) * (1.0 + et) - 0.5 * (an[3] + an[5]);
            an[6] = 0.25 * (1.0 - xi) * (1.0 + et) - 0.5 * (an[5] + an[7]);

            // Modification of functions due to degeneration
            int deg = degeneration(ind);
            if (deg > 0 && ind[1] > 0 && ind[3] > 0 && ind[5] > 0 && ind[7] > 0) {
                double delta = 0.125*(1 - xi*xi)*(1 - et*et);
                an[deg - 1] += delta;
                an[deg] -= 2.0*delta;
                an[(deg + 1)%8] += delta;
            }
        }

        // Derivatives of shape functions
        // with respect to global coordinates x and y.
        // xi, et - local coordinates;
        // ind[8] - element connectivities;
        // xy[8][2] - nodal coordinates;
        // dnxy[8][2] - derivatives of shape functions (out);
        // returns  determinant of the Jacobian matrrix
        public static double deriv(double xi, double et, int[] ind, double[][] xy, double[][] dnxy) 
         {
            // Derivatives in local coords dN/dXi, dN/dEta
            // Midside nodes
            double[][] dnxe = new double[8][]{new double[2], new double[2], new double[2], new double[2], new double[2], new double[2], new double[2], new double[2]};
            dnxe[1][0] = dnxe[1][1] = dnxe[3][0] = dnxe[3][1] = dnxe[5][0] = dnxe[5][1] = dnxe[7][0] = dnxe[7][1] = 0.0;
            if (ind[1] > 0) {
                dnxe[1][0] = -xi * (1.0 - et);
                dnxe[1][1] = -0.5 * (1.0 - xi * xi);
            }
            if (ind[3] > 0) {
                dnxe[3][0] = 0.5 * (1.0 - et * et);
                dnxe[3][1] = -et * (1.0 + xi);
            }
            if (ind[5] > 0) {
                dnxe[5][0] = -xi * (1.0 + et);
                dnxe[5][1] = 0.5 * (1.0 - xi * xi);
            }
            if (ind[7] > 0) {
                dnxe[7][0] = -0.5*(1.0-et*et);
                dnxe[7][1] = -et*(1.0-xi);
            }
            // Corner nodes
            dnxe[0][0] = -0.25 * (1.0 - et) - 0.5 * (dnxe[7][0] + dnxe[1][0]);
            dnxe[0][1] = -0.25 * (1.0 - xi) - 0.5 * (dnxe[7][1] + dnxe[1][1]);
            dnxe[2][0] = 0.25 * (1.0 - et) - 0.5 * (dnxe[1][0] + dnxe[3][0]);
            dnxe[2][1] = -0.25 * (1.0 + xi) - 0.5 * (dnxe[1][1] + dnxe[3][1]);
            dnxe[4][0] = 0.25 * (1.0 + et) - 0.5 * (dnxe[3][0] + dnxe[5][0]);
            dnxe[4][1] = 0.25 * (1.0 + xi) - 0.5 * (dnxe[3][1] + dnxe[5][1]);
            dnxe[6][0] = -0.25 * (1.0 + et) - 0.5 * (dnxe[5][0] + dnxe[7][0]);
            dnxe[6][1] = 0.25 * (1.0 - xi) - 0.5 * (dnxe[5][1] + dnxe[7][1]);

            // Modification of derivatives due to degeneration
            int deg = degeneration(ind);
            if (deg > 0 && ind[1] > 0 && ind[3] > 0 && ind[5] > 0 && ind[7] > 0) {
                double z = -0.25*xi*(1-et*et);
                double t = -0.25*(1-xi*xi)*et;
                int j = (deg + 1)%8;
                dnxe[deg - 1][0] = dnxe[deg - 1][0] + z;
                dnxe[deg - 1][1] = dnxe[deg - 1][1] + t;
                dnxe[deg][0] = dnxe[deg][0] - 2.0*z;
                dnxe[deg][1] = dnxe[deg][1] - 2.0*t;
                dnxe[j][0] = dnxe[j][0] + z;
                dnxe[j][1] = dnxe[j][1] + t;
            }

            // Jacobian matrix
            double[][] aj = new double[2][]{new double[2], new double[2]};
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 2; i++) {
                    aj[i][j] = 0.0;
                    for (int k = 0; k < 8; k++)
                        aj[i][j] += dnxe[k][j]*xy[k][i];
                }
            }
            double det = aj[0][0]*aj[1][1] - aj[0][1]*aj[1][0];
            // Zero or negative determinant
            if (det<=0) 
                UTIL.errorMsg("Negative/zero Jacobian determinant for 8N element "+(float)det);
            // Jacobian inverse
            double aj00 = aj[1][1]/det;
            aj[1][1] =  aj[0][0]/det;
            aj[0][0] =  aj00;
            aj[1][0] = -aj[1][0]/det;
            aj[0][1] = -aj[0][1]/det;

            // Derivatives in global coordinates dN/dx, dN/dy
            for (int k = 0; k < 8; k++)
                for (int i = 0; i < 2; i++)
                    dnxy[k][i] = aj[0][i]*dnxe[k][0]
                               + aj[1][i]*dnxe[k][1];
            return det;
        }

        // One-dimensional quadratic shape functions and
        //    their derivatives in local coordinates
        // xi - local coordinate;
        // kmid - index of midside node (=0 no midside node);
        // an[3] - shape functions (out);
        // dndxi[3] - derivatives of shape functions (out)
        public static void shapeDerivFace(double xi, int kmid, double[] an, double[] dndxi)
        {
            double x1 = 1 - xi;
            double x2 = 1 + xi;
            if (kmid > 0)
            {
                an[1] = x1 * x2;
                dndxi[1] = -2.0 * xi;
            }
            an[0] = 0.5 * x1 - 0.5 * an[1];
            an[2] = 0.5 * x2 - 0.5 * an[1];
            dndxi[0] = -0.5 - 0.5 * dndxi[1];
            dndxi[2] = 0.5 - 0.5 * dndxi[1];
        }
    }
}
