using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.util;

namespace FEA3D.elem
{
    // Quadratic 3D shape functions and their derivatives
    public class ShapeQuad3D
    {

        // Degeneration check.
        // The only degeneration is: 0=7=6,8=11,12=19=18
        public static int degeneration(int[] ind)
        {
            // Element should be quadratic
            if ((ind[0] == ind[7] && ind[7] == ind[6]) &&
                (ind[8] == ind[11]) &&
                (ind[12] == ind[19] && ind[19] == ind[18]))
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
            double s2 = 1 - xi * xi;
            double t2 = 1 - et * et;
            double d2 = 1 - ze * ze;
            // Midside nodes
            n[1] = n[3] = n[5] = n[7] = n[8] = n[9] = n[10] =
                   n[11] = n[13] = n[15] = n[17] = n[19] = 0;
            if (ind[1] > 0) n[1] = 0.25 * s2 * t1 * d1;
            if (ind[5] > 0) n[5] = 0.25 * s2 * t0 * d1;
            if (ind[17] > 0) n[17] = 0.25 * s2 * t0 * d0;
            if (ind[13] > 0) n[13] = 0.25 * s2 * t1 * d0;
            if (ind[7] > 0) n[7] = 0.25 * t2 * s1 * d1;
            if (ind[3] > 0) n[3] = 0.25 * t2 * s0 * d1;
            if (ind[15] > 0) n[15] = 0.25 * t2 * s0 * d0;
            if (ind[19] > 0) n[19] = 0.25 * t2 * s1 * d0;
            if (ind[8] > 0) n[8] = 0.25 * d2 * s1 * t1;
            if (ind[9] > 0) n[9] = 0.25 * d2 * s0 * t1;
            if (ind[10] > 0) n[10] = 0.25 * d2 * s0 * t0;
            if (ind[11] > 0) n[11] = 0.25 * d2 * s1 * t0;
            // Vertex nodes
            n[0] = 0.125 * s1 * t1 * d1 - 0.5 * (n[1] + n[7] + n[8]);
            n[2] = 0.125 * s0 * t1 * d1 - 0.5 * (n[1] + n[3] + n[9]);
            n[4] = 0.125 * s0 * t0 * d1 - 0.5 * (n[3] + n[5] + n[10]);
            n[6] = 0.125 * s1 * t0 * d1 - 0.5 * (n[5] + n[7] + n[11]);
            n[12] = 0.125 * s1 * t1 * d0 - 0.5 * (n[8] + n[13] + n[19]);
            n[14] = 0.125 * s0 * t1 * d0 - 0.5 * (n[9] + n[13] + n[15]);
            n[16] = 0.125 * s0 * t0 * d0 - 0.5 * (n[10] + n[15] + n[17]);
            n[18] = 0.125 * s1 * t0 * d0 - 0.5 * (n[11] + n[17] + n[19]);
            // Modification of functions due to degeneration
            if (degeneration(ind) == 1)
            {
                double dn1 = 0.0625 * s2 * t2 * d1;
                double dn2 = 0.0625 * s2 * t2 * d0;
                n[2] += dn1;
                n[3] -= 2 * dn1;
                n[4] += dn1;
                n[14] += dn2;
                n[15] -= 2 * dn2;
                n[16] += dn2;
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
            double[][] d = new double[20][];
            for(int i=0; i<20;i++)
                d[i]=new double[3];

            double s0 = 1 + xi;
            double t0 = 1 + et;
            double d0 = 1 + ze;
            double s1 = 1 - xi;
            double t1 = 1 - et;
            double d1 = 1 - ze;
            double s2 = 1 - xi*xi;
            double t2 = 1 - et*et;
            double d2 = 1 - ze*ze;
            // Midside nodes
            if (ind[1] > 0) { d[1][0] = -0.5*xi*t1*d1;
                              d[1][1] = -0.25*s2*d1;
                              d[1][2] = -0.25*s2*t1;
            } else { d[0][1] = d[1][1] = d[2][1] = 0;  }

            if (ind[5] > 0) { d[5][0] = -0.5*xi*t0*d1;
                              d[5][1] =  0.25*s2*d1;
                              d[5][2] = -0.25*s2*t0;
            } else { d[5][0] = d[5][1] = d[5][2] = 0;  }

            if (ind[17] > 0) { d[17][0] = -0.5*xi*t0*d0;
                               d[17][1] =  0.25*s2*d0;
                               d[17][2] =  0.25*s2*t0;
            } else {  d[17][0] = d[17][1] = d[17][2] = 0;  }

            if (ind[13] > 0) { d[13][0] = -0.5*xi*t1*d0;
                               d[13][1] = -0.25*s2*d0;
                               d[13][2] =  0.25*s2*t1;
            } else { d[13][0] = d[13][1] = d[13][2] = 0;  }

            if (ind[7] > 0) {  d[7][0] = -0.25*t2*d1;
                               d[7][1] = -0.5*et*s1*d1;
                               d[7][2] = -0.25*t2*s1;
            } else { d[7][0] = d[7][1] = d[7][2] = 0;  }

            if (ind[3] > 0) {  d[3][0] =  0.25*t2*d1;
                               d[3][1] = -0.5*et*s0*d1;
                               d[3][2] = -0.25*t2*s0;
            } else { d[3][0] = d[3][1] = d[3][2] = 0;  }

            if (ind[15] > 0) { d[15][0] =  0.25*t2*d0;
                               d[15][1] = -0.5*et*s0*d0;
                               d[15][2] =  0.25*t2*s0;
            } else { d[15][0] = d[15][1] = d[15][2] = 0;  }

            if (ind[19] > 0) { d[19][0] = -0.25*t2*d0;
                               d[19][1] = -0.5*et*s1*d0;
                               d[19][2] =  0.25*t2*s1;
            } else { d[19][0] = d[19][1] = d[19][2] = 0;  }

            if (ind[8] > 0) {  d[8][0] = -0.25*d2*t1;
                               d[8][1] = -0.25*d2*s1;
                               d[8][2] = -0.5*ze*s1*t1;
            } else { d[8][0] = d[8][1] = d[8][2] = 0;  }

            if (ind[9] > 0) {  d[9][0] =  0.25*d2*t1;
                               d[9][1] = -0.25*d2*s0;
                               d[9][2] = -0.5*ze*s0*t1;
            } else { d[9][0] = d[9][1] = d[9][2] = 0;  }

            if (ind[10] > 0) { d[10][0] =  0.25*d2*t0;
                               d[10][1] =  0.25*d2*s0;
                               d[10][2] = -0.5*ze*s0*t0;
            } else { d[10][0] = d[10][1] = d[10][2] = 0; }

            if (ind[11] > 0) { d[11][0] = -0.25*d2*t0;
                               d[11][1] =  0.25*d2*s1;
                               d[11][2] = -0.5*ze*s1*t0;
            } else { d[11][0] = d[11][1] = d[11][2] = 0; }
            // Vertex nodes
            d[0] [0]=-0.125*t1*d1-0.5*(d[1] [0]+d[7] [0]+d[8] [0]);
            d[0] [1]=-0.125*s1*d1-0.5*(d[1] [1]+d[7] [1]+d[8] [1]);
            d[0] [2]=-0.125*s1*t1-0.5*(d[1] [2]+d[7] [2]+d[8] [2]);

            d[2] [0]= 0.125*t1*d1-0.5*(d[1] [0]+d[3] [0]+d[9] [0]);
            d[2] [1]=-0.125*s0*d1-0.5*(d[1] [1]+d[3] [1]+d[9] [1]);
            d[2] [2]=-0.125*s0*t1-0.5*(d[1] [2]+d[3] [2]+d[9] [2]);

            d[4] [0]= 0.125*t0*d1-0.5*(d[3] [0]+d[5] [0]+d[10][0]);
            d[4] [1]= 0.125*s0*d1-0.5*(d[3] [1]+d[5] [1]+d[10][1]);
            d[4] [2]=-0.125*s0*t0-0.5*(d[3] [2]+d[5] [2]+d[10][2]);

            d[6] [0]=-0.125*t0*d1-0.5*(d[5] [0]+d[7] [0]+d[11][0]);
            d[6] [1]= 0.125*s1*d1-0.5*(d[5] [1]+d[7] [1]+d[11][1]);
            d[6] [2]=-0.125*s1*t0-0.5*(d[5] [2]+d[7] [2]+d[11][2]);

            d[12][0]=-0.125*t1*d0-0.5*(d[8] [0]+d[13][0]+d[19][0]);
            d[12][1]=-0.125*s1*d0-0.5*(d[8] [1]+d[13][1]+d[19][1]);
            d[12][2]= 0.125*s1*t1-0.5*(d[8] [2]+d[13][2]+d[19][2]);

            d[14][0]= 0.125*t1*d0-0.5*(d[9] [0]+d[13][0]+d[15][0]);
            d[14][1]=-0.125*s0*d0-0.5*(d[9] [1]+d[13][1]+d[15][1]);
            d[14][2]= 0.125*s0*t1-0.5*(d[9] [2]+d[13][2]+d[15][2]);

            d[16][0]= 0.125*t0*d0-0.5*(d[10][0]+d[15][0]+d[17][0]);
            d[16][1]= 0.125*s0*d0-0.5*(d[10][1]+d[15][1]+d[17][1]);
            d[16][2]= 0.125*s0*t0-0.5*(d[10][2]+d[15][2]+d[17][2]);

            d[18][0]=-0.125*t0*d0-0.5*(d[11][0]+d[17][0]+d[19][0]);
            d[18][1]= 0.125*s1*d0-0.5*(d[11][1]+d[17][1]+d[19][1]);
            d[18][2]= 0.125*s1*t0-0.5*(d[11][2]+d[17][2]+d[19][2]);
            // Modification of derivatives due to degeneration
            if (degeneration(ind) == 1) {
                double[] dn1 = new double[3];
                double[] dn2 = new double[3];
                dn1[0] = -0.125*xi*t2*d1;
                dn1[1] = -0.125*et*s2*d1;
                dn1[2] = -0.0625*s2*t2;
                dn2[0] = -0.125*xi*t2*d0;
                dn2[1] = -0.125*et*s2*d0;
                dn2[2] = -dn1[2];
                for (int i = 0; i < 3; i++) {
                    d[2] [i] += dn1[i];
                    d[3] [i] -= 2*dn1[i];
                    d[4] [i] += dn1[i];
                    d[14][i] += dn2[i];
                    d[15][i] -= 2*dn2[i];
                    d[16][i] += dn2[i];
                }
            }
            // Jacobian matrix ja
            double[][] ja = new double[3][]{new double[3], new double[3], new double[3]};
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    ja[j][i] = 0.0;
                    for (int k = 0; k < 20; k++)
                        ja[j][i] += d[k][i]*xy[k][j];
                }
            // Determinant of Jacobian matrix det
            double det =
                ja[0][0]*(ja[1][1]*ja[2][2] - ja[2][1]*ja[1][2])
              - ja[1][0]*(ja[0][1]*ja[2][2] - ja[0][2]*ja[2][1])
              + ja[2][0]*(ja[0][1]*ja[1][2] - ja[0][2]*ja[1][1]);
            if (det<=0) 
                UTIL.errorMsg("Negative/zero Jacobian determinant for 20Nr element "+(float)det);
            // Jacobian inverse ja1
            double[][] ja1 = new double[3][] { new double[3], new double[3], new double[3] };
            double v = 1.0/det;
            ja1[0][0] = (ja[1][1]*ja[2][2]-ja[2][1]*ja[1][2])*v;
            ja1[1][0] = (ja[0][2]*ja[2][1]-ja[0][1]*ja[2][2])*v;
            ja1[2][0] = (ja[0][1]*ja[1][2]-ja[0][2]*ja[1][1])*v;
            ja1[0][1] = (ja[1][2]*ja[2][0]-ja[1][0]*ja[2][2])*v;
            ja1[1][1] = (ja[0][0]*ja[2][2]-ja[2][0]*ja[0][2])*v;
            ja1[2][1] = (ja[0][2]*ja[1][0]-ja[0][0]*ja[1][2])*v;
            ja1[0][2] = (ja[2][1]*ja[1][0]-ja[2][0]*ja[1][1])*v;
            ja1[1][2] = (ja[0][1]*ja[2][0]-ja[0][0]*ja[2][1])*v;
            ja1[2][2] = (ja[0][0]*ja[1][1]-ja[1][0]*ja[0][1])*v;

            for (int k = 0; k < 20; k++)
                for (int i = 0; i < 3; i++) {
                    dnxy[k][i] = 0.0;
                    for (int j = 0; j < 3; j++)
                        dnxy[k][i] += ja1[i][j]*d[k][j];
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
            // Shape funcrions for midside nodes
            an[1] = an[3] = an[5] = an[7] = 0;
            if (ind[1]>0) an[1] =  x1*x2*e1*0.5;
            if (ind[3]>0) an[3] =  e1*e2*x2*0.5;
            if (ind[5]>0) an[5] =  x1*x2*e2*0.5;
            if (ind[7]>0) an[7] =  e1*e2*x1*0.5;
            // Shape functions for corner nodes
            an[0] = x1*e1*0.25 - 0.5*(an[7]+an[1]);
            an[2] = x2*e1*0.25 - 0.5*(an[1]+an[3]);
            an[4] = x2*e2*0.25 - 0.5*(an[3]+an[5]);
            an[6] = x1*e2*0.25 - 0.5*(an[5]+an[7]);

            dn[1][0] = dn[1][1] = dn[3][0] = dn[3][1] =
                dn[5][0] = dn[5][1] = dn[7][0] = dn[7][1] = 0;
            // Derivatives for midside nodes
            if (ind[1]>0) {
                dn[1][0] = -xi*e1; dn[1][1] = -0.5*x1*x2;}
            if (ind[3]>0) {
                dn[3][0] = 0.5*e1*e2;  dn[3][1] = -et*x2;}
            if (ind[5]>0) {
                dn[5][0] = -xi*e2;  dn[5][1] = 0.5*x1*x2;}
            if (ind[7]>0) {
                dn[7][0] = -0.5*e1*e2; dn[7][1] = -et*x1;}
            // Derivatives for corner nodes
            dn[0][0] = -0.25*e1 - 0.5*(dn[7][0]+dn[1][0]);
            dn[0][1] = -0.25*x1 - 0.5*(dn[7][1]+dn[1][1]);
            dn[2][0] =  0.25*e1 - 0.5*(dn[1][0]+dn[3][0]);
            dn[2][1] = -0.25*x2 - 0.5*(dn[1][1]+dn[3][1]);
            dn[4][0] =  0.25*e2 - 0.5*(dn[3][0]+dn[5][0]);
            dn[4][1] =  0.25*x2 - 0.5*(dn[3][1]+dn[5][1]);
            dn[6][0] = -0.25*e2 - 0.5*(dn[5][0]+dn[7][0]);
            dn[6][1] =  0.25*x1 - 0.5*(dn[5][1]+dn[7][1]);

            // Degeneration check
            int ig = 0;
            for (int i = 0; i < 7; i += 2) {
                if (ind[i] == ind[i + 1]) {
                    ig = (i + 5) % 8;
                    break;
                }
            }
            if (ig>0 && ind[1]>0 && ind[3]>0 && ind[5]>0 && ind[7]>0){
                double delta = 0.125*x1*x2*e1*e2;
                double z = -0.25*xi*e1*e2;
                double t = -0.25*x1*x2*et;
                int j = (ig + 1)%8;
                an[ig-1] +=  delta;
                an[ig]   -=  2.0*delta;
                an[j]    +=  delta;
                dn[ig-1][0] += z;
                dn[ig-1][1] += t;
                dn[ig][0]   -= 2*z;
                dn[ig][1]   -= 2*t;
                dn[j][0]    += z;
                dn[j][1]    += t;
            }
        }

    }
}
