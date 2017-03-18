using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;

using FEA3D.util;
//using FEA3D.model;
using FEA3D.material;

namespace FEA3D.elem
{
    class ElemEBbeam2d : Element
    {
        // Element edges (local numbers)
        private static int[][] faceInd = new int[1][] { new int[2] { 0, 1 } };
        // Shape functions
        private static double[] an = new double[4];
        // Derivative of shape function in local dN/dEta
        private static double[][] dnxe = new double[1][] { new double[4] };
        // Derivatives of shape functions
        private static double[][] dnxy = new double[1][] { new double[4] };
        // Double derivate of shape function in local d2N/dEta2
        private static double[][] d2nxe = new double[1][] { new double[4] };
        // Double derivatives of shape function
        private static double[][] d2nxy = new double[1][] { new double[4] };
        // Displacements differentiation matrix
        private static double[][] bmat = new double[1][] { new double[4] };
        // Elasticity matrix
        private static double[][] emat = new double[4][] { new double[4], new double[4], new double[4], new double[4] };
        // Thermal strains
        private static double[] ept = new double[1];
        private static double[][] tmat = new double[2][] { new double[4] { 0, 0, 0, 0 }, new double[4] { 0, 0, 0, 0 } };
        private static double L;
        // Gauss rules for stiffness matrix, thermal vector,
        // surface load and stress integration
        private static GaussRule gk = new GaussRule(1, 1);
        private static GaussRule gh = new GaussRule(2, 1);
        private static GaussRule gf = new GaussRule(2, 1);
        private static GaussRule gs = new GaussRule(1, 1);

        // Constructor for truss22  element
        public ElemEBbeam2d()
            : base("ebbm22", 2, 2)
        {
            //super ("quad8", 8, 4);
            //base ("ebbm22",2,2);
        }

        // Compute stiffness matrix
        public override void stiffnessMatrix()
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));
            double lij = (xy[1][0] - xy[0][0]) / L;
            double mij = (xy[1][1] - xy[0][1]) / L;
            tmat[0][0] = tmat[1][2] = lij;
            tmat[0][1] = tmat[1][3] = mij;

            // Zeros to stiffness matrix kmat
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    kmat[i][j] = 0.0;
            double[][] kmat_lcl = new double[6][];
            for (int i = 0; i < 6; i++)
                kmat_lcl[i] = new double[6] { 0, 0, 0, 0, 0, 0 };

            // ld = length of strain/stress vector (1)
            int ld = 1;
            // Material mat
            mat = (Material)fem.materials[matName];
            //Console.WriteLine(matName);
            if (mat == null)
                UTIL.errorMsg("Element material name: " + matName);
            //mat.elasticityMatrix(emat);

            // Gauss integration loop
            for (int ip = 0; ip < gk.nIntPoints; ip++)
            {
                // Set displacement differentiation matrix bmat
                double det = 0.0;//setBmatrix(gk.xii[ip]);
                deriv2(gk.xii[ip], xy, d2nxy);
                double dv = det * gk.wi[ip];
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                    {
                        double s = 0.0;
                        for (int k = 0; k < ld; k++)
                            for (int l = 0; l < ld; l++)
                                s += bmat[l][i] * mat.getElasticModulus() * bmat[k][j] * this.A;

                        kmat_lcl[i][j] += s * dv;
                    }
            }
        }



        public static void shape(double et, double[][] xy, double[] an)
        {

            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));
            // Shape functions of end nodes
            an[0] = 0.25 * (2.0 - 3.0 * et + et * et * et);
            an[1] = 0.125 * L * (1.0 - et - et * et + et * et * et);
            an[2] = 0.25 * (2.0 + 3.0 * et - et * et * et);
            an[3] = 0.125 * L * (-1.0 - et + et * et + et * et * et);
        }

        public static double deriv(double et, double[][] xy, double[][] dnxy)
        {
            // Jacobian matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));
            double aj = L / 2.0;
            double det = aj;
            // Zero or negative determinant
            if (det <= 0)
                UTIL.errorMsg("Negative/zero Jacobian determinant for 2N element " + (float)det);
            // Jacobian inverse
            double aj00 = 1.0 / det;


            // End nodes
            dnxe[0][0] = 0.25 * (-3.0 + 3.0 * et * et);
            dnxe[0][1] = 0.125 * L * (1.0 - 2.0 * et + 3.0 * et * et);
            dnxe[0][2] = 0.25 * (3.0 - 3.0 * et * et);
            dnxe[0][3] = 0.125 * L * (-1.0 + 2.0 * et + 3.0 * et * et);

            // Derivatives in global coordinates dN/dx
            dnxy[0][0] = aj00 * dnxe[0][0];
            dnxy[0][1] = aj00 * dnxe[0][1];
            dnxy[0][2] = aj00 * dnxe[0][2];
            dnxe[0][3] = aj00 * dnxe[0][3];

            return det;
        }

        public static void deriv2(double et, double[][] xy, double[][] d2nxy)
        {
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));

            // double derivatives in global coordinates dN/dx
            d2nxy[0][0] = 6.0 * et / (L * L);
            d2nxy[0][1] = (-1.0 + 3.0 * et) / L;
            d2nxy[0][2] = -6.0 * et / (L * L);
            d2nxy[0][3] = (1.0 + 3.0 * et) / L;


        }

        public static void shapeDerivFace(double et, double[] an, double[] dndxi)
        {
            an[0] = 0.5 * (1.0 - et);
            an[1] = 0.5 * (1.0 + et);
            dndxi[0] = -0.5;
            dndxi[1] = 0.5;
        }

    }
}




