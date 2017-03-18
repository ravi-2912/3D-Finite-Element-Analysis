using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.util;
using FEA3D.model;
using FEA3D.material;

namespace FEA3D.elem
{
    class ElemTruss3d2 : Element
    {
        // Element edges (local numbers)
        private static int[][] faceInd = new int[1][] { new int[2] { 0, 1 } };
        // Shape functions
        private static double[] an = new double[2];
        // Derivatives of shape functions
        private static double[] dnxy = new double[2] { 0, 0 }; 
        // Displacements differentiation matrix
        private static double[] bmat = new double[2] { 0, 0 };
        // Elasticity matrix
        private static double[][] emat = new double[4][] { new double[4], new double[4], new double[4], new double[4] };
        // Thermal strains
        private static double ept = 0.0;
        // 3d transformation matrix
        private static double[][] tmat = new double[2][] { new double[6] { 0, 0, 0, 0, 0, 0 }, new double[6] { 0, 0, 0, 0, 0, 0 } };
        // Truss Length
        private static double L;
        // Gauss rules for stiffness matrix, thermal vector,
        // surface load and stress integration
        private static GaussRule gk = new GaussRule(1, 1); 
        private static GaussRule gh = new GaussRule(2, 1);
        private static GaussRule gf = new GaussRule(2, 1);
        private static GaussRule gs = new GaussRule(1, 1);

        // Constructor for truss32  element
        public ElemTruss3d2()
            : base("truss32", 2, 1) 
        {
            // Constructor code
        }

        // Compute stiffness matrix
        public override void stiffnessMatrix()
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0) + Math.Pow((xy[1][2] - xy[0][2]), 2.0));            
            tmat[0][0] = tmat[1][3] = (xy[1][0] - xy[0][0]) / L; // lij
            tmat[0][1] = tmat[1][4] = (xy[1][1] - xy[0][1]) / L; // mij
            tmat[0][2] = tmat[1][5] = (xy[1][2] - xy[0][2]) / L; // nij

            double[][] kmat_lcl = new double[2][] { new double[2] { 0, 0 }, new double[2] { 0, 0 } };

            // Material mat
            mat = (Material)fem.materials[matName];
            if (mat == null)
                UTIL.errorMsg("Element material name: " + matName);
            //mat.elasticityMatrix(emat);

            // Gauss integration loop
            for (int ip = 0; ip < gk.nIntPoints; ip++)
            {
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gk.xii[ip]);
                double dv = det * gk.wi[ip];
                kmat_lcl[0][0] = bmat[0] * mat.getElasticModulus() * bmat[0] * this.A * dv;
                kmat_lcl[0][1] = bmat[0] * mat.getElasticModulus() * bmat[1] * this.A * dv;
                kmat_lcl[1][0] = bmat[1] * mat.getElasticModulus() * bmat[0] * this.A * dv;
                kmat_lcl[1][1] = bmat[1] * mat.getElasticModulus() * bmat[1] * this.A * dv;
            }

            // Apply Transformation and create upper symmetrical part of the stiffness matrix
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    kmat[i][j] = 0.0; // Zero to stiffness matrix
                    for (int k = 0; k < 2; k++)
                        for (int l = 0; l < 2; l++)
                            kmat[i][j] += tmat[l][i] * kmat_lcl[l][k] * tmat[k][j];
                }
        
        }

        // Set displacement differentiation matrix bmat.
        // xi, et - local coordinates,
        // returns  determinant of Jacobian matrix
        private double setBmatrix(double xi) 
        {
            // Derivatives of shape functions
            double det = deriv(xi, xy, dnxy);

            if (det <= 0) 
                UTIL.errorMsg("Negative/zero 2N element area");

            bmat[0] = dnxy[0];
            bmat[1] = dnxy[1];
            return det;
        }

        // Compute thermal vector
        public override void thermalVector()
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0) + Math.Pow((xy[1][2] - xy[0][2]), 2.0));
            tmat[0][0] = tmat[1][3] = (xy[1][0] - xy[0][0]) / L; // lij
            tmat[0][1] = tmat[1][4] = (xy[1][1] - xy[0][1]) / L; // mij
            tmat[0][2] = tmat[1][5] = (xy[1][2] - xy[0][2]) / L; // nij

            double[] evec_lcl = new double[2] { 0, 0 };

            // Material mat
            mat = (Material)fem.materials[matName];
            double alpha = mat.getAlpha();
            double nu = mat.getNu();

            // Gauss integration loop
            for (int ip = 0; ip < gh.nIntPoints; ip++)
            {
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gh.xii[ip]);
                // Shape functions an
                shape(gh.xii[ip], an);
                double t = an[0] * dtn[0] + an[1] * dtn[1];
                double dv = det * gh.wi[ip];
                ept = alpha * t;

                evec_lcl[0] = bmat[0] * mat.getElasticModulus() * ept * dv;
                evec_lcl[1] = bmat[1] * mat.getElasticModulus() * ept * dv;

            }

            // Apply Transform
            for (int i = 0; i < 6; i++)
                evec[i] = tmat[0][i] * evec_lcl[0] + tmat[1][i] * evec_lcl[1];

        }

        // Set nodal equivalent of distributed face load to evec.
        // surLd - object describing element face load;
        // returns loaded element face
        // or -1 (loaded nodes do not match elem face)
        public override int equivFaceLoad(ElemFaceLoad surLd)
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0) + Math.Pow((xy[1][2] - xy[0][2]), 2.0));
            tmat[0][0] = tmat[1][3] = (xy[1][0] - xy[0][0]) / L; // lij
            tmat[0][1] = tmat[1][4] = (xy[1][1] - xy[0][1]) / L; // mij
            tmat[0][2] = tmat[1][5] = (xy[1][2] - xy[0][2]) / L; // nij

            // Shape functons
            double[] an = new double[2]; 
            // Derivatives of shape functions
            double[] xin = new double[2];

            double[] evec_lcl = new double[2] { 0, 0 };

            int loadedFace = surLd.rearrange(faceInd, ind);
            if (loadedFace == -1) return -1;

            // Gauss integration loop
            for (int ip = 0; ip < gf.nIntPoints; ip++)
            {
                shapeDerivFace(gf.xii[ip], an, xin);
                double p = 0.0, xs = 0.0, ys = 0.0, zs = 0.0;
                for (int i = 0; i < 2; i++)
                {
                    p += an[i] * surLd.forceAtNodes[i];
                    int j = faceInd[loadedFace][i];
                    xs += xin[i] * xy[j][0];
                    ys += xin[i] * xy[j][1];
                    zs += xin[i] * xy[j][2];
                }
                double ds = Math.Sqrt(xs * xs + ys * ys + zs * zs);
                for (int i = 0; i < 2; i++)
                {
                    int j = faceInd[loadedFace][i];
                    evec_lcl[j] += an[i] * p * ds * gf.wi[ip];
                    //evec_lcl[j + 1] += an[i] * p * ds * gf.wi[ip];
                }
            }

            // Apply Transform
            for (int i = 0; i < 6; i++)
                evec[i] = tmat[0][i] * evec_lcl[0] + tmat[1][i] * evec_lcl[1];

            return loadedFace;
        }

        // Compute equivalent stress vector (with negative sign)
        // FOR ELASTIC-PLASTIC
        public override void equivStressVector()
        {
            Console.WriteLine("EquivStressVector - verify code");
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0) + Math.Pow((xy[1][2] - xy[0][2]), 2.0));
            tmat[0][0] = tmat[1][3] = (xy[1][0] - xy[0][0]) / L; // lij
            tmat[0][1] = tmat[1][4] = (xy[1][1] - xy[0][1]) / L; // mij
            tmat[0][2] = tmat[1][5] = (xy[1][2] - xy[0][2]) / L; // nij

            double[] evec_lcl = new double[2] { 0, 0 };

            for (int ip = 0; ip < gs.nIntPoints; ip++)
            {
                // Accumulated stress
                double s = str[ip].sStress[0] + str[ip].dStress[0];
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gs.xii[ip]);
                double dv = det * gs.wi[ip];

                evec_lcl[0] = bmat[0] * s * dv;
                evec_lcl[1] = bmat[1] * s * dv;

            }
            // Apply Transform
            for (int i = 0; i < 6; i++)
                evec[i] = tmat[0][i] * evec_lcl[0] + tmat[1][i] * evec_lcl[1];
        }


        // Extrapolate values from integration points to nodes.
        // fip [1][1] - values at integration points;
        // fn [2][3] - values at nodes (out)
        public override void extrapolateToNodes(double[][] fip, double[][] fn)
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0) + Math.Pow((xy[1][2] - xy[0][2]), 2.0));
            tmat[0][0] = tmat[1][3] = (xy[1][0] - xy[0][0]) / L; // lij
            tmat[0][1] = tmat[1][4] = (xy[1][1] - xy[0][1]) / L; // mij
            tmat[0][2] = tmat[1][5] = (xy[1][2] - xy[0][2]) / L; // nij

            // Extrapolation matrix
            double lim = 1; // for 1 intwgration point

            double[][] fn_lcl = new double[2][] { new double[1] { 0 }, new double[1] { 0 } };

            for (int i = 0; i < 2; i ++)
                for (int j = 0; j < 1; j++)
                    fn[i][j] = 0;

            fn_lcl[0][0] = lim * fip[0][0];
            fn_lcl[1][0] = lim * fip[0][0];

            // Apply Transform
            int jind = 0;
            for(int c = 0; c < 2; c++)
                for(int s = 0; s < 3; s++)
                {
                    fn[c][s] = tmat[0][jind] * fn_lcl[0][0] + tmat[1][jind] * fn_lcl[1][0];
                    jind++;
                }

        }

        // Get local node numbers for element faces.
        // returns elementFaces[nFaces][nNodesOnFace]
        public override int[][] getElemFaces()
        {
            return faceInd;
        }

        // Get strains at integration point.
        // ip - integration point number (stress);
        // returns  strain vector (ex, ey, gxy, ez)
        public override double[] getStrainsAtIntPoint(int ip)
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0) + Math.Pow((xy[1][2] - xy[0][2]), 2.0));
            tmat[0][0] = tmat[1][3] = (xy[1][0] - xy[0][0]) / L; // lij
            tmat[0][1] = tmat[1][4] = (xy[1][1] - xy[0][1]) / L; // mij
            tmat[0][2] = tmat[1][5] = (xy[1][2] - xy[0][2]) / L; // nij

            double[] evec_lcl = new double[2] { 0, 0 };

            // Apply Transform
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 6; j++)
                {
                    evec_lcl[i] += tmat[i][j] * evec[j];                    
                }
            setBmatrix(gs.xii[ip]);
            double[] strain = new double[] { 0, 0, 0, 0, 0, 0 };
            strain[0] = bmat[0] * evec_lcl[0] + bmat[1] * evec_lcl[1];
            
            return strain;
        }

        // Get temperature at integration point (stress)
        public override double getTemperatureAtIntPoint(int ip)
        {
            shape(gs.xii[ip], an);
            return (an[0] * dtn[0] + an[1] * dtn[1]);
        }

        // Shape functions.
        // xi - local coordinates;
        // an[2] - shape functions (out)
        public static void shape(double xi, double[] an)
        {
            // Shape functions of end nodes
            an[0] = 0.5 * (1.0 - xi);
            an[1] = 0.5 * (1.0 + xi) ;
        }

        // Derivatives of shape functions with respect to global coordinates x.
        // xi - local coordinates;
        // xy[2][3] - nodal coordinates;
        // dnxy[2][1] - derivatives of shape functions (out);
        // returns  determinant of the Jacobian matrrix
        public static double deriv(double xi, double[][] xy, double[] dnxy)
        {
            // Derivatives in local coords dN/dXi
            double dnxe1 = -0.5;
            double dnxe2 = 0.5;
            // Jacobian 
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0) + Math.Pow((xy[1][2] - xy[0][2]), 2.0));            
            double det = L / 2.0;
            // Zero or negative determinant
            if (det <= 0)
                UTIL.errorMsg("Negative/zero Jacobian determinant for 2N element " + (float)det);
            // Jacobian inverse
            double aj00 = 1.0 / det;

            // Derivatives in global coordinates dN/dx
            dnxy[0] = aj00 * dnxe1;
            dnxy[1] = aj00 * dnxe2;
            
            return det;
        }

        // One-dimensional linear shape functions and
        //    their derivatives in local coordinates
        // xi - local coordinate;
        // an[2] - shape functions (out);
        // dndxi[2] - derivatives of shape functions (out)
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
