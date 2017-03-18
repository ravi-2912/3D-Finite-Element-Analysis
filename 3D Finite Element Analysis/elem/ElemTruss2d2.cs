using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.util;
using FEA3D.model;
using FEA3D.material;

namespace FEA3D.elem
{
    public class ElemTruss2d2 : Element
    {
        // Element edges (local numbers)
        private static int[][] faceInd = new int[1][] { new int[2] { 0, 1 } };
        // Shape functions
        private static double[] an = new double[2];
        // Derivatives of shape functions
        private static double[][] dnxy = new double[1][] { new double[2] };      
        // Displacements differentiation matrix
        private static double[][] bmat = new double[1][] { new double[2] };
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
        public ElemTruss2d2()
            : base("truss22", 2, 1) 
        {
            //super ("quad8", 8, 4);
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
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    kmat[i][j] = 0.0;
            double[][] kmat_lcl = new double[2][] { new double[2] { 0, 0 }, new double[2] { 0, 0 } };

            // ld = length of strain/stress vector (1)
            int ld = 1;
            // Material mat
            mat = (Material)fem.materials[matName];
            //Console.WriteLine(matName);
            if (mat == null)
                UTIL.errorMsg("Element material name: " + matName);
            mat.elasticityMatrix(emat);

            // Gauss integration loop
            for (int ip = 0; ip < gk.nIntPoints; ip++)
            {
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gk.xii[ip]);
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

            // Apply Transformation and create upper symmetrical part of the stiffness matrix
            for (int i = 0; i < 4; i++)
                for (int j = i; j < 4; j++)
                    for (int k = 0; k < 2; k++)
                        for (int l = 0; l < 2; l++)
                            kmat[i][j] += tmat[l][i] * kmat_lcl[l][k] * tmat[k][j];
        
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

            bmat[0][0] = dnxy[0][0];
            bmat[0][1] = dnxy[0][1];
            return det;
        }

        // Compute thermal vector
        public override void thermalVector()
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));            
            double lij = (xy[1][0] - xy[0][0]) / L;
            double mij = (xy[1][1] - xy[0][1]) / L;
            tmat[0][0] = tmat[1][2] = lij;
            tmat[0][1] = tmat[1][3] = mij;

            // Zeros to thermal vector evec
            for (int i = 0; i < 4; i++) evec[i] = 0.0;
            double[] evec_lcl = new double[2] { 0, 0 };

            int ld = 1;
            // Material mat
            mat = (Material)fem.materials[matName];
            mat.elasticityMatrix(emat);
            double alpha = mat.getAlpha();
            double nu = mat.getNu();

            // Gauss integration loop
            for (int ip = 0; ip < gh.nIntPoints; ip++)
            {
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gh.xii[ip]);
                // Shape functions an
                shape(gh.xii[ip], an);
                double t = 0.0;
                for (int i = 0; i < 2; i++) 
                    t += an[i] * dtn[i];
                double dv = det * gh.wi[ip];
                ept[0] = alpha * t;
                
                for (int i = 0; i < 2; i++)
                {
                    double s = 0;
                    for (int j = 0; j < ld; j++)
                        for (int k = 0; k < ld; k++)
                            s += bmat[k][i] * mat.getElasticModulus() * ept[j];
                    evec_lcl[i] += s * dv;
                }
            }

            // Apply Transform
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 2; j++)
                    evec[i] += tmat[j][i] * evec_lcl[j];
        }

        // Set nodal equivalent of distributed face load to evec.
        // surLd - object describing element face load;
        // returns loaded element face
        // or -1 (loaded nodes do not match elem face)
        public override int equivFaceLoad(ElemFaceLoad surLd)
        {
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));            
            double lij = (xy[1][0] - xy[0][0]) / L;
            double mij = (xy[1][1] - xy[0][1]) / L;
            tmat[0][0] = tmat[1][2] = lij;
            tmat[0][1] = tmat[1][3] = mij;

            // Shape functons
            double[] an = new double[2]; 
            // Derivatives of shape functions
            double[] xin = new double[2];

            double[] evec_lcl = new double[2] { 0, 0 };

            for (int i = 0; i < 4; i++) evec[i] = 0.0;

            int loadedFace = surLd.rearrange(faceInd, ind);
            if (loadedFace == -1) return -1;

            // Gauss integration loop
            for (int ip = 0; ip < gf.nIntPoints; ip++)
            {
                shapeDerivFace(gf.xii[ip], an, xin);
                double p = 0.0, xs = 0.0, ys = 0.0;
                for (int i = 0; i < 2; i++)
                {
                    p += an[i] * surLd.forceAtNodes[i];
                    int j = faceInd[loadedFace][i];
                    xs += xin[i] * xy[j][0];
                    ys += xin[i] * xy[j][1];
                }
                double ds = Math.Sqrt(xs * xs + ys * ys);
                for (int i = 0; i < 2; i++)
                {
                    int j = faceInd[loadedFace][i];
                    evec_lcl[j] += an[i] * p * ds * gf.wi[ip];
                    //evec_lcl[j + 1] += an[i] * p * ds * gf.wi[ip];
                }
            }

            // Apply Transform
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 2; j++)
                    evec[i] += tmat[j][i] * evec_lcl[j];

            return loadedFace;
        }

        // Compute equivalent stress vector (with negative sign)
        // FOR ELASTIC-PLASTIC
        public override void equivStressVector()
        {
            Console.WriteLine("EquivStressVector - verify code");
            // Transformation matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));            
            double lij = (xy[1][0] - xy[0][0]) / L;
            double mij = (xy[1][1] - xy[0][1]) / L;
            tmat[0][0] = tmat[1][2] = lij;
            tmat[0][1] = tmat[1][3] = mij;

            for (int i = 0; i < 4; i++) evec[i] = 0.0;
            double[] evec_lcl = new double[2] { 0, 0 };
            int ld = 1;

            for (int ip = 0; ip < gs.nIntPoints; ip++)
            {
                // Accumulated stress
                double[] s = new double[1]; // 4 for full integration
                for (int i = 0; i < 1; i++) // 4 fpr full integration
                    s[i] = str[ip].sStress[i] + str[ip].dStress[i];
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gs.xii[ip]);
                double dv = det * gs.wi[ip];

                for (int i = 0; i < 2; i++)
                {
                    double a = 0;
                    for (int j = 0; j < ld; j++)
                        a += bmat[j][i] * s[j];
                    evec_lcl[i] -= a * dv;
                }
            }
            // Apply Transform
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 2; j++)
                    evec[i] += tmat[j][i] * evec_lcl[i];
        }


        // Extrapolate values from integration points to nodes.
        // fip [1][1] - values at integration points;
        // fn [2][1] - values at nodes (out)
        public override void extrapolateToNodes(double[][] fip, double[][] fn)
        {
            // Extrapolation matrix
            double lim = 1;
            double[][] fn_lcl = new double[2][] { new double[1] { 0 }, new double[1] { 0 } };

            for (int i = 0; i < 2; i ++)
                for (int j = 0; j < 1; j++)
                    fn[i][j] = 0;

            for (int corner = 0; corner < 2; corner++)
            {
                for (int k = 0; k < 1; k++)
                {
                    double c = 0.0;
                    for (int ip = 0; ip < 1; ip++)
                        c += lim * fip[ip][k];
                    fn_lcl[corner][k] = c;
                }                
            }
            for (int c = 0; c < 2; c++)
                Console.WriteLine(fn_lcl[c][0]);
                // Transformation matrix
                L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));
            double lij = (xy[1][0] - xy[0][0]) / L;
            double mij = (xy[1][1] - xy[0][1]) / L;
            tmat[0][0] = tmat[1][2] = lij;
            tmat[0][1] = tmat[1][3] = mij;

            // Apply Transform
            int jind = 0;
            for(int c =0; c<2; c++)
                for(int s = 0; s<2; s++)
                {
                    for (int i = 0; i < 2; i++)
                        fn[c][s] += tmat[i][jind] * fn_lcl[i][0];
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
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));            
            double lij = (xy[1][0] - xy[0][0]) / L;
            double mij = (xy[1][1] - xy[0][1]) / L;
            tmat[0][0] = tmat[1][2] = lij;
            tmat[0][1] = tmat[1][3] = mij;

            double[] evec_lcl = new double[2] { 0, 0 };

            // Apply Transform
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 4; j++)
                {
                    evec_lcl[i] += tmat[i][j] * evec[j];                    
                }
            
            // Set displacement differentiation matrix bmat
            setBmatrix(gs.xii[ip]);
            double[] strain = new double[4] { 0, 0, 0, 0 };
            for (int i = 0; i < 1; i++)
            {
                strain[i] = 0;
                for (int j = 0; j < 2; j++)
                    strain[i] += bmat[i][j] * evec_lcl[j];
            }
            
            return strain;
        }

        // Get temperature at integration point (stress)
        public override double getTemperatureAtIntPoint(int ip)
        {
            shape(gs.xii[ip], an);
            double t = 0;
            for (int i = 0; i < 2; i++)
                t += an[i] * dtn[i];
            return t;
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

        // Derivatives of shape functions
        // with respect to global coordinates x.
        // xi - local coordinates;
        // xy[2][2] - nodal coordinates;
        // dnxy[2][1] - derivatives of shape functions (out);
        // returns  determinant of the Jacobian matrrix
        public static double deriv(double xi, double[][] xy, double[][] dnxy)
        {
            // Derivatives in local coords dN/dXi
            double[][] dnxe = new double[1][] { new double[2] };

            // End nodes
            dnxe[0][0] = -0.5;
            dnxe[0][1] = 0.5;

            // Jacobian matrix
            L = Math.Sqrt(Math.Pow((xy[1][0] - xy[0][0]), 2.0) + Math.Pow((xy[1][1] - xy[0][1]), 2.0));            
            double aj = L / 2.0;
            double det = aj;
            // Zero or negative determinant
            if (det <= 0)
                UTIL.errorMsg("Negative/zero Jacobian determinant for 2N element " + (float)det);

            // Jacobian inverse
            double aj00 = 1.0 / det;

            // Derivatives in global coordinates dN/dx, dN/dy
            dnxy[0][0] = aj00 * dnxe[0][0];
            dnxy[0][1] = aj00 * dnxe[0][1];

            

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
