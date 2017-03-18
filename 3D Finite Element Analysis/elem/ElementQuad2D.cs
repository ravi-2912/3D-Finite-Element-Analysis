using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.util;
using FEA3D.model;
using FEA3D.material;

namespace FEA3D.elem
{
    // 2D quadratic isoparametric element (4-8 nodes)
    class ElementQuad2D : Element 
    {
        // Element edges (local numbers)
        private static int[][] faceInd = new int[][]{new int[]{0,1,2},new int[]{2,3,4},new int[]{4,5,6},new int[]{6,7,0}};
        // Shape functions
        private static double[] an = new double[8];
        // Derivatives of shape functions
        private static double[][] dnxy = new double[8][]{new double[2], new double [2], new double [2], new double [2], new double [2], new double [2], new double [2], new double [2]};
        // Displacements differentiation matrix
        private static double[][] bmat = new double[4][]{new double[16], new double[16], new double[16], new double[16]};
        // Elasticity matrix
        private static double[][] emat = new double[4][]{new double[4],new double[4],new double[4],new double[4]};
        // Thermal strains
        private static double[] ept = new double[4];
        // Radius in the axisymmetric problem
        private static double r;
        // Gauss rules for stiffness matrix, thermal vector,
        // surface load and stress integration
        private static GaussRule gk = new GaussRule(2,2); // from 2 to 3 - making full integation
        private static GaussRule gh = new GaussRule(3,2);
        private static GaussRule gf = new GaussRule(3,1);
        private static GaussRule gs = new GaussRule(2,2); // from 2 to 3 - making full integration

        // Constructor for 2D quadratic element
        public ElementQuad2D() : base ("quad8r", 8, 4) // from 4 to 9 making it full integration
        {
            //super ("quad8", 8, 4);
        }

        // Compute stiffness matrix
        public override void stiffnessMatrix() 
        {
            // Zeros to stiffness matrix kmat
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 16; j++)
                    kmat[i][j] = 0.0;

            // ld = length of strain/stress vector (3 or 4)
            int ld = (FeModel.stressState == FeModel.StrStates.axisym) ? 4 : 3;
            // Material mat
            mat = (Material)fem.materials[matName];
            if (mat == null) 
                UTIL.errorMsg("Element material name: " + matName);
            mat.elasticityMatrix(emat);

            // Gauss integration loop
            for (int ip = 0; ip < gk.nIntPoints; ip++)
            {
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gk.xii[ip], gk.eti[ip]);
                double dv = det * gk.wi[ip];
                if (FeModel.stressState == FeModel.StrStates.axisym)
                    dv *= 2.0 * Math.PI * r;
                // Upper symmetrical part of the stiffness matrix
                for (int i = 0; i < 16; i++) {
                    for (int j = i; j < 16; j++) {
                        double s = 0.0;
                        for (int k = 0; k < ld; k++) {
                            for (int l = 0; l < ld; l++) {
                                s += bmat[l][i] * emat[l][k] * bmat[k][j] * this.t;
                            }
                        }
                        kmat[i][j] += s * dv;
                    }
                }
            }
        }

        // Set displacement differentiation matrix bmat.
        // xi, et - local coordinates,
        // returns  determinant of Jacobian matrix
        private double setBmatrix(double xi, double et) 
        {
            // Derivatives of shape functions
            double det = ShapeQuad2D.deriv(xi, et, ind, xy, dnxy);
            if (det <= 0) 
                UTIL.errorMsg("Negative/zero 8Nr element area");
            if (FeModel.stressState == FeModel.StrStates.axisym) {
                ShapeQuad2D.shape(xi, et, ind, an);
                r = 0.0;
                for (int i = 0; i < 8; i++) r += an[i] * xy[i][0];
            }
            // Eight blocks of the displacement differentiation
            //   matrix
            for (int ib = 0; ib < 8; ib++) {
                bmat[0][2*ib]   = dnxy[ib][0];
                bmat[0][2*ib+1] = 0.0;
                bmat[1][2*ib]   = 0.0;
                bmat[1][2*ib+1] = dnxy[ib][1];
                bmat[2][2*ib]   = dnxy[ib][1];
                bmat[2][2*ib+1] = dnxy[ib][0];
                if(FeModel.stressState==FeModel.StrStates.axisym) {
                    bmat[3][2*ib]   = an[ib] / r;
                    bmat[3][2*ib+1] = 0.0;
                }
            }
            return det;
        }

        // Compute thermal vector
        public override void thermalVector() 
        {
            // Zeros to thermal vector evec
            for (int i = 0; i < 16; i++) evec[i] = 0.0;
            int ld = (FeModel.stressState
                        == FeModel.StrStates.axisym) ? 4 : 3;
            // Material mat
            mat = (Material)fem.materials[matName];
            mat.elasticityMatrix(emat);
            double alpha = mat.getAlpha();
            double nu = mat.getNu();

            // Gauss integration loop
            for (int ip = 0; ip < gh.nIntPoints; ip++) {
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gh.xii[ip], gh.eti[ip]);
                // Shape functions an
                ShapeQuad2D.shape(gh.xii[ip], gh.eti[ip], ind, an);
                double t = 0.0;
                for (int i = 0; i < 8; i++) t += an[i]*dtn[i];
                double dv = det*gh.wi[ip];
                if (FeModel.stressState==FeModel.StrStates.axisym)
                    dv *= 2.0*Math.PI*r;
                ept[0] = alpha*t;
                if(FeModel.stressState==FeModel.StrStates.plstrain)
                        ept[0] *= (1 + nu);
                ept[1] = ept[0];
                ept[2] = 0.0;
                ept[3] = ept[0];

                for (int i = 0; i < 16; i++) {
                    double s = 0;
                    for (int j = 0; j < ld; j++) {
                        for (int k = 0; k < ld; k++) {
                            s += bmat[k][i] * emat[j][k] * ept[j] * this.t;
                        }
                    }
                    evec[i] += s*dv;
                }
            }
        }

        // Set nodal equivalent of distributed face load to evec.
        // surLd - object describing element face load;
        // returns loaded element face
        // or -1 (loaded nodes do not match elem face)
        public override int equivFaceLoad(ElemFaceLoad surLd) 
        {
            // Shape functons
            double[] an = new double[3];
            // Derivatives of shape functions
            double[] xin = new double[3];

            for (int i=0; i<16; i++) evec[i] = 0.0;
            int loadedFace = surLd.rearrange(faceInd, ind);
            if (loadedFace == -1) return -1;

            // Gauss integration loop
            for (int ip=0; ip<gf.nIntPoints; ip++){
                ShapeQuad2D.shapeDerivFace(gf.xii[ip], surLd.faceNodes[1], an, xin);
                double p = r = 0.0;
                double xs = 0.0;
                double ys = 0.0;
                for (int i=0; i<3; i++){
                    p  += an[i]*surLd.forceAtNodes[i];
                    int j = faceInd[loadedFace][i];
                    r  += an[i]*xy[j][0];
                    xs += xin[i]*xy[j][0];
                    ys += xin[i]*xy[j][1];
                }
                double dl = Math.Sqrt(xs*xs+ys*ys);
                double ds = dl;
                if (FeModel.stressState==FeModel.StrStates.axisym)
                     ds *= 2.0*Math.PI*r;
                double p1, p2;
                // direction=0 - normal load, =1,2 - along axes x,y
                if (surLd.direction == 0 && ds > 0.0) {
                    p1 = p*ys/dl;
                    p2 = -p*xs/dl;
                }
                else if (surLd.direction == 1) { p1 = p; p2 = 0; }
                else                           { p1 = 0; p2 = p; }

                for (int i=0; i<3; i++){
                    int j = faceInd[loadedFace][i];
                    evec[2*j    ]  += an[i]*p1*ds*gf.wi[ip]*this.t;
                    evec[2*j + 1]  += an[i]*p2*ds*gf.wi[ip]*this.t;
                }
            }
            return loadedFace;
        }

        // Compute equivalent stress vector (with negative sign)
        public override void equivStressVector() 
        {
            for (int i = 0; i < 16; i++) evec[i] = 0.0;
            int ld = (FeModel.stressState == FeModel.StrStates.axisym) ? 4 : 3;

            for (int ip = 0; ip < gs.nIntPoints; ip++) {
                // Accumulated stress
                double[] s = new double[4]; // changef 4 to 9 for full integration
                for (int i = 0; i < 4; i++) // changef i<4 to i<9 for full integration
                    s[i] = str[ip].sStress[i] + str[ip].dStress[i];
                // Set displacement differentiation matrix bmat
                double det = setBmatrix(gs.xii[ip], gs.eti[ip]);
                double dv = det*gs.wi[ip];
                if (FeModel.stressState == FeModel.StrStates.axisym)
                    dv *= 2.0*Math.PI*r;

                for (int i=0; i<16; i++) {
                    double a = 0;
                    for (int j=0; j<ld; j++) 
                        a += bmat[j][i]*s[j]*this.t;
                    evec[i] -= a*dv;
                }
            }
        }
        

        // Extrapolate values from integration points to nodes.
        // fip [4][4] - values at integration points;
        // fn [8][4] - values at nodes (out)
        public override void extrapolateToNodes(double[][] fip, double[][] fn) 
        {
            double A = 1.0 + 0.5*Math.Sqrt(3.0),
                   B = -0.5,
                   C = 1.0 - 0.5*Math.Sqrt(3.0);
            // Extrapolation matrix
            double[][] lim = new double[][]{new double[]{A, B, B, C},
                                            new double[]{B, C, A, B},
                                            new double[]{C, B, B, A},
                                            new double[]{B, A, C, B}};

            for (int i = 1; i < 8; i+=2)
                for (int j = 0; j < 4; j++) fn[i][j] = 0;

            for (int corner = 0; corner < 4; corner++) 
            {
                int n = (corner==0) ? 7 : 2*corner-1;
                for (int k = 0; k < 4; k++) {
                    double c = 0.0;
                    for (int ip = 0; ip < 4; ip++)
                    {
                        c += lim[corner][ip] * fip[ip][k];
                        //Console.WriteLine(fip[ip][3]);
                    }
                    fn[2*corner][k] = c;  // corner node
                    fn[n][k] += 0.5*c;
                    fn[2*corner+1][k] += 0.5*c;
                }
            }
            //for (int i = 0; i < 8; i++)
            //    Console.Write("{0,7:0.00}", fn[i][3]);
            //Console.WriteLine();
        }

        // Get local node numbers for element faces.
        // returns elementFaces[nFaces][nNodesOnFace]
        public override int[][] getElemFaces() {
            return faceInd;
        }

        // Get strains at integration point.
        // ip - integration point number (stress);
        // returns  strain vector (ex, ey, gxy, ez)
        public override double[] getStrainsAtIntPoint(int ip) {
            // Set displacement differentiation matrix bmat
            setBmatrix(gs.xii[ip], gs.eti[ip]);
            double[] strain = new double[4];
            for (int i=0; i<4; i++) {
                strain[i] = 0;
                for (int j=0; j<16; j++)
                    strain[i] += bmat[i][j]*evec[j];
            }
            return strain;
        }

        // Get temperature at integration point (stress)
        public override double getTemperatureAtIntPoint(int ip) {
            ShapeQuad2D.shape(gs.xii[ip], gs.eti[ip], ind, an);
            double t = 0;
            for (int i=0; i<8; i++) t += an[i]*dtn[i];
            return t;
        }
    }
}
