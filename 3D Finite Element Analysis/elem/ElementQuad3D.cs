using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.model;
using FEA3D.material;
using FEA3D.util;

namespace FEA3D.elem
{
    // 3D 8-20 node isoparametric brick-type element.
    public class ElementQuad3D : Element 
    {

        private static int[][] faceInd = new int[6][]{
            new int[8]{ 0, 8,12,19,18,11, 6, 7},
            new int[8]{ 2, 3, 4,10,16,15,14, 9},
            new int[8]{ 0, 1, 2, 9,14,13,12, 8},
            new int[8]{ 4, 5, 6,11,18,17,16,10},
            new int[8]{ 0, 7, 6, 5, 4, 3, 2, 1},
            new int[8]{12,13,14,15,16,17,18,19} };
        private static double[] an = new double[20];
        private static double[][] dnxy = new double[20][]{new double[3], new double[3], new double[3], new double[3], new double[3], 
                                                          new double[3], new double[3], new double[3], new double[3], new double[3], 
                                                          new double[3], new double[3], new double[3], new double[3], new double[3], 
                                                          new double[3], new double[3], new double[3], new double[3], new double[3]};
        // Gauss rules for stiffness matrix, thermal vector,
        // surface load and stress integration
        private static GaussRule gk = new GaussRule(14,3);
        private static GaussRule gh = new GaussRule(3,3);
        private static GaussRule gf = new GaussRule(3,2);
        private static GaussRule gs = new GaussRule(2,3);

        // Constructor for 3D 20 node element.
        public ElementQuad3D() : base("hex20r", 20, 8)
        {
            //super ("hex20", 20, 8);
        }

        // Compute stiffness matrix
        public override void stiffnessMatrix() {

            for (int i = 0; i < 60; i++)
                for (int j = i; j < 60; j++) kmat[i][j] = 0.0;
            // Material mat
            mat = (Material)fem.materials[matName];
            if (mat == null) 
                UTIL.errorMsg("Element material name: " + matName);
            double lambda = mat.getLambda();
            double mu     = mat.getMu();
            double beta   = lambda + 2*mu;

            for (int ip = 0; ip < gk.nIntPoints; ip++) {
                double det = ShapeQuad3D.deriv(gk.xii[ip], gk.eti[ip], gk.zei[ip], ind, xy, dnxy);
                double dv = det*gk.wi[ip];
                // Upper symmetrical part of the matrix by rows
                for (int i = 0; i < 20; i++) { // i = row
                    // dNi/dx, dNi/dy, dNi/dz
                    double dix = dnxy[i][0];
                    double diy = dnxy[i][1];
                    double diz = dnxy[i][2];
                    for (int j = i; j <20; j++) { // j = column
                        // dNj/dx, dNj/dy, dNj/dz
                        double djx = dnxy[j][0];
                        double djy = dnxy[j][1];
                        double djz = dnxy[j][2];

                        kmat[i*3  ][j*3  ] += (beta*dix*djx
                            + mu*(diy*djy + diz*djz))*dv;
                        kmat[i*3  ][j*3+1] += (lambda*dix*djy
                            + mu*diy*djx)*dv;
                        kmat[i*3  ][j*3+2] += (lambda*dix*djz
                            + mu*diz*djx)*dv;

                        if (j > i) kmat[i*3+1][j*3  ]
                            += (lambda*diy*djx + mu*dix*djy)*dv;
                        kmat[i*3+1][j*3+1] += (beta*diy*djy
                            + mu*(diz*djz + dix*djx))*dv;
                        kmat[i*3+1][j*3+2] += (lambda*diy*djz
                            + mu*diz*djy)*dv;

                        if (j > i) {
                            kmat[i*3+2][j*3  ]
                              += (lambda*diz*djx + mu*dix*djz)*dv;
                            kmat[i*3+2][j*3+1]
                              += (lambda*diz*djy + mu*diy*djz)*dv;
                        }
                        kmat[i*3+2][j*3+2] += (beta*diz*djz
                                + mu*(dix*djx + diy*djy))*dv;
                    }
                }
            }
        }

        // Compute thermal vector
        public override void thermalVector() 
        {

            for (int i = 0; i < 60; i++) evec[i] = 0.0;

            mat = (Material)fem.materials[matName];
            double alpha = mat.getAlpha();
            double lambda = mat.getLambda();
            double mu = mat.getMu();
            double g = 3.0*lambda + 2.0*mu;

            for (int ip = 0; ip < gh.nIntPoints; ip++) {
                ShapeQuad3D.shape(gh.xii[ip], gh.eti[ip], gh.zei[ip], ind, an);
                // Temperature at integration point
                double t = 0;
                for (int i = 0; i < 20; i++) t += an[i]*dtn[i];
                double det = ShapeQuad3D.deriv(gh.xii[ip], gh.eti[ip], gh.zei[ip], ind, xy, dnxy);
                double dv = g*alpha*t*det*gh.wi[ip];
                for (int i=0; i<20; i++) {
                    for (int j=0; j<3; j++) {
                        evec[i*3+j] += dnxy[i][j]*dv;
                    }
                }
            }
        }

        // Set nodal equivalent of distributed face load to evec.
        // surLd - object describing element face load;
        // returns loaded element face
        // or -1 (loaded nodes does not match element face)
        public override int equivFaceLoad(ElemFaceLoad surLd) 
        {
            // Shape functons
            double[] an = new double[8];
            // Derivatives of shape functions
            double[][] xin = new double[8][]{new double[2], new double[2], new double[2], new double[2], 
                                              new double[2], new double[2], new double[2], new double[2]};
            // Tangent vectors along xi and eta
            double[][] e = new double[2][]{new double[3], new double[3]};
            // Normal vector
            double[] g = new double[3];
            double[] ps = new double[3];

            for (int i=0; i<60; i++) evec[i] = 0.0;

            int loadedFace = surLd.rearrange(faceInd, ind);
            if (loadedFace == -1) return -1;

            for (int ip=0; ip<gf.nIntPoints; ip++){
                ShapeQuad3D.shapeDerivFace(gf.xii[ip], gf.eti[ip], ind, an, xin);
                double p = 0.0;
                for (int i=0; i<8; i++)
                    p += an[i]*surLd.forceAtNodes[i];
                // Tangent vectors
                for (int i=0; i<2; i++) {
                    for (int j=0; j<3; j++) {
                        double s = 0;
                        for (int k=0; k<8; k++)
                            s += xin[k][i] * xy[faceInd[loadedFace][k]][j];
                        e[i][j] = s;
                     }
                 }
                 // Normal vector g
                 g[0] = (e[0][1]*e[1][2]-e[1][1]*e[0][2]);
                 g[1] = (e[0][2]*e[1][0]-e[1][2]*e[0][0]);
                 g[2] = (e[0][0]*e[1][1]-e[1][0]*e[0][1]);

                // Element of surface ds
                double ds = Math.Sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
                if (ds<=0) 
                    UTIL.errorMsg("Negative/zero element face");
                // Surface load components ps:
                // direction=0 - normal, x=1, y=2, z=3
                if (surLd.direction == 0) {
                    for (int i=0; i<3; i++) ps[i] = p*g[i]/ds;
                }
                else {
                    for (int i=0; i<3; i++) ps[i] = 0;
                    ps[surLd.direction-1] = p;
                }
                for (int i=0; i<8; i++) {
                    int k = faceInd[loadedFace][i];
                    for (int j=0; j<3; j++) {
                        evec[3*k + j] += an[i]*ps[j]*ds*gf.wi[ip];
                    }
                }
            }
            return loadedFace;
        }

        // Compute equivalent stress vector (with negative sign)
        public override void equivStressVector()
        {

            for (int i = 0; i < 60; i++) evec[i] = 0.0;

            for (int ip = 0; ip < gs.nIntPoints; ip++) {
                // Accumulated stress  s
                double[] s = new double[6];
                for (int i=0; i<6; i++)
                    s[i] = str[ip].sStress[i] + str[ip].dStress[i];
                double det = ShapeQuad3D.deriv(gs.xii[ip],gs.eti[ip], gs.zei[ip], ind, xy, dnxy);
                double dv = det * gs.wi[ip];

                for (int i = 0; i < 20; i++) {
                    double a0 = dnxy[i][0];
                    double a1 = dnxy[i][1];
                    double a2 = dnxy[i][2];
                    evec[i*3  ] -= (a0*s[0]+a1*s[3]+a2*s[5])*dv;
                    evec[i*3+1] -= (a1*s[1]+a0*s[3]+a2*s[4])*dv;
                    evec[i*3+2] -= (a2*s[2]+a1*s[4]+a0*s[5])*dv;
                }
            }
        }

        // Extrapolate stresses from integration points to nodes.
        // fip [8][6] - stresses at integration points;
        // fn [20][6] - stresses at nodes (out)
        public override void extrapolateToNodes(double[][] fip, double[][] fn)
        {
            // Vertices
            int[] vn = new int[8]{0, 2, 4, 6, 12, 14, 16, 18};
             // Midside nodes
            int[] mn = new int[8]{8, 9, 10, 11, 8, 9, 10, 11};
            // Extrapolation matrix
            double A =  0.25*(5.0 + 3.0*Math.Sqrt(3.0)),
                   B = -0.25*(Math.Sqrt(3.0) + 1.0),
                   C =  0.25*(Math.Sqrt(3.0) - 1.0),
                   D =  0.25*(5.0 - 3.0*Math.Sqrt(3.0));
            double[][] lim = new double[8][]{new double[8]{A, B, B, C, B, C, C, D},
                                             new double[8]{B, C, C, D, A, B, B, C},
                                             new double[8]{C, D, B, C, B, C, A, B},
                                             new double[8]{B, C, A, B, C, D, B, C},
                                             new double[8]{B, A, C, B, C, B, D, C},
                                             new double[8]{C, B, D, C, B, A, C, B},
                                             new double[8]{D, C, C, B, C, B, B, A},
                                             new double[8]{C, B, B, A, D, C, C, B}};

            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 6; j++) fn[i][j] = 0;

            for (int vertex = 0; vertex < 8; vertex++) {
                int i = vn[vertex];  // node at vertex
                int im = i - 1;
                if (i == 0) im = 7;
                if (i == 12) im = 19;
                for (int k = 0; k < 6; k++) {
                    double c = 0.0;
                    for (int j = 0 ; j < 8; j++)
                        c += fip[j][k]*lim[vertex][j];
                    fn[i][k] = c;
                    fn[im][k] += 0.5*c;
                    fn[i+1][k] += 0.5*c;
                    fn[mn[vertex]][k] += 0.5*c;
                }
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
        // returns  strain vector (ex, ey, ez, gxy, gyz, gzx)
        public override double[] getStrainsAtIntPoint(int ip)
        {

            // Derivatives of shape functions
            ShapeQuad3D.deriv(gs.xii[ip], gs.eti[ip], gs.zei[ip],ind, xy, dnxy);

            // Derivatives of displacements
            double dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz;
            dux=duy=duz=dvx=dvy=dvz=dwx=dwy=dwz = 0;
            for (int i = 0; i < 20; i++) {
                double dnx = dnxy[i][0];
                double dny = dnxy[i][1];
                double dnz = dnxy[i][2];
                double u = evec[3*i  ];
                double v = evec[3*i+1];
                double w = evec[3*i+2];
                dux += dnx*u;  duy += dny*u;  duz += dnz*u;
                dvx += dnx*v;  dvy += dny*v;  dvz += dnz*v;
                dwx += dnx*w;  dwy += dny*w;  dwz += dnz*w;
            }
            // Strains
            double[] strain = new double[6];
            strain[0] = dux;   strain[1] = dvy;   strain[2] = dwz;
            strain[3] = duy + dvx;
            strain[4] = dvz + dwy;
            strain[5] = duz + dwx;
            return strain;
        }

        // Returns temperature at integration point (stress)
        public override double getTemperatureAtIntPoint(int ip)
        {

            ShapeQuad3D.shape(gs.xii[ip], gs.eti[ip], gs.zei[ip],ind, an);
            double t = 0;
            for (int i=0; i<20; i++) t += an[i]*dtn[i];
            return t;
        }

    }
}
