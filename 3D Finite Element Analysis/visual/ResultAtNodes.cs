using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.elem;
using FEA3D.model;

namespace FEA3D.visual
{
    // Result values at nodes of the finite element model.
    public class ResultAtNodes 
    {
        //private SurfaceGeometry sg;
        double sgfmin, sgfmax, sgdeltaf;
        private FeModel fem;
        private int[] multNod;
        private double[][] stressNod;
        private int[] sgsNodes;
        public double[] sgfun;

        private double sm, psi, si, f;
        static double THIRD = 1.0/3.0;
        static double SQ3 = Math.Sqrt(3.0);

        // Constructor for results at nodes.
        // sg - geometry of model surface.
        public ResultAtNodes( FeModel fem) 
        {

            //this.sg = sg;
            this.fem = fem;
            multNod = new int[fem.nNod];
            stressNod = new double[fem.nNod][];// 2*fem.nDim];
            for (int i = 0; i < fem.nNod; i++)
                stressNod[i] = new double[2 * fem.nDim];
            sgsNodes = new int[fem.nNod];
            sgfun = new double[fem.nNod];
            feStressAtNodes();
        }

        // FE stresses at nodes: global array stressNod.
        private void feStressAtNodes() {

            double[][] elStressInt = new double[8][]{ new double[6], new double[6],new double[6],new double[6],
                                                      new double[6],new double[6],new double[6],new double[6]};
            double[][] elStressNod = new double[20][]{ new double[6],new double[6],new double[6],new double[6],new double[6],
                                                       new double[6],new double[6],new double[6],new double[6],new double[6],
                                                       new double[6],new double[6],new double[6],new double[6],new double[6],
                                                       new double[6],new double[6],new double[6],new double[6],new double[6]};

            for (int i = 0; i < fem.nNod; i++) {
                multNod[i] = 0;
                for (int j = 0; j < 2*fem.nDim; j++)
                    stressNod[i][j] = 0;
            }

            for (int iel = 0; iel < fem.nEl; iel++) {
                Element el = fem.elems[iel];
                for (int ip = 0; ip < el.str.Length; ip++)
                    for (int j = 0; j < 2*fem.nDim; j++)
                        elStressInt[ip][j] = el.str[ip].sStress[j];

                el.setElemXy();
                el.extrapolateToNodes(elStressInt, elStressNod);

                // Assemble stresses
                for (int i=0; i<fem.elems[iel].ind.Length; i++) 
                {
                    int jind = fem.elems[iel].ind[i] - 1;
                    if (jind >= 0) {
                        for (int k = 0; k < 2*fem.nDim; k++)
                           stressNod[jind][k] += elStressNod[i][k];
                        multNod[jind] += 1;
                    }
                }
            }
            // Divide by node multiplicity factor
            for (int i = 0; i < fem.nNod; i++) {
                for (int j = 0; j < 2*fem.nDim; j++)
                    stressNod[i][j] /= multNod[i];
            }
        }

        // Set array sg.fun[] containing requested value at nodes.
        // parm - requested result value.
        // displ - displacement vector.
        public void setParmAtNodes(VisData.parms parm, double[] displ) {

            sgfmin = 1.0e77;
            sgfmax = -1.0e77;

            for (int node = 0; node < fem.nNod; node++) 
            {
                //if (sgsNodes[node] >= 0) 
                {
                    if (parm == VisData.parms.s1 || parm == VisData.parms.s2 || parm == VisData.parms.s3 || parm == VisData.parms.si || parm == VisData.parms.s13)
                        setEquivalentStress(node);

                    switch (parm) 
                    {
                        case VisData.parms.ux:
                            f = displ[node*fem.nDim];          
                            break;
                        case VisData.parms.uy:
                            f = displ[node*fem.nDim + 1];      
                            break;
                        case VisData.parms.uz:
                            if (fem.nDim == 3)
                                f = displ[node*fem.nDim + 2];
                            else f = 0;
                            break;
                        case VisData.parms.sx:
                            f = stressNod[node][0];           
                            break;
                        case VisData.parms.sy:
                            f = stressNod[node][1];   
                            break;
                        case VisData.parms.sz:
                            f = stressNod[node][2];    
                            break;
                        case VisData.parms.sxy:
                            f = stressNod[node][3];      
                            break;
                        case VisData.parms.syz:
                            f = stressNod[node][4];   
                            break;
                        case VisData.parms.szx:
                            f = stressNod[node][5];      
                            break;
                        case VisData.parms.s1:
                            f = sm + 2*THIRD*si*Math.Cos(psi);                            
                            break;
                        case VisData.parms.s2:
                            f = sm - 2*THIRD*si * Math.Cos(THIRD*Math.PI + psi);   
                            break;
                        case VisData.parms.s3:
                            f = sm - 2*THIRD*si * Math.Cos(THIRD*Math.PI-psi);  
                            break;
                        case VisData.parms.si:
                            f = si;                    
                            break;
                        case VisData.parms.s13:
                            f = THIRD*si*(Math.Cos(psi) + Math.Cos(THIRD*Math.PI-psi));
                            break;
                        default:
                            break;
                    }
                    sgfun[node] = f;
                    sgfmin = Math.Min(sgfmin, f);
                    sgfmax = Math.Max(sgfmax, f);
                }
            }
            //if (!(VisData.fMin == 0.0 && VisData.fMax == 0.0)) 
            //{
            //    sg.fmax = VisData.fMax;
            //    sg.fmin = VisData.fMin;
            //}
            if (sgfmax - sgfmin < 1.0e-6) 
                sgfmax += 1.0e-6;
            sgdeltaf = sgfmax - sgfmin;
        }

        // Compute stress invariants and equivalent stress.
        private void setEquivalentStress(int node) 
        {
            // Stresses
            double sx  = stressNod[node][0];
            double sy = stressNod[node][1];
            double sz = stressNod[node][2];
            double sxy = stressNod[node][3];
            double syz, szx;
            if (fem.nDim == 3) {
                syz = stressNod[node][4];
                szx = stressNod[node][5];
            }
            else { syz = 0;  szx = 0; }
            // Mean stress
            sm = THIRD*(sx + sy + sz);
            // Deiatoric stresses
            double dx = sx - sm;
            double dy = sy - sm;
            double dz = sz - sm;
            // Second and third deviatoric invariants
            double J2 =  0.5*(dx*dx + dy*dy + dz*dz)
                      + sxy*sxy + syz*syz + szx*szx;
            double J3 = dx*dy*dz + 2*sxy*syz*szx
                      - dx*syz*syz - dy*szx*szx - dz*sxy*sxy;
            // Angle
            double temp = 1.5 * SQ3 * J3 / Math.Sqrt(J2 * J2 * J2);
            psi = THIRD*Math.Acos(Math.Min(Math.Max(temp,-1.0),1.0));
            if (psi.ToString().ToLower().Equals("nan"))
                Console.WriteLine(psi + "  " + 1.5 * SQ3 * J3 + "  " + Math.Sqrt(J2 * J2 * J2) + "  " + "  " + THIRD * Math.Acos(temp));
            // Equivalent stress
            si = Math.Sqrt(3*J2);
        }

    }
}
