using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.model;
using FEA3D.fea;
using FEA3D.util;
using FEA3D.elem;

namespace FEA3D.gener
{
    // Generate mesh of quadratic elements inside a ring secton.
    // Input: nx, ny - number of elements along circumference and thickness;
    //  ri, ro, theta1, theta2 - inner radius, outer radius, start angle, end angle;
    //  [mat] - material name.
    public class ring4
    {

        private FeModel m;
        enum vars
        {
            nc, nt, ri, ro, theta1, theta2, mat, thick, end
        }

        private vars name;

        private int nc, nt;
        String mat = "1";
        //private double[] cs;
        //private double[] ts;

        private double ri, ro, theta1, theta2;
        private double theta, dtheta, dt, t;

        public ring4()
        {
            String modelName = CSMGEN.RD.next();
            CSMGEN.PR.WriteLine("Ring: {0}", modelName);
            readData();
            printData();
            m = new FeModel(CSMGEN.RD, CSMGEN.PR);
            generateMesh();
            //CSMGEN.blocks.put(modelName, m);
            CSMGEN.blocks.Add(modelName, m);

            CSMGEN.PR.WriteLine("Mesh " + modelName + ": nEl = {0}  nNod = {1}", m.nEl, m.nNod);
        }

        private void readData()
        {
            String varName, varname;

            while (CSMGEN.RD.hasNext())
            {
                varName = CSMGEN.RD.next();
                varname = varName.ToLower();
                if (varName.Equals("#"))
                {
                    CSMGEN.RD.nextLine(); continue;
                }
                try
                {
                    //name = vars.valueOf(varname);
                    name = (vars)System.Enum.Parse(typeof(vars), varname);
                }
                catch (Exception e)
                {
                    UTIL.errorMsg("Variable name is not found: " + varName);
                }

                switch (name)
                {
                    case vars.nc:
                        nc = CSMGEN.RD.readInt();
                        break;
                    case vars.nt:
                        nt = CSMGEN.RD.readInt();
                        break;
                    case vars.ri:
                        ri = CSMGEN.RD.readDouble();
                        //xs = new double[nx + 1];
                        //for (int i = 0; i <= nx; i++)
                        //    xs[i] = CSMGEN.RD.readDouble();
                        break;
                    case vars.ro:
                        ro = CSMGEN.RD.readDouble();
                        //ys = new double[ny + 1];
                        //for (int i = 0; i <= ny; i++)
                        //    ys[i] = CSMGEN.RD.readDouble();
                        break;
                    case vars.theta1:
                        theta1 = CSMGEN.RD.readDouble();
                        break;
                    case vars.theta2:
                        theta2 = CSMGEN.RD.readDouble();
                        break;
                    case vars.mat:
                        mat = CSMGEN.RD.next();
                        break;
                    case vars.thick:
                        t = CSMGEN.RD.readDouble();
                        break;
                    case vars.end:
                        return;
                }
            }
        }

        private void printData()
        {
            theta = theta2 - theta1;
            dtheta = theta / (/*2 */ nc);
            dt = (ro - ri) / (/*2 */ nt);
            
            CSMGEN.PR.WriteLine(" nc = {0,5}", nc);
            CSMGEN.PR.WriteLine(" nt = {0,5}", nt);
            CSMGEN.PR.Write(" ri = {0,5}  ro = {1,5}  theta1 = {2,5}  theta2 = {3,5}", ri, ro, theta1, theta2);
            CSMGEN.PR.Write(" cs:  ");
            //for (int i = 0; i <= nc; i++)
            //{
            //    cs[i] = ri * Math.Sin(dtheta * i * Math.PI / 180.0);
            //    CSMGEN.PR.Write("{0,7:0.000}", cs[i]);//%7.3f
            //}
            //CSMGEN.PR.Write("\n ts:  ");
            //for (int i = 0; i <= nt; i++)
            //{
            //    ts[i] = ri + dt * i;
            //    CSMGEN.PR.Write("{0,7:0.000}", ts[i]);
            //}
            CSMGEN.PR.WriteLine();
        }

        private void generateMesh()
        {
            
            int[] ind = new int[4];
            m.nDim = 2;

            // Connectivity array
            int nh = nt;
            int nv = nc;
            m.nEl = nh * nv;
            m.elems = new Element[m.nEl];
            int n = 0;

            for (int iv = 0; iv < nv; iv++)
            {
                for (int ih = 0; ih < nh; ih++)
                {
                    m.elems[n] = Element.newElement("quad4");
                    int in0 = iv * (nh + 1) + ih + 1;
                    ind[0] = in0;
                    ind[1] = in0 + 1;
                    int in1 = (iv + 1) * (nh + 1) + ih + 1;
                    ind[2] = in1 + 1;
                    ind[3] = in1;
                    m.elems[n].setElemConnectivities(ind);
                    m.elems[n].setElemMaterial(mat);
                    n++;
                }
            }

            // Node coordinate array
            m.nNod = (nt + 1) * (nc + 1);
            m.newCoordArray();
            n = 0;
            for (int iv = 0; iv < nv + 1; iv++)
            {
                for (int ih = 0; ih < nh + 1; ih++)
                {
                    //int pc = (ic + 1) / 2;
                    //if (ic % 2 == 0 && it % 2 == 0) // corner noder - 0, 2, 4, 6
                    //{
                        m.setNodeCoord(n, 0, (ri + dt * ih) * Math.Cos((theta1 + dtheta * iv) * Math.PI / 180.0));
                        m.setNodeCoord(n, 1, (ri + dt * ih) * Math.Sin((theta1 + dtheta * iv) * Math.PI / 180.0));
                        n++;
                    //}
                    //else if (ic % 2 == 1 && it % 2 == 0) // midside node on horizontal element side - 3, 7
                    //{
                        //m.setNodeCoord(n, 0, (ri + dt * it) * Math.Cos((theta1 + dtheta * ic) * Math.PI / 180.0));
                        //m.setNodeCoord(n, 1, (ri + dt * it) * Math.Sin((theta1 + dtheta * ic) * Math.PI / 180.0));
                        //n++;
                    //}
                    //else if (ic % 2 == 0 && it % 2 == 1) // midside node on vertical element side - 1, 5
                    //{
                        //m.setNodeCoord(n, 0, (ri + dt * it) * Math.Cos((theta1 + dtheta * ic) * Math.PI / 180.0));
                        //m.setNodeCoord(n, 1, (ri + dt * it) * Math.Sin((theta1 + dtheta * ic) * Math.PI / 180.0));
                        //n++;
                    //}
                }
            }
        }

    }
}