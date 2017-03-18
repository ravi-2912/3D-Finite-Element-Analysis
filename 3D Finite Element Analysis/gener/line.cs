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
    class line
    {
        private FeModel m;
        enum vars
        {
            elem, n, xs, ys, zs, mat, ndim, area, end
        }

        private vars name;

        String elType;
        private int n;
        String mat = "1";
        private double[] xs;
        private double[] ys;
        private double[] zs;
        private double a = 1.0;
        private int ndim = 2;

        public line()
        {
            String modelName = CSMGEN.RD.next();
            CSMGEN.PR.WriteLine("line: {0}", modelName);
            
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
                    case vars.elem:
                        elType = CSMGEN.RD.next();
                        break;
                    case vars.n:
                        n = CSMGEN.RD.readInt();
                        break;
                    case vars.xs:
                        xs = new double[2];
                        for (int i = 0; i <2; i++)
                            xs[i] = CSMGEN.RD.readDouble();
                        break;
                    case vars.ys:
                        ys = new double[2];
                        for (int i = 0; i <2; i++)
                            ys[i] = CSMGEN.RD.readDouble();
                        break;
                    case vars.zs:
                        zs = new double[2];
                        for (int i = 0; i < 2; i++)
                            zs[i] = CSMGEN.RD.readDouble();
                        break;
                    case vars.mat:
                        mat = CSMGEN.RD.next();
                        break;
                    case vars.area:
                        a = CSMGEN.RD.readDouble();
                        break;
                    case vars.ndim:
                        ndim = CSMGEN.RD.readInt();
                        break;
                    case vars.end:
                        return;
                }
            }
        }

        private void printData()
        {
            CSMGEN.PR.Write("Element Type: " + elType);
            CSMGEN.PR.WriteLine(" n = {0,5}", n);
            CSMGEN.PR.WriteLine(" Start Coordinate: ({0,7:0.000}, {1,7:0.000})", xs[0], ys[0]);
            CSMGEN.PR.WriteLine(" End   Coordinate: ({0,7:0.000}, {1,7:0.000})", xs[1], ys[1]);
            CSMGEN.PR.WriteLine(" Material: " + mat);
            CSMGEN.PR.WriteLine(" Cross Section Area: {0,7:0.000}", a);
            CSMGEN.PR.WriteLine();
        }

        private void generateMesh()
        {
            int[] ind = new int[2];
            m.nDim = this.ndim;

            // Connectivity array
            m.nEl = this.n;
            m.elems = new Element[m.nEl];

            int el = 0;
            for (int i = 0; i < this.n; i ++)
            {
                m.elems[el] = Element.newElement(this.elType);
                m.elems[el].A = this.a;
                ind[0] = i + 1;
                ind[1] = i + 2;
                m.elems[el].setElemConnectivities(ind);
                m.elems[el].setElemMaterial(mat);
                el++;
            }
                

            // Node coordinate array
            m.nNod = this.n + 1;
            m.newCoordArray();
            int nid = 0;
            double L, l;
            double lij = 0, mij = 0, nij = 0;
            L = Math.Sqrt(Math.Pow(ys[1] - ys[0], 2) + Math.Pow(xs[1] - xs[0], 2) + (m.nDim == 3 ? Math.Pow(zs[1] - zs[0], 2) : 0.0));
            l = L / (double)this.n;
            lij = (xs[1] - xs[0]) / L;
            mij = (ys[1] - ys[0]) / L;
            nij = (m.nDim == 3 ? (zs[1] - zs[0]) / L : 0.0);

            for (int i = 0; i < this.n + 1; i++ )
            {
                m.setNodeCoord(nid, 0, xs[0] + i * l * lij);
                m.setNodeCoord(nid, 1, ys[0] + i * l * mij);
                if(m.nDim==3)
                    m.setNodeCoord(nid, 2, zs[0] + i * l * nij);
                nid++;
            }
            
        }
    }
}
