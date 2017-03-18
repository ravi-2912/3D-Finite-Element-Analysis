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
    // Generate mesh of linear elements inside a rectangle.
    // Input: nx, ny - number of elements along x and y;
    //  xs, ys - locations of element boundaries on x and y;
    //  [mat] - material name.
    public class rect4
    {

        private FeModel m;
        enum vars
        {
            nx, ny, xs, ys, mat, thick, end
        }

        private vars name;

        private int nx, ny;
        String mat = "1";
        private double[] xs;
        private double[] ys;
        private double t;

        public rect4()
        {
            String modelName = CSMGEN.RD.next();
            CSMGEN.PR.WriteLine("Rectangle4: {0}", modelName);
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
                    case vars.nx:
                        nx = CSMGEN.RD.readInt();
                        break;
                    case vars.ny:
                        ny = CSMGEN.RD.readInt();
                        break;
                    case vars.xs:
                        xs = new double[nx + 1];
                        for (int i = 0; i <= nx; i++)
                            xs[i] = CSMGEN.RD.readDouble();
                        break;
                    case vars.ys:
                        ys = new double[ny + 1];
                        for (int i = 0; i <= ny; i++)
                            ys[i] = CSMGEN.RD.readDouble();
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
            CSMGEN.PR.WriteLine(" nx = {0,5}", nx);
            CSMGEN.PR.WriteLine(" ny = {0,5}", ny);
            CSMGEN.PR.Write(" xs:  ");
            for (int i = 0; i <= nx; i++)
                CSMGEN.PR.Write("{0,7:0.000}", xs[i]);//%7.3f
            CSMGEN.PR.Write("\n ys:  ");
            for (int i = 0; i <= ny; i++)
                CSMGEN.PR.Write("{0,7:0.000}", ys[i]);
            CSMGEN.PR.WriteLine();
        }

        private void generateMesh()
        {
            int[] ind = new int[4];
            m.nDim = 2;

            // Connectivity array
            m.nEl = nx * ny;
            m.elems = new Element[m.nEl];

            int el = 0;
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    m.elems[el] = Element.newElement("quad4");
                    m.elems[el].t = this.t;
                    int in0 = iy * (nx + 1) + ix + 1;
                    ind[0] = in0;
                    ind[1] = in0 + 1;
                    int in1 = (iy + 1) * (nx + 1) + ix + 1;
                    ind[2] = in1 + 1;
                    ind[3] = in1;                    
                    m.elems[el].setElemConnectivities(ind);
                    m.elems[el].setElemMaterial(mat);
                    el++;
                }
            }

            // Node coordinate array
            m.nNod = (nx + 1) * (ny + 1);// (3 * nx + 2) * ny + 2 * nx + 1;
            m.newCoordArray();
            int n = 0;
            for (int iy = 0; iy < ny + 1; iy++)
                for (int ix = 0; ix < nx + 1; ix++)
                {
                    m.setNodeCoord(n, 0, xs[ix]);
                    m.setNodeCoord(n, 1, ys[iy]);
                    n++;                    
                }
            
        }

    }
}
