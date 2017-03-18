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
    // Generate mesh of quadratic elements inside a rectangle.
    // Input: nx, ny - number of elements along x and y;
    //  xs, ys - locations of element boundaries on x and y;
    //  [mat] - material name.
    public class rect8 {

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

        public rect8() {
            String modelName = CSMGEN.RD.next();
            CSMGEN.PR.WriteLine("Rectangle8: {0}", modelName);
            readData();
            printData();
            m = new FeModel(CSMGEN.RD, CSMGEN.PR);
            generateMesh();
            //CSMGEN.blocks.put(modelName, m);
            CSMGEN.blocks.Add(modelName, m);

            CSMGEN.PR.WriteLine("Mesh " + modelName + ": nEl = {0}  nNod = {1}", m.nEl, m.nNod);
        }

        private void readData() {
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
                        xs = new double[nx+1];
                        for (int i = 0; i <= nx; i++)
                            xs[i] = CSMGEN.RD.readDouble();
                        break;
                    case vars.ys:
                        ys = new double[ny+1];
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

        private void printData() {
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

        private void generateMesh() {
            int[] ind = new int[8];
            m.nDim = 2;

            // Connectivity array
            m.nEl = nx*ny;
            m.elems = new Element[m.nEl];

            int el = 0;
            for (int iy=0; iy<ny; iy++) {
                for (int ix=0; ix<nx; ix++) {
                    m.elems[el] = Element.newElement("quad8r");
                    int in0 = iy*(3*nx+2) + 2*ix;
                    ind[0] = in0 + 1;
                    ind[1] = in0 + 2;
                    ind[2] = in0 + 3;
                    int in1 = iy*(3*nx+2) + 2*nx + 1 + ix + 1;
                    ind[3] = in1 + 1;
                    ind[7] = in1;
                    int in2 = (iy+1)*(3*nx+2) + 2*ix;
                    ind[4] = in2 + 3;
                    ind[5] = in2 + 2;
                    ind[6] = in2 + 1;
                    m.elems[el].setElemConnectivities(ind);
                    m.elems[el].setElemMaterial(mat);
                    el++;
                }
            }

            // Node coordinate array
            m.nNod = (3 * nx + 2) * ny + 2 * nx + 1;
            m.newCoordArray();
            int n = 0;
            for (int iy=0; iy<2*ny+1; iy++) {
                int py = (iy+1)/2;
                for (int ix=0; ix<2*nx+1; ix++) {
                    int px = (ix+1)/2;
                    if (ix%2==0 && iy%2==0) {
                        m.setNodeCoord(n, 0, xs[px]);
                        m.setNodeCoord(n, 1, ys[py]);
                        n++;
                    }
                    else if (ix%2==1 && iy%2==0) {
                        m.setNodeCoord(n,0,0.5*(xs[px-1]+xs[px]));
                        m.setNodeCoord(n, 1, ys[py]);
                        n++;
                    }
                    else if (ix%2==0 && iy%2==1) {
                        m.setNodeCoord(n, 0, xs[px]);
                        m.setNodeCoord(n,1,0.5*(ys[py-1]+ys[py]));
                        n++;
                    }
                }
            }
        }

    }

}
