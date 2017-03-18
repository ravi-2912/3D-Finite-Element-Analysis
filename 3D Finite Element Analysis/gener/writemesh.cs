using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using FEA3D.model;
using FEA3D.fea;
using FEA3D.util;

namespace FEA3D.gener
{
    // Write mesh to file.
    // Input:  modelName - name of the finite element model;
    //  fileName - name of the file.
    public class writemesh
    {

        FeModel m;

        public writemesh()
        {

            String modelName = CSMGEN.RD.next();
            String fileName = CSMGEN.RD.next();
            CSMGEN.PR.Write("WriteMesh: {0}    {1}\n", modelName, fileName);

            StreamWriter WR = new FePrintWriter().getPrinter(fileName);

            if (CSMGEN.blocks.ContainsKey(modelName))
                m = (FeModel)CSMGEN.blocks[modelName];
            else UTIL.errorMsg("No such mesh block: " + modelName);

            WR.Write("# Model name: {0}\n", modelName);
            WR.Write("nNod = {0,5}\n", m.nNod);
            WR.Write("nEl = {0,5}\n", m.nEl);
            WR.Write("nDim = {0,5}\n", m.nDim);

            WR.Write("nodCoord\n");
            for (int i = 0; i < m.nNod; i++)
            {
                for (int j = 0; j < m.nDim; j++)
                {
                    WR.Write("{0,20:0.000000000}", m.getNodeCoord(i, j));
                }
                WR.WriteLine("");
            }

            WR.Write("\nelCon");
            for (int iel = 0; iel < m.nEl; iel++)
            {
                WR.Write("\n{0} {1,6}", m.elems[iel].name, m.elems[iel].matName);
                if (m.elems[iel].name.Equals("quad4") || m.elems[iel].name.Equals("quad8r") )//|| m.elems[iel].name.Equals("genquad4") || m.elems[iel].name.Equals("genquad8"))
                    WR.Write("{0,6}", m.elems[iel].t);
                if (m.elems[iel].name.Equals("truss22") || m.elems[iel].name.Equals("truss32"))//|| m.elems[iel].name.Equals("genquad4") || m.elems[iel].name.Equals("genquad8"))
                {
                    WR.Write("{0,12}", m.elems[iel].A);
                }
                int nind = m.elems[iel].ind.Length;
                for (int i = 0; i < nind; i++)
                    WR.Write("{0,12}", m.elems[iel].ind[i]);
            }
            WR.Write("\n\nend\n");
            WR.Close();
            CSMGEN.PR.Write("Mesh " + modelName + ": nEl = {0}  nNod = {1}\n", m.nEl, m.nNod);
        }

    }
}
