using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.model;
using FEA3D.fea;
using FEA3D.util;

namespace FEA3D.gener
{
    // Read mesh data from text file.
    // Input: modelName - name of the finite element model;
    // fileName - name of the file.
    public class readmesh
    {

        public readmesh()
        {

            String modelName = CSMGEN.RD.next();
            String fileName = CSMGEN.RD.next();
            CSMGEN.PR.WriteLine("ReadMesh:  {0}    {1}\n", modelName, fileName);
            FeScanner RD = new FeScanner(fileName);
            FeModel m = new FeModel(RD, CSMGEN.PR);
            m.readData();
            CSMGEN.blocks.Add(modelName, m);
            CSMGEN.PR.WriteLine("Mesh " + modelName + ": nEl = {0,9}  nNod = {1,9}\n", m.nEl, m.nNod);
        }

    }

}
