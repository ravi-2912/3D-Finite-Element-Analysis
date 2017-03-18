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
    // Copy model.
    // Input: modelNameA - name of the model to be copied;
    //  modelNameB - name of the resulting model.
    public class copy
    {

        private FeModel mA;

        public copy()
        {

            String modelNameA = CSMGEN.RD.next();
            String modelNameB = CSMGEN.RD.next();
            CSMGEN.PR.WriteLine("Copy: {0} -> {1}\n", modelNameA, modelNameB);
            if (modelNameA.Equals(modelNameB)) return;
            if (CSMGEN.blocks.ContainsKey(modelNameA))
                mA = (FeModel)CSMGEN.blocks[modelNameA];
            else
                UTIL.errorMsg("No such mesh block: " + modelNameA);
            FeModel mB = copyMesh();
            CSMGEN.blocks.Add(modelNameB, mB);
            CSMGEN.PR.WriteLine("Mesh " + modelNameB + ": nEl = {0,9}  nNod = {1,9}\n", mB.nEl, mB.nNod);
        }

        private FeModel copyMesh()
        {

            FeModel mB = new FeModel(CSMGEN.RD, CSMGEN.PR);
            mB.nDim = mA.nDim;

            mB.nNod = mA.nNod;
            mB.newCoordArray();
            for (int i = 0; i < mB.nNod; i++)
            {
                for (int j = 0; j < mB.nDim; j++)
                    mB.setNodeCoords(i, mA.getNodeCoords(i));
            }

            mB.nEl = mA.nEl;
            mB.elems = new Element[mB.nEl];
            for (int el = 0; el < mB.nEl; el++)
            {
                mB.elems[el] = Element.newElement(mA.elems[el].name);
                mB.elems[el].setElemConnectivities(mA.elems[el].ind);
                mB.elems[el].matName = mA.elems[el].matName;
            }

            return mB;
        }

    }
}
