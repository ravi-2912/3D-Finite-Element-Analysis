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
    // Paste two meshes.
// Input: modelNameA - name of first mesh to be pasted;
//  modelNameB - name of second mesh to be pasted;
//  modelNameC - name of resulting mesh;
//  [eps] - coordinate tolerance for joining nodes.
public class connect {

    private FeModel mA, mB;
    private double eps = 0.0001;
    private int nConnected;
    private int[] newNodesB;

    public connect() 
    {

        String modelNameA = CSMGEN.RD.next();
        String modelNameB = CSMGEN.RD.next();
        String modelNameC = CSMGEN.RD.next();

        CSMGEN.PR.WriteLine("Connect: {0} + {1} -> {2}", modelNameA, modelNameB, modelNameC);
        
        readData();
        
        if (CSMGEN.blocks.ContainsKey(modelNameA))
            mA = (FeModel) CSMGEN.blocks[modelNameA];
        else UTIL.errorMsg("No such mesh block: " + modelNameA);
        if (CSMGEN.blocks.ContainsKey(modelNameB))
            mB = (FeModel) CSMGEN.blocks[modelNameB];
        else UTIL.errorMsg("No such mesh block: " + modelNameB);
        if (mA.nDim != mB.nDim) 
            UTIL.errorMsg("Models with different nDim");

        findCoincidentNodes(); 
        FeModel mC = pasteModels();
        
        CSMGEN.blocks.Add(modelNameC, mC);
        CSMGEN.PR.Write(" {0} node pairs connected\n",nConnected);
        CSMGEN.PR.Write("Mesh " + modelNameC + ": nEl = {0}  nNod = {1}\n", mC.nEl, mC.nNod);
        
    }

    private void readData() {
        while (CSMGEN.RD.hasNext()) {
            String name = CSMGEN.RD.next().ToLower();
            if (name.Equals("#")) {
                CSMGEN.RD.nextLine();
                continue;
            }
            if (name.Equals("eps"))
                eps = CSMGEN.RD.readDouble();
            else if (name.Equals("end")) break;
            else UTIL.errorMsg("Unexpected data: " + name);
        }
        CSMGEN.PR.Write(" Coordinate error tolerance eps = {0,10:E3}\n", eps);
    }

    // Find coincident nodes in models mA and mB,
    // generate new node numbers for model mB
    private void findCoincidentNodes() 
    {
        newNodesB = new int[mB.nNod];
        int ndim = mA.nDim;
        for (int i = 0; i < mB.nNod; i++) 
            newNodesB[i] = -1;
        
        // Register coincident nodes of mesh B
        //     in array newNodesB
        
        for (int ia = 0; ia < mA.nNod; ia++) 
        {
            double[] xyA = mA.getNodeCoords(ia);
        
            for (int ib = 0; ib < mB.nNod; ib++) 
            {
                for (int j = 0; j < ndim; j++) 
                {
                    if (Math.Abs(xyA[j] - mB.getNodeCoord(ib, j)) > eps)
                    {
                        goto B;
                    }
                }
                newNodesB[ib] = ia;
            B: { }
            }

        }
        nConnected = 0;
        int n = mA.nNod;

        // New node numbers for nodes of model mB
        for (int i = 0; i < mB.nNod; i++) {
            if (newNodesB[i] == -1) newNodesB[i] = n++;
            else nConnected++;
        }
    }

    // Paste two meshes.
    // Nodes and elements of the first mesh
    // are first in the resulting mesh.
    // returns   resulting mesh after pasting.
    private FeModel pasteModels() {

        FeModel mC = new FeModel(CSMGEN.RD, CSMGEN.PR);
        mC.nDim = mA.nDim;

        // nodal coordinates of model mC
        mC.nNod = mA.nNod + mB.nNod - nConnected;
        mC.newCoordArray();
        // Copy nodes of model mA
        for (int i = 0; i < mA.nNod; i++)
            mC.setNodeCoords(i, mA.getNodeCoords(i));
        // Add nodes of model mB
        for (int i = 0; i < mB.nNod; i++)
            mC.setNodeCoords(
                    newNodesB[i], mB.getNodeCoords(i));

        // Element connectivities of model mC
        mC.nEl = mA.nEl + mB.nEl;
        mC.elems = new Element[mC.nEl];
        // Copy elements of model mA
        for (int el = 0; el < mA.nEl; el++) 
        {
            mC.elems[el] = Element.newElement(mA.elems[el].name);
            mC.elems[el].setElemConnectivities(mA.elems[el].ind);
            mC.elems[el].matName = mA.elems[el].matName;
            mC.elems[el].A = mA.elems[el].A;
            mC.elems[el].t = mA.elems[el].t;
            mC.elems[el].I = mA.elems[el].I;
        }
        // Add elements of mB with renumbered connectivities
        for (int el = 0; el < mB.nEl; el++) 
        {
            mC.elems[mA.nEl + el] = Element.newElement(mB.elems[el].name);
            int[] indel = new int[mB.elems[el].ind.Length];
            for (int i = 0; i < mB.elems[el].ind.Length; i++)
                indel[i] = newNodesB[mB.elems[el].ind[i]-1]+1;
            mC.elems[mA.nEl+el].setElemConnectivities(indel);
            mC.elems[mA.nEl+el].matName = mB.elems[el].matName;
            mC.elems[mA.nEl + el].A = mB.elems[el].A;
            mC.elems[mA.nEl + el].t = mB.elems[el].t;
            mC.elems[mA.nEl + el].I = mB.elems[el].I;
        }
        return mC;
    }

}
}
