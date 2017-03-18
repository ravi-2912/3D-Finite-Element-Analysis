using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

using FEA3D.elem;
using FEA3D.util;

namespace FEA3D.model
{
    // Finite element model data
    public class FeModelData
    {
        public static FeScanner RD;
        protected static StreamWriter PR;

        // Problem dimension =2/3
        public int nDim = 3;
        // Number of degree of freedom per node = 2/3
        public int nDf = 3;
        // Number of nodes
        public int nNod;
        // Number of elements
        public int nEl;
        // Number of degrees of freedom in the FE model
        public int nEq;
        // Elements
        public Element[] elems;
        // Materials
        public Dictionary<object, object> materials = new Dictionary<object, object>(); // Java uses Hashmap
        // Coordinates of nodes
        public double[] xyz;
        // Constrained degrees of freedom
        public List<Dof> defDs = new List<Dof>();
        public bool thermalLoading;
        protected static string varName;

        public enum StrStates
        {
            plstrain, plstress, axisym, threed
        };

        public static StrStates stressState = StrStates.threed;

        public enum PhysLaws
        {
            elastic, elplastic
        };

        public PhysLaws physLaw = PhysLaws.elastic;

        protected enum vars
        {
            nel, nnod, ndim, stressstate, physlaw, solver, elcon, nodcoord, material,
            constrdispl, boxconstrdispl, thermalloading, includefile, user, end, 
            trussarea, NONE
        };

        // Allocation of nodal coordinate array
        public void newCoordArray()
        {
            xyz = new double[nNod * nDim];
        }

        // Set coordinates of node
        public void setNodeCoords(int node, double[] xyzn)
        {
            for (int i = 0; i < nDim; i++) 
                xyz[node * nDim + i] = xyzn[i];
        }

        // Set ith coordinates of node
        public void setNodeCoord(int node, int i, double v)
        {
            xyz[node * nDim + i] = v;
        }

        // Get coordinates of node
        public double[] getNodeCoords(int node) 
        {
            double[] nodeCoord = new double[nDim];
            for (int i = 0; i < nDim; i++)
                nodeCoord[i] = xyz[node * nDim + i];            
            return nodeCoord;            
        }

        // Get coordinates of node VER 2
        public double[] getNodeCoordsV2(int node)
        {
            double[] nodeCoord = new double[nDim];
            for (int i = 0; i < nDim; i++)
                nodeCoord[i] = xyz[node * nDim + i];
            return nodeCoord;
        }

        // Get ith coordinate of node
        public double getNodeCoord(int node, int i)
        {
            return xyz[node * nDim + i];
        }
    }
}
