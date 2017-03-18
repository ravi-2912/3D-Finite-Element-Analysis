using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.model;
using FEA3D.material;
using FEA3D.fea;
using FEA3D.util;

namespace FEA3D.elem
{
    // Finite element
    public abstract class Element 
    {

        // Finite element model
        public static FeModel fem;
        // Finite element load
        public static FeLoad load;
        // Material of current element
        public static Material mat;
        // Element stiffness matrix
        public static double[][] kmat = new double[60][]; // Remaining initialization in the constructor
        // Element vector
        public static double[] evec = new double[60];
        // Element nodal coordinates
        public static double[][] xy = new double[20][]; // Remaining initialization in the constructor
        // Element nodal temperatures
        public static double[] dtn = new double[20];
        // Strain vector
        public static double[] dstrain = new double[6];

        // Element name
        public String name;
        // Element material name
        public String matName;
        // Element connectivities
        public int[] ind;
        // Stress-strain storage
        public StressContainer[] str;

        // For Truss, beam, shell and plane elements
        public double A, t, I;

        // Implemented element types
        public enum elements
        {
            quad4, hex8, 
            quad8r, hex20r,
            truss22, truss32,
            nonexistent
        };

        private static Element create( elements e)
        {
            switch (e)
            {
                case elements.truss22:
                    return new ElemTruss2d2();
                case elements.truss32:
                    return new ElemTruss3d2();

                case elements.quad4:
                    return new ElementLin2D();
                case elements.quad8r:
                    return new ElementQuad2D();
                    
                case elements.hex8:
                    return new ElementLin3D();
                case elements.hex20r:
                    return new ElementQuad3D();

                default:
                    return new NonExistentElement();
            }
        }

        // Construct new element
        // name - element name
        public static Element newElement(String name) 
        {
            elements el = elements.nonexistent;
            try 
            {
                el = (elements)System.Enum.Parse(typeof(elements), name);
            } 
            catch (Exception E) 
            {
                UTIL.errorMsg("Incorrect element type: " + name);
            }
            return create(el);
        }

        // Constructor for an element.
        // name - element name;
        // nind - number of nodes;
        // nstress - number of stress points
        public Element(String name, int nind, int nstress) 
        {
            this.name = name;
            ind = new int[nind];
            if (FE.main != FE.MGEN)
            {
                str = new StressContainer[nstress];
                for (int ip = 0; ip < nstress; ip++)
                    str[ip] = new StressContainer(fem.nDim);
            }
            for (int i = 0; i < 20; i++)
                xy[i] = new double[3];

            for (int i = 0; i < 60; i++)
                kmat[i] = new double[60];
        }

        // Compute element stiffness matrix kmat[][]
        public virtual void stiffnessMatrix() 
        {
            //Console.WriteLine("NOOOO");
        }

        // Compute element thermal vector (evec[])
        public virtual void thermalVector() { }

        // Element nodal equivalent of distributed face load
        // (evec[])
        public virtual int equivFaceLoad(ElemFaceLoad surLd) 
        {
            return -1;
        }

        // Nodal vector equivalent to stresses (evec[])
        public virtual void equivStressVector() { }

        // Get local node numbers for element faces
        // returns elementFaces[nFaces][nNodesOnFace]
        public virtual int[][] getElemFaces() 
        {
            return new int[][] { new int[] { 0 }, new int[] { 0 } };
        }

        // Get strains at integration point (stress)
        // intPoint - integration point number (stress);
        // returns  strain vector [2*ndim]
        public virtual double[] getStrainsAtIntPoint(int intPoint) 
        {
            return new double[] {0.0, 0.0};
        }

        // Get temperature at integration point (stress)
        // intPoint - integration point number (stress);
        // returns  temperature
        public virtual double getTemperatureAtIntPoint(int intPoint) 
        {
            return 0.0;
        }

        // Extrapolate quantity from integration points to nodes
        // fip [nInt][2*nDim] - values at integration points;
        // fn [nind][2*nDim] - values at nodes (out)
        public virtual void extrapolateToNodes(double[][] fip, double[][] fn) { }

        // Set element connectivities
        // indel - connectivity numbers
        // nind - number of element nodes
        public void setElemConnectivities(int[] indel, int nind) 
        {
            System.Array.ConstrainedCopy(indel, 0, ind, 0, nind);
        }

        // Set element connectivities
        // indel - connectivity numbers
        public void setElemConnectivities(int[] indel) 
        {
            System.Array.ConstrainedCopy(indel, 0, ind, 0, indel.Length);
        }

        // Set element material name
        // mat - material name
        public void setElemMaterial(String mat) 
        {
            matName = mat;
        }

        // Set element nodal coordinates xy[nind][nDim]
        public void setElemXy() 
        {
            for (int i = 0; i < ind.Length; i++)
            {
                int indw = ind[i] - 1;
                if (indw >= 0)
                {
                    xy[i] = fem.getNodeCoords(indw);
                }
            }
        }

        // Set nodal coordinates xy[nind][nDim] and
        //     temperatures dtn[nind]
        public void setElemXyT() 
        {
            for (int i = 0; i < ind.Length; i++)
            {
                int indw = ind[i] - 1;
                if (indw >= 0)
                {
                    if (fem.thermalLoading)
                        dtn[i] = FeLoad.dtemp[indw];
                    xy[i] = fem.getNodeCoords(indw);
                }
            }
         }

        // Assemble element vector.
        // elVector - element vector;
        // glVector - global vector (in/out)
        public void assembleElemVector(double[] elVector, double[] glVector) 
        {
            for (int i = 0; i < ind.Length; i++)
            {
                int indw = ind[i] - 1;
                if (indw >= 0)
                {
                    int adr = indw * fem.nDim;
                    for (int j = 0; j < fem.nDim; j++)
                        glVector[adr + j] += elVector[i * fem.nDim + j];
                }
            }
            //Console.WriteLine(Environment.NewLine);
            //for (int i = 0; i < glVector.Length; i++)
            //    Console.WriteLine(glVector[i]);
        }

        // Disassemble element vector (result in evec[]).
        // glVector - global vector
        public void disAssembleElemVector(double[] glVector) 
        {
            for (int i = 0; i < ind.Length; i++)
            {
                int indw = ind[i] - 1;
                if (indw >= 0)
                {
                    int adr = indw * fem.nDim;
                    for (int j = 0; j < fem.nDim; j++)
                        evec[i * fem.nDim + j] = glVector[adr + j];
                }
            }
        }

        // Returns element connectivities
        public int[] getElemConnectivities() 
        {
            int[] indE = new int[ind.Length];
            System.Array.ConstrainedCopy(ind, 0, indE, 0, ind.Length);
            return indE;
        }

        //  Accumulate stresses and equivalent plastic strain
        public void accumulateStress() 
        {
            for (int ip = 0; ip < str.Length; ip++)
                for (int i = 0; i < 2 * fem.nDim; i++)
                    str[ip].sStress[i] += str[ip].dStress[i];
            for (int ip = 0; ip < str.Length; ip++)
                str[ip].sEpi += str[ip].dEpi;

        }
    }
}
