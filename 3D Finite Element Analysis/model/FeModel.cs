using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

using FEA3D.elem;
using FEA3D.util;
using FEA3D.material;
using FEA3D.solver;

namespace FEA3D.model
{
    // Description of the finite element model
    public class FeModel : FeModelData 
    {

        private static int[] elCon = new int[20];
        private static double[][] box= new double[2][]{new double[3], new double[3]};
        IEnumerator<Dof> it;

        // Construct finite element model.
        // RD - data scanner, PR - print writer.
        public FeModel(FeScanner RD, StreamWriter PR) 
        {
            FeModelData.RD = RD;
            FeModelData.PR = PR;

            
        }

        // Read data for a finite element model
        public void readData() 
        {
            readDataFile(RD);
        }

        private void readDataFile(FeScanner es) 
        {
            vars name = vars.NONE;
            String s;
            Material mat;
            it = defDs.GetEnumerator();

            while (es.hasNext()) 
            {
                varName = es.next();
                String varname = varName.ToLower();

                if (varname.Equals("#")) {es.nextLine(); continue; }

                try {
                    name = (vars)System.Enum.Parse(typeof(vars), varname);//vars.valueOf(varname);
                    //Console.WriteLine("hmm");
                } catch (Exception E) {
                    UTIL.errorMsg("Variable name is not found: "+varName);
                }

                switch (name) 
                {
                    case vars.nel:    
                        nEl = es.readInt();
                        //Console.WriteLine("FEMODEL: "+nEl.ToString());
                        break;

                    case vars.nnod:   
                        nNod = es.readInt();
                        //Console.WriteLine("FEMODEL: " + nNod.ToString());
                        break;

                    case vars.ndim:   
                        nDim = es.readInt();
                        nDf = nDim;
                        //Console.WriteLine("FEMODEL: " + nDf.ToString());
                        break;

                    case vars.stressstate:
                        s = es.next().ToLower();
                        try {
                            stressState = (StrStates)System.Enum.Parse(typeof(StrStates), s); //StrStates.valueOf(s);
                        } catch (Exception E) {
                            UTIL.errorMsg("stressState has forbidden value: "+s);
                        }
                        if (stressState != StrStates.threed)
                              nDim = nDf = 2;
                        else  nDim = nDf = 3;
                        //Console.WriteLine("FEMODEL: " + stressState.ToString());
                        break;

                    case vars.physlaw:
                        s = es.next().ToLower();
                        try {
                            physLaw = (PhysLaws)System.Enum.Parse(typeof(PhysLaws), s);//PhysLaws.valueOf(s);
                        } catch (Exception E) {
                            UTIL.errorMsg("physLaw has forbidden value: "+s);
                        }
                        //Console.WriteLine("FEMODEL: " + s.ToString());
                        break;

                    case vars.solver:
                        s = es.next().ToLower();
                        try {
                            Solver.solver = (Solver.solvers)System.Enum.Parse(typeof(Solver.solvers), s); //Solver.solvers.valueOf(s);
                        } catch (Exception E) {
                            UTIL.errorMsg("solver has forbidden value: "+s);
                        }
                        //Console.WriteLine("FEMODEL: " + s.ToString());
                        break;

                    case vars.elcon:
                        readElemData(es);
                        break;

                    case vars.nodcoord:
                        if (nNod == 0 || nDim == 0)
                            UTIL.errorMsg("nNod and nDim should be specified before nodCoord");
                        nEq = nNod * nDim;
                        // Nodal coordinates
                        newCoordArray();
                        for (int i = 0; i < nNod; i++)
                            for (int j = 0; j < nDim; j++)
                                setNodeCoord(i, j, es.readDouble());
                        break;

                    case vars.material:
                        String matname = es.next();
                        mat = Material.newMaterial(physLaw.ToString(), stressState.ToString());
                        double e = es.readDouble();
                        double nu = es.readDouble();
                        double alpha = es.readDouble();
                        mat.setElasticProp(e, nu, alpha);
                        if (physLaw == PhysLaws.elplastic) {
                            double sY = es.readDouble();
                            double km = es.readDouble();
                            double mm = es.readDouble();
                            mat.setPlasticProp(sY, km, mm);
                        }
                        materials.Add(matname, mat);
                        //Console.WriteLine("FEMODEL: " + matname);
                        break;

                    case vars.constrdispl:
                        readConstrDisplacements(es);
                        break;

                    case vars.boxconstrdispl:
                        createBoxConstrDisplacements(es);
                        break;

                    case vars.thermalloading:
                        s = es.next();
                        if (s.ToLower().Equals("y"))
                            thermalLoading = true;
                        else if (s.ToLower().Equals("n"))
                            thermalLoading = false;
                        else
                            UTIL.errorMsg("thermalLoading should be y/n. Specified: " + s);
                        break;

                    case vars.includefile:
                        s = es.next();
                        FeScanner R = new FeScanner(s);
                        readDataFile(R);
                        break;

                    case vars.trussarea:
                        readTrussArea(es);
                        break;

                    case vars.end:  
                        return;
                }
            }
        }

        // Read element type, material and connectivities
        // for all elements
        private void readElemData(FeScanner es) 
        {
            if (nEl == 0) 
                UTIL.errorMsg ("nEl should be defined before elCon");
            elems = new Element[nEl];
            for (int iel = 0; iel < nEl; iel++) {
                // Element type
                String s = es.next().ToLower();
                elems[iel] = Element.newElement(s);
                // Element material
                String elMat = es.next();                
                elems[iel].setElemMaterial(elMat);
                //Element Area or Thickness if Truss or Plane 2D or Shell
                if(s.Equals(Element.elements.truss22.ToString()) || s.Equals(Element.elements.truss32.ToString()))
                    elems[iel].A = es.readDouble();
                else if(s.Equals(Element.elements.quad4.ToString()) || s.Equals(Element.elements.quad8r.ToString()))
                    elems[iel].t = es.readDouble();                 
                
                // Element connectivities
                int nind = elems[iel].ind.Length;
                for (int i = 0; i < nind; i++)
                {
                    elCon[i] = es.readInt();
                }
                elems[iel].setElemConnectivities(elCon, nind);
                
            }
        }

        // Read data for specified constrained displacements
        private void readConstrDisplacements(FeScanner es) 
        {
            String s = es.next().ToLower();
            int idf = UTIL.direction(s);
            if (idf == -1) 
                UTIL.errorMsg("constrDispl direction should be x/y/z. Specified:"+s);
            if (!es.hasNextDouble())
                UTIL.errorMsg("constrDispl value is not a double: " +es.next());
            double vd = es.nextDouble();
            defDs = es.readNumberList(defDs, idf, nDim, vd);
        }

        // Create data for constrained displacements
        // specified inside a box
        private void createBoxConstrDisplacements(FeScanner es) 
        {
            String s = es.next().ToLower();
            int idf = UTIL.direction(s);
            if (idf == -1)
                UTIL.errorMsg("boxConstrDispl direction should be x/y/z. Specified:"+s);
            if (!es.hasNextDouble())
                UTIL.errorMsg("boxConstrDispl value is not a double: " + es.next());
            double vd = es.nextDouble();
            for (int i = 0; i < 2; i++)
                for (int j=0; j<nDim; j++)
                    box[i][j] = es.readDouble();
            //NODE: 
            for (int i = 0; i < nNod; i++) 
            {
                for (int j = 0; j < nDim; j++) 
                {
                    double x = getNodeCoord(i,j);
                    if (x<box[0][j] || x>box[1][j]) 
                        goto NODE;
                }
                defDs.Add(new Dof(nDim * i + idf, vd));
                NODE: { }
            }
        }

        private void readTrussArea(FeScanner es)
        {
            int count = es.readInt();
            String s = es.next().ToLower();
            int idf = UTIL.direction(s);
            if (idf == -1)
                UTIL.errorMsg("constrDispl direction should be x/y/z. Specified:" + s);
            if (!es.hasNextDouble())
                UTIL.errorMsg("constrDispl value is not a double: " + es.next());
            double vd = es.nextDouble();
            defDs = es.readNumberList(defDs, idf, nDim, vd);
        }

    }
}
