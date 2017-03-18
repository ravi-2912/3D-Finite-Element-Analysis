using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FEA3D.elem;
using FEA3D.fea;
using FEA3D.util;
using FEA3D.visual;

namespace FEA3D.model
{
    // Load increment for the finite element model
    public class FeLoad : FeLoadData  
    {

        // Finite element model
        private static FeModel fem;
        SurfaceGeometry sg;
        private IEnumerator<Dof> itnf;
        private IEnumerator<ElemFaceLoad> itsf;

        // Construct finite element load.
        // fem - finite element model
        public FeLoad(FeModel fem) 
        {
            FeLoad.fem = fem;
            RD = FeModel.RD;

            spLoad = new double[fem.nEq];
            dpLoad = new double[fem.nEq];
            dhLoad = new double[fem.nEq];
            sDispl = new double[fem.nEq];
            dDispl = new double[fem.nEq];
            RHS    = new double[fem.nEq];

            if (fem.thermalLoading)
                dtemp = new double[fem.nNod];
            sg = new SurfaceGeometry(fem);
        }

        // Read data describing load increment.
        // returns  true if load data has been read
        public bool readData( ) 
        {
            return readDataFile(RD, true);
        }

        // Read data fragment for load increment.
        // newLoad = true - beginning of new load,
        //         = false - continuation of load.
        // returns  true if load data has been read
        private bool readDataFile(FeScanner es, bool newLoad) 
        {
            if (newLoad) {
                scaleLoad = 0;
                nodForces = new List<Dof>();
                itnf = nodForces.GetEnumerator();
                surForces = new List<ElemFaceLoad>();
                itsf = surForces.GetEnumerator();
                
                if (fem.thermalLoading) {
                    for (int i = 0; i < dtemp.Length; i++)
                        dtemp[i] = 0.0;
                }
                for (int i=0; i<dDispl.Length; i++) dDispl[i] = 0.0;
            }

            if (!es.hasNext()) 
                return false;  // No load data

            vars name = vars.NONE;
            String s;

            while (es.hasNext()) 
            {
                String varName = es.next();
                String varNameLower = varName.ToLower();

                if (varName.Equals("#")) {
                    es.nextLine(); 
                    continue; 
                }
                    
                try {
                    name = (vars)System.Enum.Parse(typeof(vars), varNameLower); //vars.valueOf(varNameLower);
                } catch (Exception E) 
                {
                    UTIL.errorMsg("Variable name is not found: " + varName);
                }

                switch (name) 
                {
                    case vars.loadstep:
                        loadStepName = es.next();
                        //Console.WriteLine("FELOAD " + loadStepName.ToString());
                        break;

                    case vars.scaleload:
                        scaleLoad = es.readDouble();
                        //Console.WriteLine("FELOAD  " + scaleLoad.ToString());
                        break;

                    case vars.residtolerance:
                        residTolerance = es.readDouble();
                        //Console.WriteLine("FELOAD  " + residTolerance.ToString());
                        break;

                    case vars.maxiternumber:
                        maxIterNumber = es.readInt();
                        //Console.WriteLine("FELOAD  " + maxIterNumber.ToString());
                        break;

                    case vars.nodforce:
                        readNodalForces(es);
                        break;

                    case vars.surforce:
                        readSurForces(es);
                        break;

                    case vars.boxsurforce:
                        createBoxSurForces(es);
                        break;

                    case vars.nodtemp:
                        dtemp = new double[fem.nNod];
                        for (int i = 0; i < fem.nNod; i++)
                            dtemp[i] = es.readDouble();
                        break;

                    case vars.includefile:
                        s = es.next().ToLower();
                        FeScanner R = new FeScanner(s);
                        readDataFile(R, false);
                        break;

                    case vars.anglesurforce:
                        createAngleSurForce(es);
                        break;

                    case vars.surforcelist:
                        createSurForceList(es);
                        break;

                    case vars.elemanglesurforce:
                        createElemAngleSurForce(es);
                        break;

                    case vars.end:
                        return true;
                }
            }
            return true;
        }


        // Read data for specified nodal forces
        private void readNodalForces(FeScanner es) 
        {
            String s = es.next().ToLower();
            int idf = UTIL.direction(s);
            if (idf == -1) 
                UTIL.errorMsg("nodForce direction should be x/y/z. Specified:"+s);

            if (!es.hasNextDouble()) 
                UTIL.errorMsg("nodForce value is not a double: " + es.next());
            double vd = es.nextDouble();

            nodForces = es.readNumberList(nodForces, idf, fem.nDim, vd);
        }

        // Read data for surface forces (element face loading):
        // direction, iel, nFaceNodes, faceNodes, forcesAtNodes.
        private void readSurForces(FeScanner es) 
        {
            String s = es.next().ToLower();
            int dir = UTIL.direction(s);
            if (dir == -1) 
                UTIL.errorMsg("surForce direction should be x/y/z/n. Specified:"+s);
            int iel = es.readInt();
            int nFaceNodes = es.readInt();
            for (int i=0; i<nFaceNodes; i++)
                iw[i] = es.readInt();
            for (int i=0; i<nFaceNodes; i++)
                dw[i] = es.readDouble();
            surForces.Add(new ElemFaceLoad(iel-1, nFaceNodes, dir, iw, dw));
        }

        // Create data for distributed surface load
        // specified inside a box
        private void createBoxSurForces(FeScanner es) 
        {
            int[][] faces;
            String s = es.next().ToLower();
            int dir = UTIL.direction(s);
            if (dir == -1)
                UTIL.errorMsg("boxSurForce direction should be x/y/z/n. Specified:" + s);

            if (!es.hasNextDouble()) 
                UTIL.errorMsg("boxSurForce value is not a double: " + es.next());
            double force = es.nextDouble();

            for (int i = 0; i < 2; i++)
                for (int j = 0; j < fem.nDim; j++)
                    box[i][j] = es.readDouble();

            for (int iel=0; iel<fem.nEl; iel++) 
            {
                Element el = fem.elems[iel];
                faces = el.getElemFaces();
                //FACE:
                foreach(int[] face in faces) 
                {
                    int nNodes = face.Length;
                    for (int inod = 0; inod < nNodes; inod++)
                        iw[inod] = 0;
                    for (int inod = 0; inod < nNodes; inod++) 
                    {
                        int iGl = el.ind[face[inod]];
                        if (iGl > 0) 
                        {
                            for (int j = 0; j < fem.nDim; j++)  
                            {
                                double x = fem.getNodeCoord(iGl-1,j);
                                if (x < box[0][j] || x > box[1][j])
                                    goto FACE;
                            }
                            iw[inod] = iGl;
                        }
                    }
                    surForces.Add(new ElemFaceLoad(iel,nNodes,dir,iw,force));
                    FACE: { }
                }
            }
        }

        // Assemble right-hand side of the global equation system
        public void assembleRHS() {

            if (scaleLoad != 0.0) {
                for (int i = 0; i < fem.nEq; i++) {
                    dpLoad[i] *= scaleLoad;
                    dhLoad[i] *= scaleLoad;
                    RHS  [i] = dpLoad[i] + dhLoad[i];
                }
                return;
            }
            for (int i = 0; i < fem.nEq; i++) {
                dpLoad[i] = 0.0;
                dhLoad[i] = 0.0;
            }
            
            // Nodal forces specified directly
            itnf = nodForces.GetEnumerator();
            Dof d;
            while (itnf.MoveNext()) {
                d = (Dof) itnf.Current;
                dpLoad[d.dofNum-1] = d.value;
            }

            

            // Surface load at element faces
            itsf = surForces.GetEnumerator();
            IEnumerator<ElemFaceLoad> itsf2 = surForces.GetEnumerator();
            //IEnumerable<ElemFaceLoad> itsf = surForces.GetEnumerator();
            ElemFaceLoad efl;
            Element elm;
            while (itsf.MoveNext())
            {
                efl = (ElemFaceLoad)itsf.Current;
                elm = fem.elems[efl.iel];
                elm.setElemXy();
                if (elm.equivFaceLoad(efl)==-1)
                    UTIL.errorMsg("surForce does not match any face of element: " + efl.iel);
                elm.assembleElemVector(Element.evec,dpLoad);
            }

            // Temperature field
            if (fem.thermalLoading) {
                for (int iel = 0; iel < fem.nEl; iel++) {
                    elm = fem.elems[iel];
                    elm.setElemXyT();
                    elm.thermalVector();
                    elm.assembleElemVector(Element.evec,dhLoad);
                }
            }

            // Right-hand side = actual load + fictitious load
            for (int i = 0; i < fem.nEq; i++)
            {
                RHS[i] = dpLoad[i] + dhLoad[i];
                
            }

            // Displacement boundary conditions for right-hand side
            IEnumerator<Dof> itdbc = fem.defDs.GetEnumerator();
            while (itnf.MoveNext()) 
            {
                d = (Dof)itdbc.Current;
                RHS[d.dofNum-1] = FE.bigValue * d.value;
            }           
        }


        void createAngleSurForce(FeScanner es)
        {
            //
            int nFaceNodes = es.readInt();
            int[] fnode = new int[nFaceNodes];
            for (int i = 0; i < nFaceNodes; i++)
                fnode[i] = es.readInt();
            double angle = es.readDouble();
            
            String s = es.next().ToLower();
            int dir = UTIL.direction(s);
            if (dir == -1)
                UTIL.errorMsg("surForce direction should be x/y/z/n. Specified:" + s);
            double force = es.readDouble();

            sg.getSurfacebyAngle(fnode, angle);
            IEnumerator<int> iase = sg.anglSurfaceElems.GetEnumerator();
            IEnumerator<int[]> iasf = sg.angleSurface.GetEnumerator();
            //Console.WriteLine(sg.angleSurface.Count);
            for (int i = 0; i < sg.angleSurface.Count; i++) 
            {
                iase.MoveNext();
                iasf.MoveNext();
                //Console.WriteLine(iase.Current + " : " + iasf.Current[0] + " " + iasf.Current[1] + " " + iasf.Current[2] + " " + iasf.Current[3] + " ");
                surForces.Add(new ElemFaceLoad(iase.Current, iasf.Current.Length, dir, iasf.Current, force));
            }
        }

        void createSurForceList(FeScanner es)
        {
            int list_count = es.readInt();
            for (int j = 0; j < list_count; j++)
                readSurForces(es);
            
        }

        void createElemAngleSurForce(FeScanner es)
        {
            String s = es.next().ToLower();
            int dir = UTIL.direction(s);
            if (dir == -1)
                UTIL.errorMsg("surForce direction should be x/y/z/n. Specified:" + s);

            double force = es.readDouble();
            double angle = es.readDouble();

            int nFaceNodes = es.readInt();
            int[] f = new int[nFaceNodes];
            for (int i = 0; i < nFaceNodes; i++)
                f[i] = es.readInt();

            int nElem = es.readInt();
            int[] e = new int[nElem];
            for (int i = 0; i < nElem; i++)
                e[i] = es.readInt()-1;

            sg.getSurfacebyAngle(f, angle);
            IEnumerator<int> iase = sg.anglSurfaceElems.GetEnumerator();
            IEnumerator<int[]> iasf = sg.angleSurface.GetEnumerator();

            //foreach(int[] fss in sg.angleSurface)
            //{
            //    iase.MoveNext();
            //    Console.Write(iase.Current + " : ");
            //    for (int k = 0; k < fss.Length; k++)
            //        Console.Write("{0,6}", fss[k]);
            //    Console.WriteLine();
            //}


            foreach (int el in sg.anglSurfaceElems)
            {
                iasf.MoveNext();
                for (int j = 0; j < e.Length; j++)
                {
                    
                    if (el == e[j])
                    {

                        surForces.Add(new ElemFaceLoad(el, iasf.Current.Length, dir, iasf.Current, force));
                    }
                }
                
            }
        }
    }
}
