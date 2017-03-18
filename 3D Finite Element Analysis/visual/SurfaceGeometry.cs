using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.model;
using FEA3D.util;

namespace FEA3D.visual
{
    // Geometry: surface faces, edges and nodes
    class SurfaceGeometry {

        private static FeModel fem;

        // Numbers of surface element faces, edges and nodes
        int nFaces, nEdges, nsNodes, nElements;
        // Surface element faces and edges, surface nodes
        public List<int[]> listFaces;
        public List<int> listSurfaceElement;
        public List<int[]> listEdges;
        public List<int[]> listSharedFaces;
        int[] sNodes;

        double fmin, fmax, deltaf;
        double[] fun;
        private static double[] xyzmin = new double[3];
        private static double[] xyzmax = new double[3];
        static double sizeMax;

        public List<int[]> angleSurface = new List<int[]>();
        public List<int> anglSurfaceElems = new List<int>();

        public SurfaceGeometry(FeModel f) {

            SurfaceGeometry.fem = f;
            listFaces = new List<int[]>();
            listEdges = new List<int[]>();
            listSurfaceElement = new List<int>();
            listSharedFaces = new List<int[]>();

            sNodes = new int[fem.nNod];

            // Create element faces located at the surface
            createFaces();
            nFaces = listFaces.Count;
            nElements = listSurfaceElement.Count;
            // Create element edges located at the surface
            createEdges();
            nEdges = listEdges.Count();


            // Create nodes located at the surface
            createNodes();

            //if (VisData.drawContours) 
            {
                //fun = new double[fem.nNod];
                //ResultAtNodes ran = new ResultAtNodes(this, fem);
                //ran.setParmAtNodes(VisData.parm, VisData.displ);
            }

            //modifyNodeCoordinates();
        }

        // Create linked list listFaces containing element faces
        // located on the model surface. 2D case: element = face
        void createFaces() 
        {
            if (fem.nDim == 3) 
            {  // 3D mesh
                for (int iel = 0; iel < fem.nEl; iel++)
                {
                    int[][] elemFaces  = fem.elems[iel].getElemFaces();
                    foreach (int[] elemFace in elemFaces) 
                    {
                        int nNodes = elemFace.Length;
                        int[] faceNodes = new int[nNodes];
                        for (int i = 0; i < nNodes; i++) 
                        {
                            faceNodes[i] = fem.elems[iel].ind[elemFace[i]];
                        }
                        // Zero area degenerated 8-node face
                        //if (nNodes == 8 && (faceNodes[3] == faceNodes[7] || faceNodes[1] == faceNodes[5]))
                        //    continue;
                        IEnumerator<int[]> f = listFaces.GetEnumerator();
                        IEnumerator<int> e = listSurfaceElement.GetEnumerator();
                        bool faceFound = false;
                        while (f.MoveNext())
                        {
                            e.MoveNext();
                            int[] faceNodesA = f.Current;
                            if (equalFaces(faceNodes,faceNodesA)) 
                            {
                                listSharedFaces.Add(f.Current);
                                listFaces.Remove(f.Current);
                                listSurfaceElement.Remove(e.Current);
                                f.Dispose();
                                e.Dispose();
                                faceFound = true;
                                break;
                            }
                        }
                        if (!faceFound)
                        {
                            listFaces.Add(faceNodes);
                            listSurfaceElement.Add(iel);
                        }
                    }
                }                
            }
            else 
            {  // 2D - faces = elements
                IEnumerator<int[]> f = listFaces.GetEnumerator();
                for (int iel = 0; iel < fem.nEl; iel++) 
                {
                    listFaces.Add(fem.elems[iel].ind);
                    //surfaceElement.Add(iel);
                }
            }
        }

        // Compare two element faces.
        // Surface has 8 or 4 nodes, corners are compared.
        // f1 - first face connectivities.
        // f2 - second face connectivities.
        // returns  true if faces are same.
        bool equalFaces(int[] f1, int[] f2) {

            // Quadratic elements or linear elements
            int step = (f1.Length > 4) ? 2 : 1;

            for (int j = 0; j < f1.Length; j += step) 
            {
                int n1 = f1[j];
                bool nodeFound = false;
                for (int i = 0; i < f2.Length; i += step) 
                {
                    if (f2[i] == n1) 
                    {
                        nodeFound = true;
                        break;
                    }
                }
                if (!nodeFound) return false;
            }
            return true;
        }

        // Create linked list listEdges containing element edges
        // located on the model surface
        void createEdges() {

            for (int iFace = 0; iFace < nFaces; iFace++) {

                int[] faceNodes = (int[]) listFaces[iFace];
                int nFaceNodes = faceNodes.Length;
                int step = (nFaceNodes > 4) ? 2 : 1;

                for (int inod=0; inod < nFaceNodes; inod += step) {
                    int[] edgeNodes = new int[step + 1];
                    for (int i = inod, k = 0; i <= inod+step;
                         i++, k++)
                        edgeNodes[k] = faceNodes[i%nFaceNodes];

                    IEnumerator<int[]> ea = listEdges.GetEnumerator();
                    bool edgeFound = false;
                    while (ea.MoveNext()) {
                        int[] edgeNodesA = (int[]) ea.Current;
                        if (equalEdges(edgeNodes, edgeNodesA)) {
                            edgeFound = true;
                            break;
                        }
                    }
                    if (!edgeFound) listEdges.Add(edgeNodes);
                }
            }
        }

        // Compare two element edges.
        // e1 - first edge connectivities.
        // e2 - second edge connectivities.
        // returns  true if edges have same node numbers at ends
        bool equalEdges(int[] e1, int[] e2) {

            int len = e1.Length - 1;
            return (e1[0] == e2[0] && e1[len] == e2[len]) ||
                    (e1[0] == e2[len] && e1[len] == e2[0]);
        }

        // Fill out array of surface nodes sNodes (0/1).
        void createNodes() {

            for (int i = 0; i < sNodes.Length; i++) sNodes[i] = 0;

            IEnumerator<int[]> e = listEdges.GetEnumerator();

            for (int iEdge = 0; iEdge < nEdges; iEdge++) {
                e.MoveNext();
                int[] edgeNodes = (int[]) e.Current;
                int nEdgeNodes = edgeNodes.Length;
                for (int i = 0; i < nEdgeNodes; i++)
                    sNodes[edgeNodes[i] - 1] = 1;
            }
            nsNodes = 0;
            foreach (int sNode in sNodes)
                if (sNode > 0) nsNodes++;
        }

        // Add scaled displacements to nodal coordinates and
        //  center finite element mesh
        void modifyNodeCoordinates() {

            // Deformed shape: add scaled displacements
            // to nodal coordinates
            if (VisData.showDeformShape) {
                setBoundingBox();
                double displMax = 0;
                for (int i = 0; i < fem.nNod; i++) {
                    double d = 0;
                    for (int j = 0; j < fem.nDim; j++) {
                        double s =  VisData.displ[i*fem.nDim+j];
                        d += s*s;
                    }
                    displMax = Math.Max(d, displMax);
                }
                displMax = Math.Sqrt(displMax);
                // Scale for visualization of deformed shape
                double scaleD =
                        sizeMax*VisData.deformScale/displMax;
                for (int i = 0; i < fem.nNod; i++) {
                    for (int j = 0; j < fem.nDim; j++)
                        fem.setNodeCoord(i, j,
                            fem.getNodeCoord(i, j) +
                            scaleD*VisData.displ[i*fem.nDim+j]);
                }
            }

            setBoundingBox();
            // Translate JFEM model to have the bounding
            //  box center at (0,0,0).
            double[] xyzC = new double[3];
            for (int j = 0; j < 3; j++)
                xyzC[j] = 0.5*(xyzmin[j] + xyzmax[j]);
            for (int i = 0; i < fem.nNod; i++)
                for (int j = 0; j < fem.nDim; j++)
                    fem.setNodeCoord(i, j,
                        fem.getNodeCoord(i, j) - xyzC[j]);
        }

        // Set min-max values of xyz coordinates of JFEM model
        // xyzmin[] and xyzmax[].
        void setBoundingBox() {

            for (int j = 0; j < fem.nDim; j++) {
                xyzmin[j] = fem.getNodeCoord(0, j);
                xyzmax[j] = fem.getNodeCoord(0, j);
            }
            for (int i = 1; i < fem.nNod; i++) {
                if (sNodes[i] >= 0) {
                    for (int j = 0; j < fem.nDim; j++) {
                        double c = fem.getNodeCoord(i, j);
                        xyzmin[j] = Math.Min(xyzmin[j], c);
                        xyzmax[j] = Math.Max(xyzmax[j], c);
                    }
                }
            }
            if (fem.nDim == 2) {
                xyzmin[2] = -0.01;
                xyzmax[2] =  0.01;
            }
            sizeMax = 0;
            for (int i = 0; i < 3; i++) {
                double s = xyzmax[i] - xyzmin[i];
                sizeMax = Math.Max(s, sizeMax);
            }
        }

        // Compute scale for the finite element model.
        // returns  scale value.
        double getScale() {

            if (sizeMax > 0) return 0.8/sizeMax;
            else return 1.0;
        }

        public double getDotProduct(double[] v1, double[] v2)
        {
            if(v1.Length==v2.Length)
            {
                double vd = 0;
                for (int i = 0; i < v1.Length; i++)
                {
                    vd += v1[i] * v2[i];
                }
                return vd;
            }
            else
            {
                Console.WriteLine("Unequal vector length for [ getDotProduct ].");
                return 0;
            }
        }

        public double getAngleBwVector(double[] v1, double[] v2)
        {
            double vd = getDotProduct(v1, v2);
            double v1d = getDotProduct(v1, v1);
            double v2d = getDotProduct(v2, v2);
            return Math.Acos(vd / (v1d * v2d));
        }

        public double[] getCrossProduct(double[] u, double[] v)
        {
            if (u.Length == v.Length && u.Length == 3)
            {
                double[] v3 = new double[3] { 0, 0, 0 };
                v3[0] = u[1] * v[2] - u[2] * v[1];
                v3[1] = u[2] * v[0] - u[0] * v[2];
                v3[2] = u[0] * v[1] - u[1] * v[0];
                return v3;
            }
            else
            {
                Console.WriteLine("Unequal vector length or length not equal to 3 for [ getCrossProduct ].");
                return new double[] { -1, -1, -1 };
            }
        }

        public double[] getSumVector(double[] u, double[] v)
        {
            double[] x = new double[] { 0, 0, 0 };
            x[0] = u[0] + v[0];
            x[1] = u[1] + v[1];
            x[2] = u[2] + v[2];
            return x;
        }

        public double[] normalizeVector(double[] u)
        {
            if (u.Length == 3) 
            {
                double ud = Math.Sqrt(getDotProduct(u, u));
                u[0] /= ud;
                u[1] /= ud;
                u[2] /= ud;
                return u;
            }
            else
            {
                Console.WriteLine("Vector length not equal to 3 for [ getVectorNormal ].");
                return new double[] { -1, -1, -1 };
            }
        }

        public double[] getFaceNormalVector(int[] faceNodes)
        {
            double[] fn = new double[3] { 0, 0, 0 };
            for (int i = 0; i < faceNodes.Length; i++) 
            {
                double[] current = fem.getNodeCoords(faceNodes[i]-1);
                //Console.WriteLine(current[0] + " " + current[1] + " " + current[2]);
                double[] next = fem.getNodeCoords(faceNodes[(i+1)%faceNodes.Length]-1);
                fn[0] += (current[1] - next[1]) * (current[2] + next[2]);
                fn[1] += (current[2] - next[2]) * (current[0] + next[0]);
                fn[2] += (current[0] - next[0]) * (current[1] + next[1]);

            }
            return normalizeVector(fn);
        }

        public bool sharedEdge(int[] f1, int[] f2)
        {
            int step1 = (f1.Length > 4) ? 2 : 1;
            int step2 = (f2.Length > 4) ? 2 : 1;
            int f = 0;
            for (int i = 0; i < f1.Length; i += step1)
            {
                for (int j = 0; j < f2.Length; j += step2)
                {
                    if (f1[i] == f2[j])
                        f++;
                    if (f == 2)
                        return true;
                }
            }
            return false;
        }

        public void getSurfacebyAngle(int[] face, double angle)
        {
            List<int[]> surfaces = listFaces;
            List<int> surfacesElems = listSurfaceElement;
            IEnumerator<int> ilfe = surfacesElems.GetEnumerator();

            foreach(int[] sf in surfaces)
            {
                ilfe.MoveNext();
                if (sf.Intersect(face).Count() >= 4)
                {
                    angleSurface.Add(sf);
                    anglSurfaceElems.Add(ilfe.Current);
                    if (!surfaces.Remove(sf))
                        Console.WriteLine("Could not remove input face in [ getSurfacebyAngle ]");
                    if (!surfacesElems.Remove(ilfe.Current))
                        Console.WriteLine("Could not remove input face's element in [ getSurfacebyAngle ]");
                    
                    break;
                }
            }

            List<int[]> faceSharingEdge = new List<int[]>();
            ilfe.Dispose();            

        FACE: { }
            ilfe.Dispose();
            ilfe = surfacesElems.GetEnumerator();
            foreach(int[] f in surfaces)
            {
                ilfe.MoveNext();
                faceSharingEdge.Clear();
                foreach (int[] asf in angleSurface)
                    if (sharedEdge(f, asf))
                        faceSharingEdge.Add(asf);

                foreach (int[] fse in faceSharingEdge)
                {
                    double[] fsen = getFaceNormalVector(fse);
                    double[] fn = getFaceNormalVector(f);
                    double ang = (180 / Math.PI) * (getAngleBwVector(fn, fsen));
                    if (ang <= angle)
                    {
                        angleSurface.Add(f);
                        anglSurfaceElems.Add(ilfe.Current);
                        surfaces.Remove(f);
                        surfacesElems.Remove(ilfe.Current);
                        goto FACE;
                    }
                }                
            } 
        }

        

    }
}

/*
Begin Function CalculateSurfaceNormal (Input Polygon) Returns Vector

   Set Vertex Normal to (0, 0, 0)

   Begin Cycle for Index in [0, Polygon.vertexNumber)

      Set Vertex Current to Polygon.verts[Index]
      Set Vertex Next    to Polygon.verts[(Index plus 1) mod Polygon.vertexNumber]

      Set Normal.x to Sum of Normal.x and (multiply (Current.y minus Next.y) by (Current.z plus Next.z))
      Set Normal.y to Sum of Normal.y and (multiply (Current.z minus Next.z) by (Current.x plus Next.x))
      Set Normal.z to Sum of Normal.z and (multiply (Current.x minus Next.x) by (Current.y plus Next.y))

   End Cycle

   Returning Normalize(Normal)

End Function
*/