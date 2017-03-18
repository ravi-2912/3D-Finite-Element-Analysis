using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using FEA3D.fea;
using FEA3D.elem;
using FEA3D.model;
using FEA3D.util;

namespace FEA3D.visual
{
    public class CsfeaToVtk
    {
        static FeModel fem;
        static FeStress stress;
        static FeScanner RD = null;
        static FeScanner fes = null;
        static StreamWriter WR = null;
        static string meshFile = null;
        static string resultFile = null;
        static string outFile = null;
        static double[] displ;
        static ResultAtNodes resAtNod;

        static int[] vtkArrInd_hex20 = new int[20] { 0, 2, 4, 6, 12, 14, 16, 18, 1, 3, 5, 7, 13, 15, 17, 19, 8, 9, 10, 11 };
        static int[] vtkArrInd_hex8 = new int[8] { 0, 1, 2, 3, 4, 5, 6, 7 };
        static int[] vtkArrInd_quad8 = new int[8] { 0, 2, 4, 6, 1, 3, 5, 7 };
        static int[] vtkArrInd_quad4 = new int[4] { 0, 1, 2, 3 };
        static int[] vtkArrInd_truss2 = new int[2] { 0, 1 };

        public static void CsfeaToVtk_Main(String[] args) 
        {
            Console.Clear();
            Console.WriteLine();
            Console.WriteLine("***************************************************");
            Console.WriteLine("*   C# Finite Element Analysis - Post-Processor   *");
            Console.WriteLine("*   -------------------------------------------   *");
            Console.WriteLine("*   Copyright © Ravinder Singh, 2017.             *");
            Console.WriteLine("*   License: Apache v2.0 License.                 *");
            Console.WriteLine("*   For further information on this software      *");
            Console.WriteLine("*   email : ravi_29_12@hotmail.com.               *");
            Console.WriteLine("***************************************************");
            Console.WriteLine();

            string[] w;
            switch(args.Length)
            {
                case 1:
                    w = args[0].Split('.');
                    meshFile = args[0];
                    outFile = w[0] + "_msh" + ".vtk";
                    break;
                case 2:
                    w = args[1].Split('.');
                    meshFile = args[0];
                    resultFile = args[1];
                    outFile = w[0] + "_" + w[1] + ".vtk";
                    break;
                default:
                    break;
            }
            
            
            fes = new FeScanner(meshFile);
            fem = new FeModel(fes, null);
            Element.fem = fem;
            fem.readData();

            WR = new FePrintWriter().getPrinter(outFile);
            WR.WriteLine("# vtk DataFile Version 3.1");
            WR.WriteLine(outFile);
            WR.WriteLine("ASCII\nDATASET UNSTRUCTURED_GRID");
            WR.WriteLine("POINTS {0} DOUBLE", fem.nNod);

            int ElemConnectNodeCount = 0;
            for (int iel = 0; iel < fem.nEl; iel++)
            {
                if (fem.elems[iel].name.Equals("hex20") || fem.elems[iel].name.Equals("hex20r"))
                {
                    ElemConnectNodeCount += 21;
                }
                if (fem.elems[iel].name.Equals("hex8") || fem.elems[iel].name.Equals("hex8r"))
                {
                    ElemConnectNodeCount += 9;
                }
                if (fem.elems[iel].name.Equals("quad8") || fem.elems[iel].name.Equals("quad8r"))
                {
                    ElemConnectNodeCount += 9;                    
                }
                if (fem.elems[iel].name.Equals("quad4") || fem.elems[iel].name.Equals("quad4r"))
                {
                    ElemConnectNodeCount += 5;
                }
                if (fem.elems[iel].name.Equals("truss22") || fem.elems[iel].name.Equals("truss32"))
                {
                    ElemConnectNodeCount += 3;
                }
            }

            for (int i = 0; i < fem.nNod; i++)
            {
                for (int j = 0; j < fem.nDim; j++)
                    WR.Write("{0,20:0.000000000}", fem.getNodeCoord(i, j));
                if (fem.nDim < 3)
                    WR.Write("{0,20:0.000000000}", 0.0);
                if (fem.nDim < 2)
                    WR.Write("{0,20:0.000000000}", 0.0);
                WR.WriteLine();
            }

            WR.WriteLine();

            WR.Write("\nCELLS {0} ", fem.nEl); 
            WR.Write("{0} \n", ElemConnectNodeCount);

            

            for (int iel = 0; iel < fem.nEl; iel++)
            {
                if (fem.elems[iel].name.Equals("hex20") || fem.elems[iel].name.Equals("hex20r"))
                {
                    WR.Write("20 ");
                    int nind = fem.elems[iel].ind.Length;
                    for (int i = 0; i < vtkArrInd_hex20.Length; i++)
                        WR.Write("{0,6}", fem.elems[iel].ind[vtkArrInd_hex20[i]] - 1);                    
                }

                if (fem.elems[iel].name.Equals("hex8") || fem.elems[iel].name.Equals("hex8r"))
                {
                    WR.Write("8 ");
                    int nind = fem.elems[iel].ind.Length;
                    for (int i = 0; i < vtkArrInd_hex8.Length; i++)
                        WR.Write("{0,6}", fem.elems[iel].ind[vtkArrInd_hex8[i]] - 1);
                }

                if (fem.elems[iel].name.Equals("quad8") || fem.elems[iel].name.Equals("quad8r"))
                {
                    WR.Write("8 ");
                    int nind = fem.elems[iel].ind.Length;
                    for (int i = 0; i < vtkArrInd_quad8.Length; i++)
                        WR.Write("{0,6}", fem.elems[iel].ind[vtkArrInd_quad8[i]] - 1);
                }

                if (fem.elems[iel].name.Equals("quad4") || fem.elems[iel].name.Equals("quad4r"))
                {
                    WR.Write("4 ");
                    int nind = fem.elems[iel].ind.Length;
                    for (int i = 0; i < vtkArrInd_quad4.Length; i++)
                        WR.Write("{0,6}", fem.elems[iel].ind[vtkArrInd_quad4[i]] - 1);
                }

                if (fem.elems[iel].name.Equals("truss22") || fem.elems[iel].name.Equals("truss32"))
                {
                    WR.Write("2 ");
                    int nind = fem.elems[iel].ind.Length;
                    for (int i = 0; i < vtkArrInd_truss2.Length; i++)
                        WR.Write("{0,6}", fem.elems[iel].ind[vtkArrInd_truss2[i]] - 1);
                }

                WR.WriteLine();
            }
            WR.WriteLine("\nCELL_TYPES {0}", fem.nEl);
            for (int i = 0; i < fem.nEl; i++)
            {
                if (fem.elems[i].name.Equals("hex20") || fem.elems[i].name.Equals("hex20r"))
                    WR.Write("{0,3}", 25);

                if (fem.elems[i].name.Equals("hex8") || fem.elems[i].name.Equals("hex8r"))
                    WR.Write("{0,3}", 12);

                if (fem.elems[i].name.Equals("quad8") || fem.elems[i].name.Equals("quad8r"))
                    WR.Write("{0,3}", 23);

                if (fem.elems[i].name.Equals("quad4") || fem.elems[i].name.Equals("quad4r"))
                    WR.Write("{0,3}", 9);

                if (fem.elems[i].name.Equals("truss22") || fem.elems[i].name.Equals("truss32"))
                    WR.Write("{0,3}", 3);
            }
                

            WR.WriteLine("\n");

            if(resultFile != null)
            {
                
                displ = new double[fem.nNod * fem.nDf];
                stress = new FeStress(fem);
                stress.readResults(resultFile, displ);
                resAtNod = new ResultAtNodes(fem);

                WR.Write("POINT_DATA {0}", fem.nNod);


                
                writeScalars("Si");
                writeScalars("Sx");
                
                writeScalars("Sy");
                writeScalars("Sz");
                writeScalars("Sxy");
                
                if (fem.nDim == 3)
                {
                    writeScalars("Syz");
                    writeScalars("Szx");
                }

                //writeScalars("S1");
                //writeScalars("S2");
                //if (fem.nDim == 3)
                //{
                //    writeScalars("S3");
                //    writeScalars("S13");
                //}

                writeScalars("Ux");
                writeScalars("Uy");

                if (fem.nDim == 3)
                {
                    writeScalars("Uz");
                }

            }
            
            WR.Flush();
            WR.Close();
        }

        private static void writeScalars(string s, string s2=null)
        {
            VisData.parms vd;
            if (s2 == null)
                vd = (VisData.parms)System.Enum.Parse(typeof(VisData.parms), s.ToLower());
            else vd = (VisData.parms)System.Enum.Parse(typeof(VisData.parms), s2.ToLower());
            WR.WriteLine("\n");
            WR.WriteLine("SCALARS " + s + " double 1");
            WR.WriteLine("LOOKUP_TABLE default");
            WR.WriteLine();
            resAtNod.setParmAtNodes(vd, displ);
            for (int i = 0; i < resAtNod.sgfun.Length; i++)
            {
                WR.WriteLine("{0,20:0.000000000}", resAtNod.sgfun[i]);
            }
            //Console.WriteLine(vd.ToString() + "  :  " + s + "  " + s2);
        }
        public CsfeaToVtk()
        {

        }

    }
}
