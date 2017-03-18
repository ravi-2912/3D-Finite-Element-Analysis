using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

using FEA3D.elem;
using FEA3D.model;
using FEA3D.solver;
using FEA3D.util;

namespace FEA3D.fea
{
    public class CSFEM
    {

        private static FeScanner RD;
        private static StreamWriter PR;
        public static String fileOut;

        public static void CSFEM_Main(String[] args)
        {
            Console.Clear();
            Console.WriteLine();
            Console.WriteLine("***************************************************");
            Console.WriteLine("*   C# Finite Element Analysis - Pre-Processor    *");
            Console.WriteLine("*                              & Linear Solver    *");
            Console.WriteLine("*   -------------------------------------------   *");
            Console.WriteLine("*   Copyright © Ravinder Singh, 2017.             *");
            Console.WriteLine("*   License: Apache v2.0 License.                 *");
            Console.WriteLine("*   For further information on this software      *");
            Console.WriteLine("*   email : ravi_29_12@hotmail.com.               *");
            Console.WriteLine("***************************************************");

            
            FE.main = FE.FEM;

            RD = new FeScanner(args[0]);
            string[] file = args[0].Split('.');
            //fileOut = (args.Length == 1) ? file[0] + "_fem" : args[1];
            fileOut = file[0] + "_fem";
            PR = new FePrintWriter().getPrinter(fileOut);


            PR.WriteLine();
            PR.WriteLine("***************************************************");
            PR.WriteLine("*   C# Finite Element Analysis - Pre-Processor    *");
            PR.WriteLine("*                              & Linear Solver    *");
            PR.WriteLine("*   -------------------------------------------   *");
            PR.WriteLine("*   Copyright © Ravinder Singh, 2017.             *");
            PR.WriteLine("*   License: Apache v2.0 License.                 *");
            PR.WriteLine("*   For further information on this software      *");
            PR.WriteLine("*   email : ravi_29_12@hotmail.com.               *");
            PR.WriteLine("***************************************************");

            PR.WriteLine();
            PR.WriteLine("FEA Model & Solver. Data file: " + args[0]);
            Console.WriteLine("\nFEA Model & Solver. Data file: " + args[0]);

            new CSFEM();
            PR.Close();
            
        }

        public CSFEM()
        {

            UTIL.printDate(PR);

            FeModel fem = new FeModel(RD, PR);
            Element.fem = fem;

            fem.readData();
            
            //Console.WriteLine(fem.nEl);
            PR.Write(Environment.NewLine + "Number of elements    nEl = {0,9}" + Environment.NewLine +
                       "Number of nodes      nNod = {1,9}" +Environment.NewLine +
                       "Number of dimensions nDim = {2,9}" +Environment.NewLine,
                      fem.nEl, fem.nNod, fem.nDim);
            
            long t0 = Environment.TickCount;

            Solver solver = Solver.newSolver(fem);
            
            solver.assembleGSM();

            PR.WriteLine("Memory for global matrix: {0,9:0.00} MB" + Environment.NewLine, Solver.lengthOfGSM * 8.0e-6);

            FeLoad load = new FeLoad(fem);
            Element.load = load;

            FeStress stress = new FeStress(fem);

            // Load step loop
            while (load.readData()) {
                load.assembleRHS();
                int iter = 0;
                // Equilibrium iterations
                do {
                    iter++;
                    
                    int its = solver.solve(FeLoad.RHS);
                    
                    if (its > 0)
                        PR.Write(Environment.NewLine + "Solver: {0} iterations" + Environment.NewLine, its);
                    
                    stress.computeIncrement();

                } while (!stress.equilibrium(iter));

                stress.accumulate();
                stress.writeResults();

                PR.WriteLine("Loadstep {0}", FeLoad.loadStepName);
                if (iter > 1)
                    PR.WriteLine("{0,5} iterations, " + "Relative residual norm = {1,9:E6}", iter, FeStress.relResidNorm);
                PR.Write(Environment.NewLine);
            }

            PR.Write(Environment.NewLine + "Solution time = {0:0.00} s" + Environment.NewLine, (Environment.TickCount - t0) * 0.001);
            
        }

    }

}