using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using FEA3D.util;
using FEA3D.model;

namespace FEA3D.fea
{
    // Main class of the mesh generator
    public class CSMGEN
    {

        public static FeScanner RD;
        public static StreamWriter PR;
        //public static HashMap blocks;
        public static Dictionary<String, FeModel> blocks;
        public static String fileOut;

        public static void CSMGEN_Main(String[] args) 
        {
            Console.Clear();
            Console.WriteLine();
            Console.WriteLine("***************************************************");
            Console.WriteLine("*   C# Finite Element Analysis - Mesh Generator   *");
            Console.WriteLine("*   -------------------------------------------   *");
            Console.WriteLine("*   Copyright © Ravinder Singh, 2017.             *");
            Console.WriteLine("*   License: Apache v2.0 License.                 *");
            Console.WriteLine("*   For further information on this software      *");
            Console.WriteLine("*   email : ravi_29_12@hotmail.com.               *");
            Console.WriteLine("***************************************************");

         
            FE.main = FE.MGEN;

            RD = new FeScanner(args[0]);
            string[] file = args[0].Split('.');
            //fileOut = (args.Length == 1) ? file[0] + "_fem" : args[1];
            fileOut = file[0] + "_gen";
            PR = new FePrintWriter().getPrinter(fileOut);

            PR.WriteLine();
            PR.WriteLine("***************************************************");
            PR.WriteLine("*   C# Finite Element Analysis - Mesh Generator   *");
            PR.WriteLine("*   -------------------------------------------   *");
            PR.WriteLine("*   Copyright © Ravinder Singh, 2017.             *");
            PR.WriteLine("*   License: Apache v2.0 License.                 *");
            PR.WriteLine("*   For further information on this software      *");
            PR.WriteLine("*   email : ravi_29_12@hotmail.com.               *");
            PR.WriteLine("***************************************************");
            PR.WriteLine();

            PR.WriteLine("Mesh Generator. Data file: " + args[0]);
            PR.WriteLine();
            Console.WriteLine("\nMesh generator. Data file: " + args[0]);
            Console.WriteLine();

            new CSMGEN();
        }

        public CSMGEN()
        {
            UTIL.printDate(PR);

            // Hash table for storing mesh blocks
            //blocks = new HashMap();
            blocks = new Dictionary<string, model.FeModel>();

            while (RD.hasNext())
            {
                
                String name = RD.next().ToLower();
                if (name.Equals("#")) 
                { 
                    RD.nextLine(); continue; 
                }
                PR.WriteLine("------------------------------------");

                try
                {
                    Activator.CreateInstance(null, "FEA3D.gener." + name);
                    //Class.forName("gener." + name).newInstance();                    
                }
                catch (Exception e)
                {
                    UTIL.errorMsg("Class name not found: " + name + "\n\n\n" + e.ToString());
                    
                }
            }
            PR.Close();
        }

    }
}
