
using System;
using System.IO;

using FEA3D.fea;
using FEA3D.visual;

namespace FEA3D
{
    public class FEA3D
    {
        enum EXT
        { 
            MGEN, FEM, VIS
        }

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                Console.WriteLine("\nUsage: fea3d [FileIn].gen");
                Console.WriteLine(" OR    fea3d [FileIn].fem");
                Console.WriteLine(" OR    fea3d [MeshFile].msh [FileIn]_fem.[LoadStep]");
                return;
            }
            string ext = Path.GetExtension(args[0]);
            EXT e = EXT.MGEN;
            if (args.Length == 1)
            {
                if (string.Equals(ext, ".gen"))
                    e = EXT.MGEN;
                if (string.Equals(ext, ".fem"))
                    e = EXT.FEM;
                if (string.Equals(ext, ".msh"))
                    e = EXT.VIS;
            }
            if (args.Length == 2)
                e = EXT.VIS;


            switch(e)
            {
                case EXT.MGEN:
                    CSMGEN.CSMGEN_Main(args);
                    break;

                case EXT.FEM:
                    CSFEM.CSFEM_Main(args);
                    break;

                case EXT.VIS:
                    CsfeaToVtk.CsfeaToVtk_Main(args);
                    break;
            }

        }

        
        
    }
}