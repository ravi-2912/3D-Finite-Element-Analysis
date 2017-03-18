using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using FEA3D.elem;
using FEA3D.model;
using FEA3D.solver;
using FEA3D.util;

namespace FEA3D.visual
{
    public class VisData 
    {

        static FeModel fem;
        // Global vectot of nodal displacements
        public static double[] displ;

        // Input data names
        public enum vars {
            meshfile, resultfile, parm, showedges, shownodes,
            ndivmin, ndivmax, fmin, fmax, ncontours, deformscale,
            end, none
        }

        // Parameters that can be visualized: displacements,
        // stresses, principal stresses and equivalent stress
        public enum parms {
            ux, uy, uz, sx, sy, sz, sxy, syz, szx,
            s1, s2, s3, si, s13, none
        }

        static String meshFile = null, resultFile = null;
        static parms parm = parms.none;
        public static bool showEdges = true, showNodes = false,
            showDeformShape = false, drawContours = false;
        public static double deformScale = 0.0;
        static int nDivMin = 2, nDivMax = 16;
        static double fMin = 0, fMax = 0;
        static int nContours = 256;

        static float offset = 500.0f;
        static float offsetFactor = 1.0f;


        // Size of the color gradation strip
        static int textureSize = 256;
        // Coefficient for curvature: n = 1 + C*ro
        static double Csub = 15;
        // Coefficient for contours: n = 1 + F*abs(df)/deltaf
        static double Fsub = 20;

        public static void readData(FeScanner RD) {

            readDataFile(RD);

            FeScanner fes = new FeScanner(meshFile);
            fem = new FeModel(fes, null);
            Element.fem = fem;
            fem.readData();

            if (resultFile != null) {
                displ = new double[fem.nNod*fem.nDf];
                FeStress stress = new FeStress(fem);
                stress.readResults(resultFile, displ);
                if (deformScale > 0) showDeformShape = true;
                drawContours = VisData.parm != VisData.parms.none;

            }
        }

        static void readDataFile(FeScanner RD) {

            vars name = vars.none;

            while (RD.hasNext()) {

                String varName = RD.next();
                String varNameLower = varName.ToLower();
                if (varName.Equals("#")) {
                    RD.nextLine();    continue;
                }
                try 
                {
                    name = (vars)System.Enum.Parse(typeof(vars), varNameLower);
                } 
                catch (Exception e) 
                {
                    UTIL.errorMsg("Variable name is not found: " + varName);
                }

                switch (name) 
                {

                    case vars.meshfile:
                        meshFile = RD.next();
                        break;
                    case vars.resultfile:
                        resultFile = RD.next();
                        break;
                    case vars.parm:
                        try 
                        {
                          varName = RD.next();
                          parm = (parms)System.Enum.Parse(typeof(parms), varName.ToLower());
                        } 
                        catch (Exception e) 
                        { UTIL.errorMsg("No such result parameter: " + varName); }
                        break;
                    case vars.showedges:
                        showEdges = RD.next().Equals("y");
                        break;
                    case vars.shownodes:
                        showNodes = RD.next().Equals("y");
                        break;
                    case vars.ndivmin:
                        nDivMin = RD.readInt();
                        break;
                    case vars.ndivmax:
                        nDivMax = RD.readInt();
                        break;
                    case vars.fmin:
                        fMin = RD.readDouble();
                        break;
                    case vars.fmax:
                        fMax = RD.readDouble();
                        break;
                    case vars.ncontours:
                        nContours = RD.readInt();
                        break;
                    case vars.deformscale:
                        deformScale = RD.readDouble();
                        break;
                    case vars.end:
                        return;
                    default:
                        return;
                }
            }
        }

    }
}
