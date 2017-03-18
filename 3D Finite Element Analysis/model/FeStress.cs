using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

using FEA3D.elem;
using FEA3D.material;
using FEA3D.util;

namespace FEA3D.model
{
    // Stress increment due to displacement increment
    public class FeStress
    {

        public static double relResidNorm;
        private FeModel fem;

        // Constructor for stress increment.
        // fem - finite element model
        public FeStress(FeModel fem)
        {
            this.fem = fem;
        }

        // Compute stress increment for the finite element model
        public void computeIncrement()
        {

            // Accumulate solution vector in displacement increment
            for (int i = 0; i < fem.nEq; i++)
            {
                FeLoad.dDispl[i] += FeLoad.RHS[i];
            }

            // Compute stresses at reduced integration points
            for (int iel = 0; iel < fem.nEl; iel++)
            {
                Element elm = fem.elems[iel];
                elm.setElemXyT();
                elm.disAssembleElemVector(FeLoad.dDispl);

                for (int ip = 0; ip < elm.str.Length; ip++)
                {
                    Material mat = (Material)fem.materials[elm.matName];//.kget(elm.matName);
                    mat.strainToStress(elm, ip);
                }
            }
        }

        // Check equilibrium and assemble residual vector.
        // iter - number of iterations performed
        public bool equilibrium(int iter)
        {
            
            if (fem.physLaw == FeModel.PhysLaws.elastic || iter == FeLoad.maxIterNumber) 
                return true;
            // ALL BELOW FOR ELASTIC-PLASTIC
            // Assemble residual vector to right-hand side
            for (int i = 0; i < fem.nEq; i++)
                FeLoad.RHS[i] = FeLoad.spLoad[i] + FeLoad.dpLoad[i];
            Element elm;
            for (int iel = 0; iel < fem.nEl; iel++)
            {
                elm = fem.elems[iel];
                elm.setElemXy();
                elm.equivStressVector();
                elm.assembleElemVector(Element.evec, FeLoad.RHS);
            }
            // Displacement boundary conditions
            IEnumerator<Dof> it = fem.defDs.GetEnumerator();
            while (it.MoveNext())
            {
                Dof d = (Dof)it.Current;
                FeLoad.RHS[d.dofNum - 1] = 0;
            }
            // Relative residual norm
            double dpLoadNorm = vectorNorm(FeLoad.dpLoad);
            if (dpLoadNorm < 1e-30)
                dpLoadNorm = vectorNorm(FeLoad.dhLoad);
            relResidNorm = vectorNorm(FeLoad.RHS) / dpLoadNorm;
            return relResidNorm < FeLoad.residTolerance;
        }

        // Returns norm of a vector v
        double vectorNorm(double[] v) 
        {

            double norm = 0;
            foreach (double aV in v) 
                norm += aV * aV;
            return Math.Sqrt(norm);
        }

        // Accumulate loads, temperature and stresses
        public void accumulate()
        {

            for (int i = 0; i < fem.nEq; i++)
            {
                FeLoad.spLoad[i] += FeLoad.dpLoad[i];
                FeLoad.sDispl[i] += FeLoad.dDispl[i];
            }
            for (int iel = 0; iel < fem.nEl; iel++)
                fem.elems[iel].accumulateStress();
        }

        // Write results to a file.
        public void writeResults() 
        {

            String fileResult = fea.CSFEM.fileOut + "."+ FeLoad.loadStepName;
            StreamWriter PR = new FePrintWriter().getPrinter(fileResult);

            PR.WriteLine("Displacements" + Environment.NewLine);
            if (fem.nDim == 2)
                PR.WriteLine(" Node             ux             uy");
            else
                PR.WriteLine(" Node             ux             uy             uz");
            for (int i = 0; i < fem.nNod; i++) 
            {
                PR.Write( "{0,5}", i + 1);
                for (int j = 0; j < fem.nDim; j++)
                {
                    PR.Write("{0,15:E6}", FeLoad.sDispl[fem.nDim * i + j]);
                    //Console.Write("{0,15:E6}", FeLoad.sDispl[fem.nDim * i + j]);
                }
                PR.WriteLine();
            }

            PR.WriteLine(Environment.NewLine + Environment.NewLine + "Stresses" );
            for (int iel = 0; iel < fem.nEl; iel++) 
            {
                if (fem.nDim == 2)
                    PR.Write(Environment.NewLine+"El {0,4}     sxx            syy            sxy            szz            epi", iel+1);
                else
                    PR.Write(Environment.NewLine+"El {0,4}     sxx            syy            szz            sxy            syz            szx            epi", iel+1);
                
                foreach (StressContainer aStr in fem.elems[iel].str) 
                {                        
                    PR.Write(Environment.NewLine);
                    for (int i = 0; i < 2 * fem.nDim; i++)
                    {
                        PR.Write("{0,20:0.00000000}", aStr.sStress[i]);
                        //PR.Write("{0,15:0.00000000}", aStr.sStress[i]);
                    }
                    PR.Write("{0,20:0.00000000}", aStr.sEpi);
                }
            }
            PR.Close();
        }

        // Read results from a file.
        // displ - displacements for the finite element model (out)
        public void readResults(String resultFile, double[] displ) 
        {

            if (resultFile==null) return;

            FeScanner RD = new FeScanner(resultFile);
            // Read displacements
            RD.moveAfterLineWithWord("node");
            for (int i = 0; i < fem.nNod; i++) 
            {
                RD.readInt();
                for (int j = 0; j < fem.nDim; j++)
                    displ[fem.nDim*i + j] = RD.readDouble();
            }
            // Read stresses
            for (int iel = 0; iel < fem.nEl; iel++) 
            {
                RD.moveAfterLineWithWord("el");
                foreach (StressContainer aStr in fem.elems[iel].str)
                {
                    for (int i = 0; i < 2 * fem.nDim; i++)
                    {
                        if (fem.nDim == 3)
                            aStr.sStress[i] = RD.readDouble();
                        else
                        {
                            if (i == 2)
                                aStr.sStress[i + 1] = RD.readDouble();
                            else if (i == 3)
                                aStr.sStress[i - 1] = RD.readDouble();
                            else aStr.sStress[i] = RD.readDouble();
                        }
                    }
                        
                    aStr.sEpi = RD.readDouble();
                }
            }
            RD.close();
        }

    }
}
