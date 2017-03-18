using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FEA3D.util
{
    // Finite element printer to file
    public class FePrintWriter
    {
        StreamWriter PR;
        public StreamWriter getPrinter(String fileOut)
        {
            try
            {
                PR = new StreamWriter(new BufferedStream(new FileStream(fileOut, FileMode.Create,FileAccess.Write)));
                //PR.AutoFlush = true;
            }
            catch(Exception E)
            {
                UTIL.errorMsg("Cannot open output file: " + fileOut);
            }
            return PR;
        }
    }
}
