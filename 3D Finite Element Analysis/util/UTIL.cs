using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Globalization;
using System.Runtime.InteropServices;

namespace FEA3D.util
{
    public class UTIL
    {
        // Print error message and exit.
        // message - error message that is printed.
        public static void errorMsg(String message)
        {
            //System.
            Console.WriteLine("=== ERROR: " + message);
            Environment.Exit(1);
        }

        // Transform text directio to integer.
        // s - directio x/y/z/n.
        // returns integer direction 1/2/3/0, error -1.
        public static int direction(String s)
        {
            if (s.Equals("x")) return 1;
            else if (s.Equals("y")) return 2;
            else if (s.Equals("z")) return 3;
            else if (s.Equals("n")) return 0;
            else return -1;
        }

        // Print date and time.
        // PR - PrintWriter for listing file
        public static void printDate(StreamWriter PR)
        {
            Calendar c = new GregorianCalendar();
            DateTime dt =  DateTime.Now;

            PR.WriteLine("Date: {0}-{1}-{2}, Time: {3}:{4}:{5}", 
                c.GetYear(dt), c.GetMonth(dt), c.GetDayOfMonth(dt), 
                c.GetHour(dt), c.GetMinute(dt), c.GetSecond(dt));
        }
    }
}
