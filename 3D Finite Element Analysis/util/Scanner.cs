using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FEA3D.util
{
    public class Scanner
    {
        private StreamReader es;
        private char[] delimiters;
        private String line;
        private String[] str;
        private String fileName;
        private List<String> strL;
        private IEnumerator<String> it;
        private bool nextHasReadLine = false;

        public Scanner(String fileIn)
        {
            fileName = fileIn;
            try
            {
                es = new StreamReader(new BufferedStream(new FileStream(fileName, FileMode.Open, FileAccess.Read)));
            }
            catch (Exception E)
            {
                UTIL.errorMsg("File not found: " + fileIn);
            }
            read();
        }

        private void read()
        {
            try
            {
                do
                {
                    line = es.ReadLine();
                    str = line.Split(delimiters, StringSplitOptions.RemoveEmptyEntries);
                } while (str.Length <= 0 && es.Peek() >= 0);
                strL = new List<String>(str);
                it = strL.GetEnumerator();
                it.MoveNext();
            }
            catch (Exception E)
            {
                Console.WriteLine("cannot read file " + fileName);
            }

        }

        public void useDelimiters(char[] delims)
        {
            delimiters = delims;
        }

        public bool hasNext()
        {
            return it.Current != null;
        }

        public String next()
        {
            try
            {
                return it.Current;
            }
            finally
            {
                if (!it.MoveNext() && es.Peek() >= 0)
                {
                    read();
                    nextHasReadLine = true;
                }
                else nextHasReadLine = false;
            }
        }

        public void readLine()
        {
            if (!nextHasReadLine)
                read();
        }

        public bool hasNextDouble()
        {
            double dummy;
            return double.TryParse(it.Current, out dummy);
        }

        public bool hasNextInt()
        {
            int dummy;
            return int.TryParse(it.Current, out dummy);
        }

        public int nextInt()
        {
            try
            {
                return int.Parse(it.Current);
            }
            finally
            {
                if (!it.MoveNext() && es.Peek() >= 0)
                {
                    read();
                    nextHasReadLine = true;
                }
                else nextHasReadLine = false;
            }
        }

        public double nextDouble()
        {
            try
            {
                return double.Parse(it.Current);
            }
            finally
            {
                if (!it.MoveNext() && es.Peek() >= 0)
                {
                    read();
                    nextHasReadLine = true;
                }
                else nextHasReadLine = false;
            }
        }

        public void close()
        {
            es.Close();
        }


    }
}
