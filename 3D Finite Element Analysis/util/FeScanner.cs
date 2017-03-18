using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//using System.IO;

using FEA3D.model;

namespace FEA3D.util
{
    // FE data scanner. Delimiters: blank, =.
    public class FeScanner 
    {
        private Scanner es;


        // Construct FE data scanner.
        // fileIn - name of the file containing data.
        public FeScanner(String fileIn)
        {
            try
            {
                es = new Scanner(fileIn);
            }
            catch (Exception E)
            {
                UTIL.errorMsg("Input file not found: " + fileIn);
            }
            es.useDelimiters(new char[] { ' ', '=', '\t' });
        }

        // Returns  true if another token is available.
        public bool hasNext() { return es.hasNext(); }

        // Returns  true if double is next in input.
        public bool hasNextDouble() { return es.hasNextDouble(); }

        // Returns  true if int is next in input.
        public bool hasNextInt() { return es.hasNextInt(); }

        // Gives the next token from this scanner.
        public String next() { return es.next(); }

        // Gives the next double from this scanner.
        public double nextDouble() { return es.nextDouble(); }

        // Gives the next int from this scanner.
        public double nextInt() { return es.nextInt(); }

        // Reads the next integer.
        // Generates an error if next token is not integer.
        public int readInt()
        {
            if (!es.hasNextInt()) 
                UTIL.errorMsg("Expected integer. Instead: " + es.next());
            return es.nextInt();
        }

        // Reads the next double.
        // Generates an error if next token is not double.
        public double readDouble()
        {
            if (!es.hasNextDouble())
                UTIL.errorMsg("Expected double. Instead: " + es.next());
            return es.nextDouble();
        }

        // Advances the scanner past the current line.
        public void nextLine() { es.readLine(); }

        // Moves to line which follows a line with the word.
        public void moveAfterLineWithWord(String word)
        {

            while (es.hasNext())
            {
                String varname = es.next().ToLower();
                if (varname.Equals("#"))
                {
                    es.readLine();
                    continue;
                }
                if (varname.Equals(word))
                {
                    es.readLine();
                    return;
                }
            }
            UTIL.errorMsg("moveAfterLineWithWord cannot find: " + word);
        }

        // Method reads < nNumbers numbers > and places resulting
        // degrees of freedom in a List data structure.
        // Here numbers is a sequence of the type n1 n2 -n3 ...
        // where n2 -n3 means from n2 to n3 inclusive.
        // it - list iterator.
        // dir - direction (1,2,3).
        // nDim - problem dimension (2/3).
        // sValue - specified value.
        // returns - modified list iterator it.
        public List<Dof> readNumberList(List<Dof> it, int dir, int ndim, double sValue)
        {
            // number of items in the list
            int ndata = readInt();
            int i1, i2;
            i1 = i2 = readInt();
            for (int i = 1; i < ndata; i++)
            {
                i2 = readInt();
                if (i2 > 0 && i1 >= 0)
                {
                    if (i1 > 0)
                    {
                        it.Add(new Dof(ndim * (i1 - 1) + dir, sValue));
                    }
                    i1 = i2;
                }
                else if (i2 < 0)
                {
                    for (int j = i1; j <= (-i2); j++)
                    {
                        it.Add(new Dof(ndim * (j - 1) + dir, sValue));
                    }
                    i1 = 0;
                    i2 = 0;
                }
            }
            if (i2 > 0)
            {
                it.Add(new Dof(ndim * (i2 - 1) + dir, sValue));
            }
            return it;
        }

        // Closes this scanner.
        public void close() { es.close(); }

    }
}
