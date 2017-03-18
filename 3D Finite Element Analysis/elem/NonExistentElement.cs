using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEA3D.elem
{
    public class NonExistentElement : Element
    {
        // Constructor for non existent element
        public NonExistentElement() : base(" ", 0, 0) 
        {
            //super (" ", 0, 0);
            // i.e return an error
        }
    }
}
