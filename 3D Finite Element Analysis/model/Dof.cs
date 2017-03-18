using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEA3D.model
{
    public class Dof
    {
        public int dofNum;
        public double value;

        public Dof(int dofNum, double value)
        {
            this.dofNum = dofNum;
            this.value = value;
        }
    }
}
