using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Article1
{
    static class Fields
    {
        public static double Gx, Gy, R;
        public static double To;

        public static void SetThermalField(double Toi, double Ri, double Gxi, double Gyi)
        {
            Gx = Gxi; Gy = Gyi; R = Ri; To = Toi;
        }

        public static double SetGx(double Tright, double Tleft, double Lx)
        {
            return (Tright - Tleft) / Lx;
        }

        public static double SetGy(double Tup, double Tdown, double Ly)
        {
            return (Tup - Tdown) / Ly;
        }

        public static double SetR(double Tfinal, double Toi, double tmax)
        {
            return (Tfinal - Toi) / tmax;
        }

        public static double ThermalField( double t, double x, double y)
        {
            return To + R * t + Gx * x + Gy * y;
        }
    }
}
