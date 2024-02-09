using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Article1.PiLinearAlgebra;

namespace Article1
{
    public class RunParameter
    {
        public double dtmax;
        public double tmax;
        public double RefreshTimeInterval;
        public int mx;
        public int my;
        public double Lx;
        public double Ly;

        public RunParameter(string stmax, string sdtmax, string rti, string mxi, string myi, string Lxi, string Lyi)
        {
            dtmax = Convert.ToDouble(sdtmax);
            tmax = Convert.ToDouble(stmax);
            RefreshTimeInterval = Convert.ToDouble(rti);
            mx = Convert.ToInt16(mxi);
            my = Convert.ToInt16(myi);
            Lx = Convert.ToDouble(Lxi);
            Ly = Convert.ToDouble(Lyi);
        }

    };


    public class Model
    {
        public RappazGandinCACell[,] CA;
        public RunParameter r;

        public double esi = 1.0;
        public double Lo;

        public class Al
        {
            public double esi = 1.0;

            public double pl = 1900.0;
            public double ps = 2000.0;
            public double Dl = 8.0e-9;
            //public double ve = 1.0e-4;
            public double kpart = 0.13;
            public double ml = -6.7;
            public double Co = 7.0;
            //public double Cd = 9.0;
            public double del_C = 1.0e-5;
            public double f2 = 4.0 * Math.Sqrt(3.0);
            public double f3 = 3.0 / 4.0;


            public double Cons;

            public Al()
            { }

            double AtoL(double A)
            {
                return Math.Sqrt(A / 2.0);
            }

            public void RefreshEsi(double dTdt, double T, double A, double ve, double dt)
            {
                double L = AtoL(A);
                double Cd = Co + (T - Alloys.Tliq) / ml;
                double dCdt = dTdt / ml;
                esi = esi + dt*( pl*(1.0-esi)/ps/(1-kpart)/Cd*dCdt -esi*3.0/L*ve + pl*Dl*f2/(ps*del_C*f3*L)*(Cd-Co)/Cd/(1.0 - kpart));
            }

            public double getA( RappazGandinCACell[,] CA, RunParameter r){
                double A = 0.0;
                for (int i = 0; i < r.my; i++)
                {
                    for (int j = 0; j < r.mx; j++)
                    {
                        if (CA[i, j].state != 0)
                        {
                            A += r.Lx / r.mx * r.Lx / r.mx;
                        }
                    }
                }
                return A;
            }

        }

        public Model(RunParameter c) { 

            
            r = c;
            CA = new RappazGandinCACell[r.my, r.mx];

            for( int i = 0 ; i < r.my ; i++){
                for( int j = 0 ; j< r.mx ; j++){
                    this.CA[i,j] = new RappazGandinCACell(r);
                    this.CA[i, j].setPos(j, i, r);
                }
            }

            CA[(int) (r.my/2.0), (int) (r.mx/2.0)].state = 1;
            CA[(int)(r.my / 2.0), (int)(r.mx / 2.0)].theta.deg = 30.0;
        }

        public void Initialize(int mx, int my)
        {
            
            
        }
        public void Run(double t)
        {
            for (int i = 0; i < r.my; i++)
            {
                for (int j = 0; j < r.mx; j++)
                {
                    this.CA[i, j].Growth(i, j, t, r, this.CA);
                }
            }
        }

        class RappazGandinCA
        {

           
        }

        class WangBeckermann1994
        {

        }
        public string Test()
        {
            string tmp;

            Vector2 v1 = new Vector2(new double[] {1.0, 2.0});
            Vector2 v2 = new Vector2(new double[] {3.0, 4.0});

            v1 = v1 + v2;
            return Convert.ToString(v1[0]);
        }

        public static bool CheckTimeInterval(double t, double dt, double ProfilesTimeInterval)
        {
            if ((t >= 0.0 && t < 0.0 + dt) || t * (1.0 / ProfilesTimeInterval) - (int)(t * (1.0 / ProfilesTimeInterval)) >= 0.0 && t * (1.0 / ProfilesTimeInterval) - (int)(t * (1.0 / ProfilesTimeInterval)) < dt * (1.0 / ProfilesTimeInterval))
                return true;
            else
            {
                return false;
            }
        }

        public void Test1()
        {
            double t;
            double dt = 0.001;

            for (t = 0.0; t < 10.0; t += dt)
            {

                if (CheckTimeInterval(t, dt, 1.0))
                {

                }
            }

        }

    }
}
