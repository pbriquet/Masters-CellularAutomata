using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Article1.PiLinearAlgebra;

namespace Article1
{
    public class RappazGandinCACell
    {
        public static double dxCA;

        public double[] L = new double[4];
        public Angle theta = new Angle();
        public Vector2 dC = new Vector2();
        public double fs = new double();
        public int state = new int();
        public Vector2 pos = new Vector2();

        double dT = new double();

        public RappazGandinCACell(RunParameter r){
            dxCA = r.Lx / r.mx;
            L = new double[4];
            L[0] = 0.0; L[1] = 0.0; L[2] = 0.0; L[3] = 0.0;
            theta = new Angle();
            dC = new Vector2();
            fs = 0.0;
            state = 0;
            pos = new Vector2();

            dT = 0.0;
        }

        public void setPos(int j, int i, RunParameter r)
        {
            pos.x = dxCA * (j + 1 / 2);
            pos.y = r.Ly - dxCA * (i + 1 / 2);
        }

        public void Grow(RunParameter r, double t)
        {
            for (int k = 0; k < 4; k++)
            {

                Kinectics.ThermalParameter tmp = Kinectics.Param(new string[] { Convert.ToString(Alloys.Tliq - Fields.ThermalField(t, 0.0, 0.0)), Convert.ToString(Alloys.a2), Convert.ToString(Alloys.a3) }); ;
                this.L[k] += Convert.ToDouble(Kinectics.V(new object[] { tmp })) * r.dtmax;
            }
        }

        public void Growth(int i, int j, double t, RunParameter r, RappazGandinCACell[,] Mesh)
        {
            if (this.state > 0)
            {
                this.Grow(r, t);

                
                if (i != 0)
                    this.Capture(Mesh[i - 1, j]);
                if (i != r.my - 1)
                    this.Capture(Mesh[i + 1, j]);
                if (j != 0)
                    this.Capture(Mesh[i, j - 1]);
                if (j != r.mx - 1)
                    this.Capture(Mesh[i, j + 1]);
                /*
                if (i > 0 && j > 0)
                    this.Capture(Mesh[i - 1, j - 1]);
                if (i < r.my - 1 && j < r.mx - 1)
                    this.Capture(Mesh[i + 1, j + 1]);
                if (j > 0 && i < r.my - 1)
                    this.Capture(Mesh[i + 1, j - 1]);
                if (j < r.mx - 1 && i > 0)
                    this.Capture(Mesh[i - 1, j + 1]);
                */
                Inactivate(i, j, r, Mesh);

            }

        }

        public void Inactivate(int i, int j, RunParameter r, RappazGandinCACell[,] Mesh)
        {
            int counter = 0;
            if (i != 0) { if (Mesh[i - 1, j].state != 0) counter++; }
            else {counter++;}
            if (i != r.my - 1 ){if (Mesh[i + 1, j].state != 0) counter++;}
            else{counter++;}
            if (j != 0){if (Mesh[i, j - 1].state != 0)counter++;}
            else{counter++;}
            if (j != r.mx - 1){if (Mesh[i, j + 1].state != 0)counter++;}
            else{counter++;}
            
            /*
            if (i != 0 && j != 0) { if (Mesh[i - 1, j - 1].state != 0) counter++; }
            else { counter++; }
            if (i != 0 && j != r.mx - 1) { if (Mesh[i - 1, j + 1].state != 0) counter++; }
            else { counter++; }
            if (i != r.my - 1 && j != 0) { if (Mesh[i + 1, j - 1].state != 0) counter++; }
            else { counter++; }
            if (i != r.my - 1 && j != r.mx - 1) { if (Mesh[i + 1, j + 1].state != 0) counter++; }
            else { counter++; }
            */
            if (counter >= 4)
                Mesh[i, j].state = -Mesh[i, j].state;

        }
        public void Capture(RappazGandinCACell neigh)
        {
            if (neigh.state == 0)
            {
                Vector2 dr = new Vector2();

                dr = neigh.pos - this.pos;
                dr = dr - this.dC;

                dr.InverseRotate(this.theta);

                Vector2[] L = new Vector2[4];
                for (int k = 0; k < 4; k++)
                    L[k] = new Vector2();
                double[] l = new double[4];
                for (int k = 0; k < 4; k++)
                    l[k] = new double();
                int[] vecs = new int[2];
                for (int k = 0; k < 2; k++)
                    vecs[k] = 0;
                int nvec = 0;
                bool captured = false;

                for (int k = 0; k < 4; k++)
                {
                    L[k].x = ( (k+1)%2)*Math.Pow( -1.0, ( (int)( (k/2) )))*this.L[k] ;
                    L[k].y = ( k%2)*Math.Pow( -1.0, ( (int)( (k/2) )))*this.L[k];

                    l[k] = (L[k] * dr) / (L[k] * L[k]);

                    if (l[k] >= 0.0)
                    {
                        if (nvec < 2)
                        {
                            vecs[nvec] = k;
                            nvec++;
                        }
                        else
                        {
                            nvec = nvec;
                        }
                    }
                }
                if (nvec == 2)
                {
                    if (l[vecs[0]] + l[vecs[1]] < 1.0)
                        captured = true;
                }
                if (captured == true)
                {
                    double beta = new double();
                    double[] IL = new double[2];
                    double diffL1L2new = new double();
                    double diffL1L2 = new double();

                    Vector2[] newL = new Vector2[4];
                    Vector2[] newdC = new Vector2[2];

                    double alfa = new double();

                    diffL1L2 = (L[vecs[0]] - L[vecs[1]]).norm();

                    beta = ( (dr - L[vecs[1]] )*( L[vecs[0]] - L[vecs[1]] ) )/( Math.Pow(diffL1L2, 2.0) );
                    IL[0] = (1.0 - beta) * diffL1L2;
                    IL[1] = (beta) * diffL1L2;


                    diffL1L2new = (Real.Min(IL[0], Math.Sqrt(2.0) * dxCA) + Real.Min(IL[1], Math.Sqrt(2.0) * dxCA))/2.0;
                    alfa = diffL1L2new / diffL1L2;

                    for (int k = 0; k < 4; k++){
                        neigh.L[k] = alfa * this.L[k];
                        newL[k] = alfa*L[k];
                    }

                    neigh.state = this.state;
                    neigh.theta = this.theta;

                   
                    dr.Rotate(this.theta);
                    newdC[0] = (L[vecs[0]] - newL[vecs[0]]);
                    newdC[0].Rotate(this.theta);
                    newdC[1] = (L[vecs[1]] - newL[vecs[1]]);
                    newdC[1].Rotate(this.theta);

                    newdC[0] -= dr;
                    newdC[1] -= dr;

                    if ( newdC[0].norm() < newdC[1].norm())
                        neigh.dC = newdC[0];
                    else
                        neigh.dC = newdC[1];


                    
 
                }
            }

        }

    }
}
