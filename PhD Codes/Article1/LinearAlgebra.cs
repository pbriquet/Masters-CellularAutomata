using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace Article1
{
    
    namespace PiLinearAlgebra
    {

        class Const
        {
            public const int maxdim = 5;
        }

        public class Real
        {
            public static double Min(double x, double y)
            {
                if (x <= y)
                    return x;
                else
                    return y;
            }
        }

        public class Angle
        {
            double theta;

            public Angle()
            {
                this.theta = 0.0;
            }
            public Angle(double t, string measure)
            {
                if (measure == "rad")
                    this.theta = t;
                else if (measure == "deg")
                    this.deg = t;
                else
                    throw new ArgumentException("Deg or Rad ?");
            }
            public double deg{
                get { return theta * 180.0 / Math.PI; }
                set { this.theta = value*Math.PI/180.0; }
            }
            public double rad
            {
                get { return theta; }
                set { this.theta = value; }
            }

            public static Angle operator +(Angle x, Angle y){
                Angle tmp = new Angle();
                tmp.rad = x.rad + y.rad;
                return tmp;
            }

            public static Angle operator -(Angle x, Angle y)
            {
                Angle tmp = new Angle();
                tmp.rad = x.rad - y.rad;
                return tmp;
            }
            public static Angle operator -(Angle x)
            {
                Angle tmp = new Angle();
                tmp.rad = -x.rad;
                return tmp;
            }
            public static Angle operator *(Angle x, Angle y)
            {
                Angle tmp = new Angle();
                tmp.rad = x.rad * y.rad;
                return tmp;
            }
            public static Angle operator *(double l, Angle x)
            {
                Angle tmp = new Angle();
                tmp.rad = l * x.rad;
                return tmp;
            }
            public static double operator *(Angle x, double l)
            {
                return l * x.rad;
            }
            public static Angle operator /(Angle x, Angle y)
            {
                Angle tmp = new Angle();
                tmp.rad = x.rad / y.rad;
                return tmp;
            }

            public double cos
            {
                get{ return Math.Cos(this.rad); }
                set{ this.theta = Math.Acos(value); }
            }
            public double sin
            {
                get { return Math.Sin(this.rad); }
                set { this.theta = Math.Asin(value); }
            }


        }

        public class Tensor
        {
            

            int[] dim = new int[3];
            double[,,] X = new double[Const.maxdim,Const.maxdim,Const.maxdim];

            public Tensor()
            {

            }

            public Tensor(int[] Dim, double[,,] value)
            {
                for (int i = 0; i < 3; i++)
                {
                    dim[i] = Dim[i];
                }
                for (int i = 0; i < this.dim[0]; i++)
                {
                    for (int j = 0; j < this.dim[1]; j++)
                    {
                        for (int k = 0; k < this.dim[2]; k++)
                        {
                            this.X[i,j,k] = value[i,j,k];
                        }
                    }
                }
            }

            public Tensor(int[] Dim)
            {
                for (int i = 0; i < 3; i++)
                {
                    this.dim[i] = Dim[i];
                }
            }

            public double this[int i, int j, int k]
            {
                get { return this.X[i,j,k]; }
                set { this.X[i,j,k] = value; }
            }

            public double[,,] m
            {
                get { return this.X; }
                set { this.X = value; }
            }

            public void setDim(int[] Dim)
            {
                this.dim = Dim;
            }

            public static bool CompareDim(Tensor x, Tensor y)
            {
                for (int i = 0; i < 3; i++)
                {
                    if (x.dim[i] != y.dim[i])
                        return false;
                }
                return true;
            }

        }

        public class Matrix
        {
            int[] dim = new int[2];
            double[,] X;
            public Matrix(){}

            public Matrix(int[] idim )
            {
                this.X = new double[idim[0], idim[1]];
            }

            public double this[int i, int j]
            {
                get { return this.X[i,j]; }
                set { this.X[i,j] = value; }
            }

            public double[,] v
            {
                get { return this.X; }
                set { this.X = value; }
            }

            
        }

        public class SqMatrix
        {
            int dim;
            double[,] X;

            public SqMatrix(int idim)
            {
                this.dim = idim;
                this.X = new double[idim, idim];
            }

            public double this[int i, int j]
            {
                get { return this.X[i, j]; }
                set { this.X[i, j] = value; }
            }

            public double[,] v
            {
                get { return this.X; }
                set { this.X = value; }
            }
            
        }

        public class SqMatrix3
        {
            const int dim = 3;
            double[,] X = new double[3, 3];

            public SqMatrix3(){}

            public double this[int i, int j]
            {
                get { return this.X[i, j]; }
                set { this.X[i, j] = value; }
            }

            public double[,] v
            {
                get { return this.X; }
                set { this.X = value; }
            }
            
        }

        public class SqMatrix2
        {
            const int dim = 2;
            double[,] X = new double[2, 2];
            
            public double this[int i, int j]
            {
                get { return this.X[i, j]; }
                set { this.X[i, j] = value; }
            }

            public double[,] v
            {
                get { return this.X; }
                set { this.X = value; }
            }

            public SqMatrix2(){}

            public static SqMatrix2 RotationMatrix(Angle angle)
            {
                SqMatrix2 M = new SqMatrix2();
                M[0, 0] = angle.cos;
                M[0, 1] = -angle.sin;
                M[1, 1] = M[0, 0];
                M[1, 0] = -M[0, 1];

                return M;
            }


        }

        public class Vector
        {
            int dim;
            double[] X = new double[Const.maxdim];

            public int Dim
            {
                get { return this.dim; }
                set { this.dim = value; }
            }
            public static bool CompareDim(Vector v1, Vector v2)
            {
                if (v1.dim != v2.dim)
                {
                    return false;
                }
                else
                    return true;
            }
            public Vector()
            {
                this.dim = Const.maxdim;
                this.v = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 };
            }

            public Vector( int idim )
            {
                this.dim = idim;
                this.v = new double[] {0.0, 0.0, 0.0, 0.0, 0.0};
            }

            public double this[int i]
            {
                get { return this.X[i]; }
                set { this.X[i] = value; }
            }

            public double[] v
            {
                get { return this.X; }
                set { this.X = value; }
            }

            public Vector(int idim, double[] X)
            {
                this.dim = idim;
                this.v = X;
            }

            public static Vector operator +(Vector v1, Vector v2)
            {

                if (CompareDim(v1, v2))
                {
                    Vector tmp = new Vector();
                    for (int i = 0; i < v1.dim; i++)
                    {
                        tmp[i] = v1[i] + v2[i];
                    }
                    return tmp;
                }
                else
                {
                    throw new ArgumentException("Different dims");
                }

                
            }
            
        }

        public class Vector2
        {
            const int dim = 2;

            double[] X = new double[dim];

            public Vector2()
            {
                this.v = new double[] { 0.0, 0.0 };
            }

            public Vector2(double[] x)
            {
                this.v = x;
            }

            public Vector2(double mod, Angle theta)
            {
                this.x = mod * theta.cos;
                this.y = mod * theta.sin;
            }

            public double this[int i]
            {
                get { if (i == 0 || i == 1) { return this.X[i]; } else { throw new ArgumentException("Out or range"); } }
                set { if (i == 0 || i == 1) { this.X[i] = value; } else { throw new ArgumentException("Out or range"); } }
            }

            public double[] v
            {
                get { return this.X; }
                set { this.X = value; }
            }

            public double x
            {
                get { return this.X[0]; }
                set { this.X[0] = value; }
            }
            public double y
            {
                get { return this.X[1]; }
                set { this.X[1] = value; }
            }



            public int Dim
            {
                get { return dim; }
            }

            public static Vector2 operator +(Vector2 v1, Vector2 v2)
            {
                Vector2 tmp = new Vector2();
                for (int i = 0; i < 2; i++)
                {
                    tmp[i] = v1[i] + v2[i];
                }
                return tmp;
            }

            public static Vector2 operator -(Vector2 v1, Vector2 v2)
            {
                Vector2 tmp = new Vector2();
                for (int i = 0; i < 2; i++)
                {
                    tmp[i] = v1[i] - v2[i];
                }
                return tmp;
            }

            public static double operator *(Vector2 v1, Vector2 v2)
            {
                double tmp = 0.0;
                for (int i = 0; i < 2; i++)
                {
                    tmp += v1[i] * v2[i];
                }
                return tmp;
            }

            public static Vector2 operator *(SqMatrix2 M, Vector2 x)
            {
                Vector2 tmp = new Vector2();
                for (int i = 0; i < x.Dim; i++)
                {
                    for (int j = 0; j < x.Dim; j++)
                    {
                        tmp[i] += M[i, j] * x[j];
                    }
                }
                return tmp;
            }

            public static Vector2 operator *(double l, Vector2 x)
            {
                Vector2 tmp = new Vector2();
                for (int i = 0; i < x.Dim; i++)
                {
                    tmp[i] = l * x[i];
                }
                return tmp;
            }

            public double norm()
            {
                return Math.Sqrt(this * this);
            }

            public void Rotate(Angle theta)
            {
                SqMatrix2 R = new SqMatrix2();
                Vector2 tmp = new Vector2();
                R = SqMatrix2.RotationMatrix(theta);
                tmp = R * this;
                this.v = tmp.v;
            }

            public void InverseRotate(Angle theta)
            {
                this.Rotate(-theta);
            }

        }

        public class Vector3 : Vector
        {

        }

        /*
        public class Tensor
        {
            int[] dim = new int[3];
            double[][][] X = new double[3][][];

            public static void SetElem(int[] i, Tensor q, double value)
            {
                q.X[i[0]][i[1]][i[2]] = value;
            }

            public Tensor(int[] Dim, double[][][] value)
            {
                for (int i = 0; i < 3; i++)
                {
                    dim[i] = Dim[i];
                }
                for (int i = 0; i < this.dim[0]; i++)
                {
                    for (int j = 0; j < this.dim[1]; j++)
                    {
                        for (int k = 0; k < this.dim[2]; k++)
                        {
                            this.X[i][j][k] = value[i][j][k];
                        }
                    }
                }
            }

            public Tensor(int[] Dim)
            {
                for (int i = 0; i < 3; i++)
                {
                    this.dim[i] = Dim[i];
                }
            }

            public Tensor()
            {

            }

            public static bool CompareDim(Tensor x, Tensor y)
            {
                for (int i = 0; i < 3; i++)
                {
                    if (x.dim[i] != y.dim[i])
                        return false;
                }
                return true;
            }

            public static Tensor Permutation;

            public void MakePermutationTensor(int Dim)
            {
                for (int i = 0; i < Dim; i++)
                {
                    for (int j = 0; j < Dim; j++)
                    {
                        for (int k = 0; k < Dim; k++)
                        {
                            Permutation.X[i][j][k] = (i - j) * (j - k) * (k - i) / 2.0;
                        }
                    }
                }
            }

            public static Tensor Equal(Tensor y)
            {
                Tensor tmp = new Tensor();
                for (int i = 0; i < 3; i++)
                {
                    tmp.dim[i] = y.dim[i];
                }
                for (int i = 0; i <= y.dim[2]; i++)
                {
                    for (int j = 0; j <= y.dim[1]; j++)
                    {
                        for (int k = 0; k <= y.dim[0]; k++)
                        {
                            tmp.X[i][j][k] = y.X[i][j][k];
                        }
                    }
                }
                return tmp;
            }

            public static Tensor operator +(Tensor x, Tensor y)
            {
                if (CompareDim(x, y))
                {
                    throw new ArgumentException("Different Dims");
                }

                Tensor tmp = new Tensor();

                for (int i = 0; i <= y.dim[2]; i++)
                {
                    for (int j = 0; j <= y.dim[1]; j++)
                    {
                        for (int k = 0; k <= y.dim[0]; k++)
                        {
                            tmp.X[i][j][k] = x.X[i][j][k] + y.X[i][j][k];
                        }
                    }
                }
                return tmp;
            }
            public static Tensor operator -(Tensor x, Tensor y)
            {
                Tensor tmp = new Tensor();

                for (int i = 0; i <= y.dim[2]; i++)
                {
                    for (int j = 0; j <= y.dim[1]; j++)
                    {
                        for (int k = 0; k <= y.dim[0]; k++)
                        {
                            tmp.X[i][j][k] = x.X[i][j][k] - y.X[i][j][k];
                        }
                    }
                }
                return tmp;
            }

        }




        public class Matrix : Tensor
        {
            public int[] dim = new int[2];
            double[][] X = new double[3][];
            public Matrix(int[] n)
            {
                this.dim[0] = n[0];
                this.dim[1] = n[1];
            }

            public Matrix(int n, int m)
            {
                this.dim[0] = n;
                this.dim[1] = m;
            }

            public Matrix(int n, int m, double[][] x)
            {
                this.dim[0] = n;
                this.dim[1] = m;
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < m; j++)
                    {
                        this.X[i][j] = x[i][j];
                    }
                }
            }

            public double this[int i, int j]
            {
                get { return this.X[i][j]; }
                set { this.X[i][j] = value; }
            }

            public double[][] m
            {
                get { return this.X; }
                set { this.X = value; }
            }

            public void setDim(int[] Dim)
            {
                this.dim = Dim;
            }

            public static Matrix operator +(Matrix A, Matrix B)
            {
                Matrix tmp = new Matrix(A.dim);
                for (int i = 0; i < A.dim[0]; i++)
                {
                    for (int j = 0; j < A.dim[1]; j++)
                    {
                        tmp[i, j] = A[i, j] + B[i, j];
                    }
                }
                return tmp;
            }

            public static Matrix operator -(Matrix A, Matrix B)
            {
                Matrix tmp = new Matrix(A.dim);
                for (int i = 0; i < A.dim[0]; i++)
                {
                    for (int j = 0; j < A.dim[1]; j++)
                    {
                        tmp[i, j] = A[i, j] - B[i, j];
                    }
                }
                return tmp;
            }

            public static Matrix operator *(Matrix A, Matrix B)
            {
                Matrix tmp = new Matrix(A.dim[0], B.dim[1]);
                for (int i = 0; i < A.dim[0]; i++)
                {
                    for (int j = 0; j < B.dim[1]; j++)
                    {
                        for (int k = 0; k < A.dim[1]; k++)
                        {
                            tmp[i, j] = A[i, k] + B[k, j];
                        }
                    }
                }
                return tmp;
            }

            public static Vector operator *(Matrix A, Vector x)
            {
                Vector tmp = new Vector(x.dim);
                for (int i = 0; i < A.dim[0]; i++)
                {
                    for (int j = 0; j < x.dim; j++)
                    {
                        tmp[i] = A[i, j] * x[j];
                    }
                }
                return tmp;
            }
            public static Matrix operator *(double l, Matrix A)
            {
                Matrix tmp = new Matrix(A.dim);
                for (int i = 0; i < A.dim[0]; i++)
                {
                    for (int j = 0; j < A.dim[1]; j++)
                    {
                        tmp[i, j] = A[i, j] * l;
                    }
                }
                return tmp;
            }

            public double[] GetCol(int i)
            {
                return this.X[i];
            }
            public double[] GetLin(int j)
            {
                Matrix tmp = this.Transpose();
                return tmp.X[j];
            }
            public void SetCol(int j, double[] col)
            {
                for (int i = 0; i < this.dim[0]; i++)
                {
                    this.X[i][j] = col[i];
                }
            }

            public void SetLin(int i, double[] lin)
            {
                for (int j = 0; j < this.dim[1]; j++)
                {
                    this.X[i][j] = lin[j];
                }
            }

            public static Matrix LineVectorsMatrix(Vector[] v)
            {
                Matrix tmp = new Matrix(v.Length, v[0].dim);
                for (int i = 0; i < tmp.dim[0]; i++)
                {
                    tmp.SetLin(i, v[i].v);
                }
                return tmp;
            }

            public static Matrix ColVectorsMatrix(Vector[] v)
            {
                Matrix tmp = new Matrix(v[0].dim, v.Length);
                for (int j = 0; j < tmp.dim[1]; j++)
                {
                    tmp.SetCol(j, v[j].v);
                }
                return tmp;
            }

            public Matrix Transpose()
            {
                Matrix tmp = new Matrix(this.dim[1], this.dim[0]);
                for (int i = 0; i < this.dim[0]; i++)
                {
                    for (int j = 0; j < this.dim[1]; j++)
                    {
                        tmp.X[j][i] = this.X[i][j];
                    }
                }

                return tmp;
            }

            public Matrix[] LU(Matrix A)
            {
                Matrix L = new Matrix(A.dim);
                Matrix U = new Matrix(A.dim);
                Matrix[] tmp = new Matrix[2];



                tmp[0].setDim(A.dim);
                tmp[1].setDim(A.dim);

                tmp[0].m = L.m;
                tmp[1].m = U.m;

                return tmp;

            }
            public static double Det(Matrix A)
            {
                double tmp = new double();
                int[] i = new int[A.dim[0]];
                for (int k = 0; k < i.Length; k++)
                {

                }
                return tmp;
            }

            public static Matrix Inv(Matrix A)
            {
                Matrix tmp = new Matrix(A.dim);

                return tmp;
            }

            public static Vector SolveLinearSystem(Matrix A, Vector b)
            {
                Vector x = new Vector(b.dim);

                return x;
            }
        }





        public class Vector : Tensor
        {

            double[] X = new double[3];
            public int dim { get; set; }

            public Vector(int Dim)
            {
                this.dim = Dim;
            }

            public Vector()
            {

            }

            public Vector(int Dim, double[] x)
            {
                this.dim = Dim;
                for (int i = 0; i < Dim; i++)
                {
                    this.X[i] = x[i];
                }
            }

            public static bool CompareDim(Vector x, Vector y)
            {
                if (x.dim == y.dim)
                    return true;
                else
                    return false;
            }


            public double[] v
            {
                get { return this.X; }
                set { this.X = value; }
            }

            public double x
            {
                get { return this.X[0]; }
                set { this.X[0] = value; }
            }

            public double y
            {
                get { return this.X[1]; }
                set { this.X[1] = value; }
            }

            public double z
            {
                get { return this.X[2]; }
                set { this.X[2] = value; }
            }

            public double this[int i]
            {
                get { return this.X[i]; }
                set { this.X[i] = value; }
            }

 
            public static Vector operator +(Vector x, Vector y)
            {
                if (!CompareDim(x, y))
                {
                    throw new ArgumentException("Different Dims");
                }

                Vector tmp = new Vector();
                tmp.dim = x.dim;

                for (int k = 0; k < y.dim; k++)
                {
                    tmp[k] = x[k] + y[k];
                }

                return tmp;
            }

 
            public static Vector operator -(Vector x, Vector y)
            {
                if (!CompareDim(x, y))
                {
                    throw new ArgumentException("Different Dims");
                }
                Vector tmp = new Vector(x.dim);


                for (int k = 0; k < y.dim; k++)
                {
                    tmp[k] = x[k] - y[k];
                }

                return tmp;
            }

            public static double operator *(Vector x, Vector y)
            {
                if (!CompareDim(x, y))
                {
                    throw new ArgumentException("Different Dims");
                }

                double tmp = new double();
                tmp = 0.0;
                for (int i = 0; i < x.dim; i++)
                {
                    tmp += x[i] + y[i];
                }
                return tmp;
            }

            public static Vector operator *(double l, Vector x)
            {
                Vector tmp = new Vector(x.dim);

                for (int i = 0; i < x.dim; i++)
                {
                    tmp[i] = l * x[i];
                }

                return tmp;
            }

 
            public static Vector operator ^(Vector x, Vector y)
            {
                if (!CompareDim(x, y))
                {
                    throw new ArgumentException("Different Dims");
                }

                Vector tmp = new Vector(x.dim, new double[] { 0.0, 0.0, 0.0 });

                for (int i = 1; i <= x.dim; i++)
                {
                    for (int j = 1; j <= x.dim; j++)
                    {
                        for (int k = 1; k <= x.dim; k++)
                        {
                            tmp.X[i - 1] += (((i - j) * (j - k) * (k - i) / 2.0)) * x.X[j - 1] * y.X[k - 1];

                        }
                    }
                }


                return tmp;
            }


            public static Vector GetHyperPlane(Vector x, Vector y, Vector z)
            {
                Vector tmp = new Vector(x.dim);

                return tmp;
            }

            public static double DistanceFromPlane(Vector x, Vector Plane)
            {
                double tmp = new double();

                return tmp;
            }

            public override string ToString()
            {
                string tmp;

                tmp = "{";
                for (int i = 0; i < this.dim; i++)
                {
                    tmp += Convert.ToString(this.X[i]);
                    if (i != this.dim - 1)
                        tmp += ",";
                }
                tmp += "}";

                return tmp;
            }
        }
        */ 
         

    }
    
}
