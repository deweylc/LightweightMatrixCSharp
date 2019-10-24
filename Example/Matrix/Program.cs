using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LightweightMatrixCSharp;
namespace _Matrix_
{
    class Program
    {

        static void Main(string[] args)
        {
            Matrix m = new Matrix();
            Matrix m_ = new Matrix();
            Matrix b = new Matrix();
            Matrix x = new Matrix();
            double w = 1.3;
            b = Matrix.Parse("7.2\r\n" +
                "8.3\r\n" +
                "4.2\r\n");
           
            m =Matrix.Parse("10 -1 -2\r\n" +
                "-1 10 -2\r\n" + 
                "-1 -1 5\r\n");
               
            x=Jacobi(m,b);
            Console.WriteLine(x.ToString());
            x = G_S(m, b);
            Console.WriteLine(x.ToString());
            x = SOR(m, b, w);
            Console.WriteLine(x.ToString());
            //     Console.WriteLine(m.ToString());
            //   Matrix x = new Matrix(3,1);
            //    x= Gauss(m, b);
            //   Console.WriteLine(x.ToString());
            //     LU(ref m);
            // m_ = m.Invert();//求逆

            Console.Read();
        }
        static Matrix Gauss(Matrix m,Matrix b)
        {
            Matrix x = new Matrix(m.rows, 1);
            double det = 1;
            for(int k= 0;k<m.cols-1;k++)
            {
                double p = 0;
                int flag = k;
                for(int i=k;i<m.cols;i++)//寻找最大列主元
                {
                    if (Math.Abs(m[i, k]) > p) 
                    {
                        p = Math.Abs(m[i, k]);
                        flag = i;
                    }
                }
                if(p == 0)//最大列主元=0，行列式为0,
                {
                    throw new Exception("error");
                }
                if(flag!=k)
                {
                    m = m.ChangeRow(flag, k);
                    b = b.ChangeRow(flag, k);
                    det = -1*det;
                }

                for (int i=k+1;i<m.cols;i++)
                {
                    double mik = m[i, k] / m[k, k];
                    for (int j = k; j < m.cols; j++)
                    {
                        m[i, j] = m[i, j] - mik * m[k, j];
                    }
                    b[i, 0] = b[i, 0] - mik * b[k, 0];
                }
               for(int i=0;i<m.rows;i++)
                {
                    for(int j=0;j<m.cols;j++)
                    {
                        m[i, j] = Math.Round(m[i, j], 10);
                    }
                    b[i, 0] = Math.Round(b[i, 0], 10);
                }
            }

            int n = x.rows - 1;

            for(int k=n;k>=0;k--)
            {
                if(k==n)
                {
                    x[k, 0] = Math.Round( b[n,0] / m[n,n],4);
                }
                else
                {
                    double s = 0;
                    for(int i=k+1;i<=n;i++)
                    {
                        s = s + m[k, i] * x[i, 0];
                    }
                    x[k, 0] =(b[k, 0] - s) / m[k, k];
                }
              
            }
            for(int i=0;i<x.rows;i++)
            {
                x[i, 0] = Math.Round(x[i, 0], 5);
            }
            return x;
        }//高斯，（列主元
        static void LU(ref Matrix m)
        {
            Console.WriteLine(m.ToString());
            m.L = Matrix.IdentityMatrix(m.rows, m.cols);
            m.U = Matrix.ZeroMatrix(m.rows, m.cols);
            int n = m.rows;
            for(int k=0;k<n;k++)
            {
                for(int i=k;i<n;i++)
                {
                    double s = 0;
                    for(int j=0;j<=k-1;j++)
                    {
                        s = s + m.L[k, j] * m.U[j, i];
                    }
                    m.U[k, i] =Math.Round( m[k, i] - s,4);
                }
                
                for(int i=k+1;i<n;i++)
                {
                    double sum = 0;
                    for (int j = 0; j <= k - 1; j++)
                    {
                        sum = sum + m.L[i, j] * m.U[j, k];
                    }
                    m.L[i, k] = Math.Round((m[i, k] - sum) / m.U[k, k], 4);
                }

                
                Console.WriteLine("L:");
                Console.WriteLine(m.L.ToString());
                Console.WriteLine("U:");
                Console.WriteLine(m.U.ToString());
            }
            Matrix a = new Matrix(m.rows, m.cols);
            a = m.L * m.U;
            Console.WriteLine(a.ToString());
        }//LU分解 A=LU
        static void LUD(ref Matrix m)
        {
            m.L = Matrix.ZeroMatrix(m.rows, m.cols);
            m.U = Matrix.ZeroMatrix(m.rows, m.cols);
            m.D = Matrix.ZeroMatrix(m.rows, m.cols);
            int n = m.rows;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        m.D[i, j] = m[i, j];
                    }
                    else if (i > j)
                    {
                        m.L[i, j] = -1 * m[i, j];
                    }
                    else if (i < j)
                    {
                        m.U[i, j] = -1 * m[i, j];
                    }
                }
            }
        }//LUD分解 A=D-L-U
        static Matrix Jacobi(Matrix _m,Matrix b)
        {
            Matrix m = new Matrix();
            m = _m;
            LUD(ref m);
            int n = m.rows;
            Matrix B = new Matrix();
            B = m.D.Invert() * (m.L + m.U);
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<n;j++)
                {
                    B[i, j] = Math.Round(B[i, j], 4);
                }
            }
            
            Matrix f = new Matrix(n, 1);
            f = m.D.Invert() * b;
           

            Matrix x = Matrix.ZeroMatrix(n, 1);
            Matrix xk =new Matrix(n, 1);
            int k = 1;
            while(true)
            {
                
                xk = B * x + f;
               
                x = xk;
                k++;
                if (k > 12)
                    break;
            }
            return xk;
           
          
        }//雅可比迭代
        static Matrix G_S(Matrix _m,Matrix b)
        {
            Matrix m = new Matrix();
            m = _m;
            LUD(ref m);
            int n = m.rows;
            Matrix B = new Matrix();
            B = (m.D-m.L).Invert() * m.U;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    B[i, j] = Math.Round(B[i, j], 4);
                }
            }

            Matrix f = new Matrix(n, 1);
            f = (m.D-m.L).Invert() * b;


            Matrix x = Matrix.ZeroMatrix(n, 1);
            Matrix xk = new Matrix(n, 1);
            int k = 1;
            while (true)
            {

                xk = B * x + f;

                x = xk;
                k++;
                if (k > 12)
                    break;
            }
            return xk;
        }//Gauss–Seidel method 高斯-赛德尔迭代
        static Matrix SOR(Matrix _m,Matrix b,double w)
        {
            Matrix m = new Matrix();
            m = _m;
            LUD(ref m);
            int n = m.rows;
            Matrix B = new Matrix();
            B = (m.D - w*m.L).Invert() * ((1-w)*m.D+w*m.U);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    B[i, j] = Math.Round(B[i, j], 4);
                }
            }

            Matrix f = new Matrix(n, 1);
            f = w*(m.D - w*m.L).Invert() * b;


            Matrix x = Matrix.ZeroMatrix(n, 1);
            Matrix xk = new Matrix(n, 1);
            int k = 1;
            while (true)
            {

                xk = B * x + f;

                x = xk;
                k++;
                if (k > 12)
                    break;
            }
            return xk;
        }//SOR超松弛迭代
    }
}
