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
            Matrix b = new Matrix();
            b = Matrix.Parse("1\r\n" +
                "2\r\n" +
                "3\r\n");
            Console.WriteLine(b.ToString());
            m =Matrix.Parse("0.001 2.000 3.000\r\n" +
                "-1.000 3.712 4.623\r\n" + 
                "-2.000 1.072 5.643\r\n");
            Console.WriteLine(m.ToString());
            Matrix x = new Matrix(3,1);
            x=function(m, b);
            Console.WriteLine(x.ToString());
            Console.Read();
        }
        static Matrix function(Matrix m,Matrix b)
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
            Console.WriteLine(m.ToString());
            Console.WriteLine(b.ToString());
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
        }


    }
}
