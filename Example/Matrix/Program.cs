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
            m=Matrix.Parse("1 0 0\r\n" +
                "1 1 0\r\n" + 
                "1 1 5\r\n");
            Console.WriteLine(m.ToString());
            Console.WriteLine(m.Det());
            Console.Read();
        }
    }
}
