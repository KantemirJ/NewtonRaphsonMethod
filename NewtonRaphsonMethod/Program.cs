using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

class Program
{
    static void Main()
    {
        // Initial guess
        var initialGuess = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(new double[] { -1, 9 });

        // Solve the system
        var solution = NewtonRaphson(initialGuess);

        // Print the solution
        Console.WriteLine("Solution: " + solution);
    }

    static MathNet.Numerics.LinearAlgebra.Vector<double> Equations(MathNet.Numerics.LinearAlgebra.Vector<double> x)
    {
        // Define the system of equations
        var f1 = Math.Pow(x[0], 3) + Math.Pow(x[1], 3) - 65;
        var f2 = Math.Pow(x[0], 2) * x[1] + x[0] * Math.Pow(x[1], 2)- 20;

        return MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(new double[] { f1, f2 });
    }

    static Matrix<double> Jacobian(MathNet.Numerics.LinearAlgebra.Vector<double> x)
    {
        // Calculate the Jacobian matrix
        var j11 = 3 * Math.Pow(x[0], 2);
        var j12 = 3 * Math.Pow(x[1], 2);
        var j21 = 2 * x[0] * x[1] + Math.Pow(x[1], 2);
        var j22 = Math.Pow(x[0], 2) + 2 * x[0] * x[1];

        return Matrix<double>.Build.DenseOfArray(new double[,] { { j11, j12 }, { j21, j22 } });
    }

    static MathNet.Numerics.LinearAlgebra.Vector<double> NewtonRaphson(MathNet.Numerics.LinearAlgebra.Vector<double> x0, double tol = 1e-6, int maxIter = 100)
    {
        for (int i = 0; i < maxIter; i++)
        {
            var f = Equations(x0);
            var J = Jacobian(x0);
            var delta_x = J.Solve(-f);
            x0 += delta_x;

            if (delta_x.L2Norm() < tol)
                break;
        }

        return x0;
    }
}
