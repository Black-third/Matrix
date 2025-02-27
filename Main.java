import java.util.*;

// Class representing a matrix with various operations
class Matrix {
    double[][] a; // 2D array to store matrix elements
    int n, m;     // Dimensions of the matrix (n = rows, m = columns)

    // Constructor to initialize a matrix with given dimensions
    Matrix(int x, int y) {
        n = x;
        m = y;
        a = new double[n][m]; // Allocate memory for the matrix
    }

    // Method to read matrix elements from a Scanner input
    void r(Scanner s) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                a[i][j] = s.nextDouble(); // Read each element
            }
        }
    }

    // Method to print the matrix
    void p() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                System.out.print(a[i][j] + " "); // Print each element
            }
            System.out.println(); // Move to the next line after each row
        }
    }

    // Method to add two matrices
    Matrix add(Matrix o) {
        if (o.n != n || o.m != m) {
            System.out.println("Matrices must have the same dimensions to add.");
            return new Matrix(0, 0); // Return an empty matrix if dimensions don't match
        }
        Matrix r = new Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                r.a[i][j] = a[i][j] + o.a[i][j]; // Add corresponding elements
            }
        }
        return r;
    }

    // Method to multiply two matrices
    Matrix mul(Matrix o) {
        if (o.n != m) {
            System.out.println("Number of columns of first matrix must equal number of rows of second matrix");
            return new Matrix(0, 0); // Return an empty matrix if dimensions don't match
        }
        Matrix r = new Matrix(n, o.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < o.m; j++) {
                for (int k = 0; k < m; k++) {
                    r.a[i][j] += a[i][k] * o.a[k][j]; // Perform matrix multiplication
                }
            }
        }
        return r;
    }

    // Method to transpose the matrix
    Matrix t() {
        Matrix r = new Matrix(m, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                r.a[j][i] = a[i][j]; // Swap rows and columns
            }
        }
        return r;
    }

    // Method to raise the matrix to a given power
    Matrix pow(int p) {
        if (m != n) {
            System.out.println("Matrix must be square to be exponentiated");
            return new Matrix(0, 0); // Return an empty matrix if not square
        }
        Matrix r = I(n); // Initialize result as identity matrix
        Matrix b = this;
        while (p > 0) {
            if (p % 2 == 1) r = r.mul(b); // Multiply if power is odd
            b = b.mul(b); // Square the matrix
            p /= 2; // Halve the power
        }
        return r;
    }

    // Method to create an identity matrix of given size
    Matrix I(int s) {
        Matrix r = new Matrix(s, s);
        for (int i = 0; i < s; i++) r.a[i][i] = 1; // Set diagonal elements to 1
        return r;
    }

    // Method to compute the inverse of the matrix
    Matrix inv() {
        if (m != n) {
            System.out.println("Matrix must be square to have an inverse");
            return new Matrix(0, 0); // Return an empty matrix if not square
        }
        Matrix r = I(n); // Initialize result as identity matrix
        for (int i = 0; i < n; i++) {
            int k = i;
            for (int j = i + 1; j < n; j++) if (Math.abs(a[j][i]) > Math.abs(a[k][i])) k = j;
            double[] t = a[i];
            a[i] = a[k];
            a[k] = t;
            t = r.a[i];
            r.a[i] = r.a[k];
            r.a[k] = t;
            double d = a[i][i];
            for (int j = 0; j < n; j++) {
                a[i][j] /= d;
                r.a[i][j] /= d;
            }
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double f = a[j][i];
                    for (int k1 = 0; k1 < n; k1++) {
                        a[j][k1] -= f * a[i][k1];
                        r.a[j][k1] -= f * r.a[i][k1];
                    }
                }
            }
        }
        return r;
    }

    // Method to solve a system of linear equations
    void solve(double[] b) {
        if (m != n || b.length != n) {
            System.out.println("Matrix must be square and the number of equations should be the same as the number of rows of the matrix");
            return;
        }
        double[][] aug = new double[n][n + 1];
        // Ensuring solve is correct
        Matrix a1 = new Matrix(m, n);
        Matrix a2 = new Matrix(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                a1.a[i][j] = a[j][i];
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                a2.a[i][j] = a1.a[j][i];
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                a[i][j] = a2.a[i][j];
            }
        }
        for (int i = 0; i < n; i++) {
            System.arraycopy(a[i], 0, aug[i], 0, n);
            aug[i][n] = b[i];
        }
        for (int i = 0; i < n; i++) {
            int k = i;
            for (int j = i + 1; j < n; j++) if (Math.abs(aug[j][i]) > Math.abs(aug[k][i])) k = j;
            double[] t = aug[i];
            aug[i] = aug[k];
            aug[k] = t;
            double d = aug[i][i];
            for (int j = i; j <= n; j++) aug[i][j] /= d;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    double f = aug[j][i];
                    for (int k1 = i; k1 <= n; k1++) aug[j][k1] -= f * aug[i][k1];
                }
            }
        }
        for (int i = 0; i < n; i++) {
            System.out.print(aug[i][n] + " "); // Print the solution
        }
        System.out.println();
    }

    // Method to compute the trace of the matrix
    double tr() {
        if (n != m) {
            System.out.println("Matrix must be square to have a trace");
            return 0;
        }
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += a[i][i]; // Sum the diagonal elements
        }
        return sum;
    }

    // Method to compute the determinant of the matrix
    double det() {
        if (m != n) {
            System.out.println("Matrix must be square to have a determinant");
            return -1;
        }
        double d = 1;
        Matrix c = new Matrix(n, n);
        for (int i = 0; i < n; i++) c.a[i] = Arrays.copyOf(a[i], n);
        for (int i = 0; i < n; i++) {
            int k = i;
            for (int j = i + 1; j < n; j++) if (Math.abs(c.a[j][i]) > Math.abs(c.a[k][i])) k = j;
            if (i != k) {
                double[] t = c.a[i];
                c.a[i] = c.a[k];
                c.a[k] = t;
                d = -d;
            }
            d *= c.a[i][i];
            for (int j = i + 1; j < n; j++) {
                double f = c.a[j][i] / c.a[i][i];
                for (int k1 = i; k1 < n; k1++) c.a[j][k1] -= f * c.a[i][k1];
            }
        }
        return d;
    }

    // Method to compute the rank of the matrix
    int rank() {
        int rank = 0;
        Matrix temp = new Matrix(n, m);
        for (int i = 0; i < n; i++) temp.a[i] = Arrays.copyOf(a[i], m);
        for (int i = 0; i < m; i++) {
            int k = rank;
            while (k < n && temp.a[k][i] == 0) k++;
            if (k == n) continue;
            double[] t = temp.a[k];
            temp.a[k] = temp.a[rank];
            temp.a[rank] = t;
            for (int j = rank + 1; j < n; j++) {
                double f = temp.a[j][i] / temp.a[rank][i];
                for (int l = i; l < m; l++) temp.a[j][l] -= f * temp.a[rank][l];
            }
            rank++;
        }
        return rank;
    }
}

// Class to compute Fibonacci numbers using matrix exponentiation
class FibonacciMatrix {
    private static final int MOD = 1000000007; // Modulo to prevent overflow

    // Method to multiply two 2x2 matrices
    private static long[][] multiply(long[][] A, long[][] B) {
        long[][] result = new long[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                result[i][j] = 0;
                for (int k = 0; k < 2; k++) {
                    result[i][j] = (result[i][j] + A[i][k] * B[k][j] % MOD) % MOD;
                }
            }
        }
        return result;
    }

    // Method to raise a 2x2 matrix to a given power
    private static long[][] matrixPower(long[][] A, int n) {
        long[][] result = {{1, 0}, {0, 1}}; // Initialize result as identity matrix
        while (n > 0) {
            if (n % 2 == 1) {
                result = multiply(result, A); // Multiply if power is odd
            }
            A = multiply(A, A); // Square the matrix
            n /= 2; // Halve the power
        }
        return result;
    }

    // Method to compute the nth Fibonacci number using matrix exponentiation
    public static long fibonacci(int n) {
        if (n == 0) return 0;
        if (n == 1) return 1;

        long[][] F = {{1, 1}, {1, 0}}; // Fibonacci matrix
        long[][] result = matrixPower(F, n - 1); // Raise matrix to (n-1)th power
        return result[0][0] % MOD; // Return the Fibonacci number
    }
}

// Main class to test the matrix operations
public class Main {
    public static void main(String[] args) {
        Scanner s = new Scanner(System.in);
        while (true) {
            System.out.println("Choose an operation:");
            System.out.println("1. Add Matrices");
            System.out.println("2. Multiply Matrices");
            System.out.println("3. Transpose Matrix");
            System.out.println("4. Matrix Power");
            System.out.println("5. Inverse Matrix");
            System.out.println("6. Rank of Matrix");
            System.out.println("7. Solve System of Equations");
            System.out.println("8. Matrix Determinant");
            System.out.println("9. Matrix Trace");
            System.out.println("10. Fibonacci Number using Matrix Exponentiation");
            System.out.println("11. Exit");
            int choice = s.nextInt();

            if (choice == 11) break; // Exit the loop if user chooses to exit

            System.out.print("Enter dimensions of matrix: ");
            int n = s.nextInt(), m = s.nextInt();
            Matrix a = new Matrix(n, m);
            a.r(s); // Read the matrix

            if (choice == 1) {
                System.out.print("Enter dimensions of second matrix: ");
                int p = s.nextInt(), q = s.nextInt();
                Matrix b = new Matrix(p, q);
                b.r(s);
                Matrix result = a.add(b); // Add matrices
                result.p(); // Print the result
            } else if (choice == 2) {
                System.out.print("Enter dimensions of second matrix: ");
                int p = s.nextInt(), q = s.nextInt();
                Matrix b = new Matrix(p, q);
                b.r(s);
                Matrix result = a.mul(b); // Multiply matrices
                result.p(); // Print the result
            } else if (choice == 3) {
                a.t().p(); // Transpose the matrix and print
            } else if (choice == 4) {
                System.out.print("Enter power: ");
                int pow = s.nextInt();
                Matrix result = a.pow(pow); // Raise matrix to power
                result.p(); // Print the result
            } else if (choice == 5) {
                Matrix result = a.inv(); // Compute inverse
                result.p(); // Print the result
            } else if (choice == 6) {
                System.out.println("Rank: " + a.rank()); // Compute and print rank
            } else if (choice == 7) {
                System.out.print("Enter the right-hand side equation values: ");
                double[] b = new double[n];
                for (int i = 0; i < n; i++) b[i] = s.nextDouble();
                a.solve(b); // Solve system of equations
            } else if (choice == 8) {
                double x = a.det(); // Compute determinant
                if (n == m) {
                    System.out.println("Determinant: " + x); // Print determinant
                }
            } else if (choice == 9) {
                double x = a.tr(); // Compute trace
                if (n == m) {
                    System.out.println("Trace: " + x); // Print trace
                }
            } else if (choice == 10) {
                System.out.println("Enter n for Fibonacci calculation:");
                int x = s.nextInt();
                if (x >= 0) System.out.println("Fibonacci number F(" + x + ") = " + FibonacciMatrix.fibonacci(x));
                else {
                    if (x % 2 == 0) System.out.println("Fibonacci number F(" + x + ") = " + FibonacciMatrix.fibonacci(-x));
                    else System.out.println("Fibonacci number F(" + x + ") = " + -FibonacciMatrix.fibonacci(-x));
                }
            } else {
                System.out.println("Invalid choice. Please try again.");
            }
        }
        s.close(); // Close the scanner
    }
}