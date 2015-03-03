public class Main {

    final static int N = 3;

    public static double[] gauss(double[][] system, double[] free) {
        for (int i = 0; i < N - 1; ++i) {
            for (int j = i + 1; j < N; ++j) {
                system[j][i] = -system[j][i] / system[i][i];
                for (int k = i + 1; k < N; ++k) {
                    system[j][k] += system[j][i] * system[i][k];
                }
                free[j] += system[j][i] * free[i];
            }
        }
        double[] res = new double[N];
        res[N - 1] = free[N - 1] / system[N - 1][N - 1];
        for (int i = N - 2; i >= 0; --i) {
            double temp = free[i];
            for (int j = i + 1; j < N; ++j) {
                temp -= res[j] * system[i][j];
            }
            res[i] = temp / system[i][i];
        }
        return res;
    }

    public static double[] jacobi(double[][] a, double[] f) {
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a.length;++j) {
                if (i != j) {
                    a[i][j] /= -a[i][i];
                }
            }
            f[i] /= a[i][i];
        }

        double[] currentX = new double[N];
        double[] nextX = new double[N];
        double cond = 10;
        final double EPS = 0.0000001;
        while (cond > EPS) {
            for (int i = 0; i < N; ++i) {
                double temp = f[i];
                for (int j = 0; j < N; ++j) {
                    if (i != j) {
                        temp += a[i][j] * currentX[j];
                    }
                }
                nextX[i] = temp;
            }
            double max = Double.MIN_VALUE;
            for (int i = 0; i < N; ++i) {
                if (max < Math.abs(currentX[i] - nextX[i])) {
                    max = Math.abs(currentX[i] - nextX[i]);
                }
            }
            cond = max;
            System.arraycopy(nextX, 0, currentX, 0, currentX.length);
        }
        return nextX;
    }

    public static void main(String[] args) {

        double[][] matrix = {{8, 7, 3}, {3, 5, 1}, {3, -2, 10}};
        double[] free = {10, 5, 4};
        double[] resJacobi = jacobi(matrix, free);
        double[] resGauss = gauss(matrix, free);
        for (int i = 0; i < N; ++i) {
            System.out.println(resJacobi[i] + " " + resGauss[i]);
        }
    }
}
