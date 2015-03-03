import java.util.ArrayList;

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
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                System.out.print(system[i][j] + " ");
            }
            System.out.println(free[i]);
        }
        return res;
    }

    public static void main(String[] args) {

        double[][] matrix = {{8, 7, 3}, {-7, -4, -4}, {-6, 5, -4}};
        double[] free = {18, -11, -15};
        double[] res = gauss(matrix, free);
        for (int i = 0; i < N; ++i) {
            System.out.println(res[i]);
        }
    }
}
