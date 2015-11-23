import java.util.Random;
import java.util.NoSuchElementException;

public class Jacobi {
    private static Matrix xApprox;

    public static Matrix jacobi_iter(Matrix xk, double tolerance, int maxIterations) {
        boolean toContinue = true;
        int iterations = 0;
        Matrix a = new Matrix(3, 3);
        a.set(0, 0, 1);
        a.set(1, 0, (double) 1 / 2);
        a.set(2, 0, (double) 1 / 3);
        a.set(0, 1, (double) 1 / 2);
        a.set(1, 1, 1);
        a.set(2, 1, (double) 1 / 4);
        a.set(0, 2, (double) 1 / 3);
        a.set(1, 2, (double) 1 / 4);
        a.set(2, 2, 1);
        Matrix c = new Matrix(3, 1);
        c.set(0, 0, 0.1);
        c.set(1, 0, 0.1);
        c.set(2, 0, 0.1);
        Matrix l = new Matrix(3, 3);
        Matrix d = new Matrix(3, 3);
        Matrix u = new Matrix(3 ,3);
        Matrix current = xk;
        for (int x = 0; x < 3; x++) {
            for (int y = 0; y < 3; y++) {
                if (x < y) {
                    l.set(y, x, a.get(y, x));
                } else if (x == y) {
                    d.set(y, x, a.get(y, x));
                } else if (x > y){
                    u.set(y, x, a.get(y, x));
                }
            }
        }
        System.out.print("X0:");
        current.print(8,8);
        while (toContinue && iterations < maxIterations) {
            d = inverseDiag(d);
            Matrix t = multiply(d, (l.plus(u)).uminus());
            Matrix previous = current;
            current = multiply(t, current).plus(c);
            iterations++;
            if (normInfinity(current.minus(previous)) < tolerance) {
                toContinue = false;
            }
        }
        if (iterations == maxIterations) {
            throw new NoSuchElementException("Cannot find an accurate solution.");
        }
        xApprox = current;
        System.out.print("XApprox");
        xApprox.print(8,8);
        System.out.println("Total Iterations: " + iterations);
        return xApprox;
    }

    public static Matrix calculateRandom() {
        Random randNum = new Random();
        Matrix approxTotal = new Matrix(3,1);
        for (int count = 0; count < 100; count++) {
            Matrix random = new Matrix(3, 1);
            for (int i = 0; i < 3; i++) {
                random.set(i, 0, randNum.nextDouble() - randNum.nextDouble());
            }
            approxTotal.plus(jacobi_iter(random, 0.00005, 100));
        }
        return approxTotal;
    }

    public static void iterativeError(Matrix approx, Matrix actual) {
        approx.timesEquals((double) 1 / 100);
        approx.minusEquals(actual);
        double norm = normInfinity(approx);
        System.out.println("Iterative Error: " + norm);
    }

    public static Matrix inverseDiag(Matrix a) {
        int x = 0;
        int y = 0;
        while (x < 3) {
            a.set(y, x, 1 / a.get(y, x));
            x++;
            y++;
        }
        return a;
    }

    public static double normInfinity(Matrix m) {
        double[] ans = new double[m.getRowDimension()];
        for(int i = 0; i < ans.length; i++) {
            for(int j = 0; j < m.getColumnDimension(); j++) {
                ans[i]+= Math.abs(m.get(i, j));
            }
        }
        double answer = ans[0];
        for(int i = 1; i < ans.length; i++) {
            if(ans[i] > answer){
                answer = ans[i];
            }
        }
        return answer;
    }

    public static Matrix multiply(Matrix a, Matrix b) {
        if (a.getColumnDimension() == b.getRowDimension()) {
            Matrix ans = new Matrix(a.getRowDimension(), b.getColumnDimension());
            for (int i = 0; i < ans.getRowDimension(); i++) {
                for (int j = 0; j < ans.getColumnDimension(); j++) {
                    for (int k = 0; k < a.getColumnDimension(); k++) {
                        ans.set(i, j, ans.get(i, j) + (a.get(i, k) * b.get(k, j)));
                    }
                }
            }
            return ans;
        } else {
            return null;
        }
    }
    
    public static void main(String[] args) {
        Matrix approxTotal = Jacobi.calculateRandom();
        Matrix actualTotal = new Matrix(3,1);
        actualTotal.set(0, 0, (double) 9 / 190);
        actualTotal.set(1, 0, (double) 28 / 475);
        actualTotal.set(2, 0, (double) 33 / 475);
        System.out.println("");
        iterativeError(approxTotal, actualTotal);
    }
}

