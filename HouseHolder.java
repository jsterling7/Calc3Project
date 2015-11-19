import Jama.*;
import java.util.Arrays;

public class HouseHolder {

    public static void main(String[] args) {
        double[][] test = {{1, 1, 1, 1}, {1, 2, 3, 4}, {1, 3, 6, 10}, {1, 4, 10, 20}};
        Matrix a = new Matrix(test);
        qr_fact_househ(a);
    }

    public static class TheResult {
        Matrix q;
        Matrix r;
        long error;

        public Matrix getQ() {
            return q;
        }

        public void setQ(Matrix q) {
            this.q = q;
        }

        public Matrix getR() {
            return r;
        }

        public void setR(Matrix r) {
            this.r = r;
        }

        public double getError() {
            return error;
        }

        public void setError() {
            this.error = error;
        }
    }


    public static TheResult qr_fact_househ(Matrix a) {
        TheResult theResult = new TheResult();
        Matrix r = a;
        Matrix q = new Matrix(a.getRowDimension(), a.getColumnDimension());
        int count = 0;
        for (int i = 0; i < a.getColumnDimension(); i++) {
            Matrix subMatrix = a.getMatrix(i, a.getRowDimension() - 1, i, i);
            Matrix subsubMatrix = subMatrix.getMatrix(1, subMatrix.getRowDimension() - 1, 0, 0);
            double[] checkForZeros = subsubMatrix.getColumnPackedCopy();
            boolean zeros = true;
            for (double item : checkForZeros) {
                if (item != 0) {
                    zeros = false;
                }
            }
            if (!zeros) {
                double firstElement = subMatrix.get(0,0);
                double[] elements = subMatrix.getColumnPackedCopy();
                double mag = 0;
                for (double item : elements) {
                    mag += item * item;
                }
                mag = Math.sqrt(mag);
                firstElement += mag;
                subMatrix.set(0, 0, firstElement);
                Matrix identity = Matrix.identity(subMatrix.getRowDimension(), subMatrix.getRowDimension());
                Matrix normSqr = subMatrix;
                subMatrix = subMatrix.times(subMatrix.transpose());
                subMatrix.timesEquals(2);
                double[] elements2 = normSqr.getColumnPackedCopy();
                double mag2 = 0;
                double one = 1;
                for (double item : elements2) {
                    mag2 += item * item;
                }
                mag2 = Math.sqrt(mag2);
                mag2 *= mag2;
                subMatrix = subMatrix.times(one / mag2);
                Matrix h = identity.minus(subMatrix);
                Matrix identityH = Matrix.identity(a.getRowDimension(), a.getColumnDimension());
                identityH.setMatrix(i, a.getRowDimension() - 1, i, a.getRowDimension() - 1, h);
                System.out.println("H" + (i + 1));
                identityH.print(2, 2);
                System.out.println("R");
                r.print(2, 2);
                r = identityH.times(r);
                if (count == 0) {
                    q = identityH;
                } else {
                    q = q.times(identityH);
                }
                count++;
            }
        }
        theResult.setR(r);
        theResult.setQ(q);
        System.out.println("R Final");
        r.print(2, 2);
        System.out.println("Q Final");
        q.print(2, 2);
        return theResult;

    }
}
