import Jama.*;
import java.util.Arrays;

public class HouseHolder {

    public static void main(String[] args) {
        double[][] test = {{1, 1, 1, 1}, {1, 2, 3, 4}, {1, 3, 6, 10}, {1, 4, 10, 20}};
        Matrix a = new Matrix(test);
        TheResult results = qr_fact_househ(a);
        System.out.println("The Q Matrix");
        results.getQ().print(2, 2);
        System.out.println("The R Matrix");
        results.getR().print(2, 2);
        System.out.println("ERROR: " + results.getError());

    }

    public static class TheResult {
        Matrix q;
        Matrix r;
        double error;

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

        public void setError(double error) {
            this.error = error;
        }
    }


    public static TheResult qr_fact_househ(Matrix a) {
        TheResult theResult = new TheResult();
        Matrix r = a;
        Matrix q = new Matrix(a.getRowDimension(), a.getColumnDimension());
        Matrix error = new Matrix(a.getRowDimension(), a.getColumnDimension());
        int count = 0;
        for (int i = 0; i < a.getColumnDimension(); i++) {
            Matrix subMatrix = r.getMatrix(i, a.getRowDimension() - 1, i, i);
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
                subMatrix = multiply(subMatrix, subMatrix.transpose());
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
                r = multiply(identityH, r);
                if (count == 0) {
                    q = identityH;
                } else {
                    q = multiply(q, identityH);
                }
                count++;
            }
        }
        theResult.setR(r);
        theResult.setQ(q);
        error = multiply(q, r);
        error.minusEquals(a);
        double theError = normInfinity(error);
        theResult.setError(theError);
        return theResult;


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
        if(a.getColumnDimension() == b.getRowDimension()) {
            Matrix ans = new Matrix(a.getRowDimension(), b.getColumnDimension());
            for (int i = 0; i < ans.getRowDimension(); i++) {
                for (int j = 0; j < ans.getColumnDimension(); j++) {
                    for (int k = 0; k < a.getColumnDimension(); k++) {
                        ans.set(i, j, ans.get(i, j) + (a.get(i, k) * b.get(k, j)));
                    }
                }
            }
            return ans;
        }
        else {
            return null;
        }
    }
}
