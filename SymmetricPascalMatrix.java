import Jama.Matrix;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;


public class SymmetricPascalMatrix {
    public static void main(String[] args) throws IOException {
        System.out.println("-------------------------------------------------");
        System.out.println("THE LU DECOMPOSITION");
        System.out.println("-------------------------------------------------");
        luRun();
        System.out.println("-------------------------------------------------");
        System.out.println("THE HOUSEHOLDER QR FACTORIZATION");
        System.out.println("-------------------------------------------------");
        HouseHolderRun();
        System.out.println("-------------------------------------------------");
        System.out.println("THE GIVENS QR FACTORIZATION");
        System.out.println("-------------------------------------------------");
        givensRun();
        System.out.println("-------------------------------------------------");
        System.out.println("Solving Ax = b");
        System.out.println("-------------------------------------------------");
        Matrix augMatrix = readDat("Enter the .dat file name for the augumented matrix to solve using LU and QR: ");
        double[][] aug = augMatrix.getArrayCopy();
        double[][] bArray = new double[aug.length][1];
        double[][] aArray = new double[aug.length][aug[0].length - 1];
        for (int i = 0; i < aug.length; i++) {
            for (int j = 0; j < aug[i].length; j++) {
                if (j == aug[i].length - 1) {
                    bArray[i][0] = aug[i][j];
                } else {
                    aArray[i][j] = aug[i][j];
                }
            }
        }
        Matrix error = new Matrix(0, 0);
        Matrix error2 = new Matrix(0, 0);
        Matrix error3 = new Matrix(0, 0);
        Matrix b = new Matrix(bArray);
        Matrix a = new Matrix (aArray);
        Matrix b2 = new Matrix(b.getArrayCopy());
        Matrix a2 = new Matrix (a.getArrayCopy());
        TheResult resultsLU = solve_lu_b(a, b);
        TheResult resultsQR = solve_qr_b(a2, b2);
        System.out.println("LU\n");
        System.out.println("Xsol =");
        resultsLU.getX().print(2, 3);
        error = multiply(a, resultsLU.getX());
        error.minusEquals(b);
        System.out.println("ERROR: ||AXsol - b||∞ = " + normInfinity(error) + "\n");
        System.out.println("HouseHolder\n");
        System.out.println("Xsol =");
        resultsQR.getXHouse().print(2, 3);
        error2 = multiply(a, resultsQR.getXHouse());
        error2.minusEquals(b);
        System.out.println("ERROR: ||PXsol - b||∞ = " + normInfinity(error2) + "\n");
        System.out.println("Givens\n");
        System.out.println("Xsol =");
        resultsQR.getXGivens().print(2, 3);
        error3 = multiply(a, resultsQR.getXGivens());
        error3.minusEquals(b);
        System.out.println("ERROR: ||PXsol - b||∞ = " + normInfinity(error3) + "\n");
        pxb();
    }

    // ------------------------------------------------------------
    // LU
    // ------------------------------------------------------------
    public static Matrix[] lu_fact(Matrix m) {
		int n = m.getRowDimension();
		double[][] l = new double[n][n];// lower
		double[][] u = m.getArrayCopy();// upper
		int r = 1;// row
		int c = 0;// column
		for (int k = 0; k < n; k++) {// set the ones of L
			l[k][k] = 1;
		}
		while (c < n) {
			while (r < n) {
				if (u[r][c] != 0) {
					double x = -u[r][c] / u[c][c];
					l[r][c] = -x;
					int cc = c;
					while (cc < n) {// goes across changing all values of the
									// modified row
						u[r][cc] = (u[c][cc] * x) + u[r][cc];
						cc++;
					}
				}
				r++;
			}
			c++;
			r = c + 1;
		}
		Matrix la = new Matrix(l);
		Matrix ua = new Matrix(u);
		Matrix[] ans = { la, ua };
		return ans;
	}

	public static void luRun() throws IOException {
        Matrix m = readDat("Enter the .dat file name for the LU FACTORIZATION: ");
        Matrix[] ans = lu_fact(m);
        System.out.println("L = ");
		ans[0].print(2, 2);
		System.out.println("U = ");
		ans[1].print(2, 2);
		System.out.println("ERROR: ||LU - A||∞ = "
				+ normInfinity(multiply(ans[0], ans[1]).minus(m)));
	}
    // ------------------------------------------------------------
    // HouseHolder
    // ------------------------------------------------------------

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

    public static void HouseHolderRun() throws IOException {
        Matrix a = readDat("Enter the .dat file name for the HouseHolder QR Decomposition: ");
        TheResult results = qr_fact_househ(a);
            System.out.println("Q =");
            results.getQ().print(2, 2);
            System.out.println("R =");
            results.getR().print(2, 2);
            System.out.println("ERROR: ||QR−A||∞ = " + results.getError());
    }

    // ------------------------------------------------------------
    // Givens
    // ------------------------------------------------------------

    public static TheResult qr_fact_givens(Matrix a) {
        TheResult theResult = new TheResult();
        double[][] copy = a.getArrayCopy();
        Matrix q = Matrix.identity(a.getRowDimension(),a.getColumnDimension());
        Matrix r = a;
        for (int currentColumn = 0; currentColumn < a.getColumnDimension(); currentColumn++) {
            for (int currentRow = currentColumn; currentRow < a.getRowDimension() - 1; currentRow++) {
                Matrix newMatrix = new Matrix(a.getRowDimension(),a.getColumnDimension());
                if (copy[currentRow + 1][currentColumn] != 0.0) {
                    newMatrix.set(currentColumn,currentColumn,copy[currentColumn][currentColumn]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2)
                            + Math.pow(copy[currentRow + 1][currentColumn], 2))));
                    newMatrix.set(currentRow + 1, currentColumn, -(copy[currentRow + 1][currentColumn]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2)
                            + Math.pow(copy[currentRow + 1][currentColumn], 2)))));
                    newMatrix.set(currentColumn, currentRow + 1, copy[currentRow + 1][currentColumn]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2)
                            + Math.pow(copy[currentRow + 1][currentColumn], 2))));
                    newMatrix.set(currentRow + 1, currentRow + 1, copy[currentColumn][currentColumn]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2)
                            + Math.pow(copy[currentRow + 1][currentColumn], 2))));
                    for (int x = 0; x < a.getColumnDimension(); x++) {
                        for (int y = 0; y < a.getRowDimension(); y++) {
                            if ((x != currentColumn && y != currentColumn) &&
                                    (x != currentColumn && y != currentRow + 1) &&
                                    (x != currentRow + 1 && y != currentColumn) &&
                                    (x != currentRow + 1 && y!= currentRow + 1))
                            {
                                if (x==y) {
                                    newMatrix.set(y, x, 1.0);
                                } else {
                                    newMatrix.set(y, x, 0.0);
                                }
                            }
                        }
                    }
                    q = multiply (q, newMatrix.transpose());
                    r = multiply(newMatrix, r);
                    copy = r.getArrayCopy();
                }
            }
        }
        theResult.setQ(q);
        theResult.setR(r);
        Matrix errorMatrix = multiply(q, r);
        errorMatrix.minusEquals(a);
        double error = normInfinity(errorMatrix);
        theResult.setError(error);
        return theResult;
    }


    public static void givensRun() throws IOException {
        Matrix test = readDat("Enter the .dat file name for the Givens QR Decomposition: ");
        TheResult answers = qr_fact_givens(test);
        System.out.println("Q =");
        answers.getQ().print(2, 2);
        System.out.println("R =");
        answers.getR().print(2, 2);
        System.out.println("ERROR: ||QR−A||∞ = " + answers.getError());
    }

    // ------------------------------------------------------------
    // Solve Ax = b
    // ------------------------------------------------------------

    public static TheResult solve_lu_b(Matrix m, Matrix b) {
        TheResult theResult = new TheResult();
		Matrix[] ans = lu_fact(m);
		Matrix l = ans[0];
		Matrix u = ans[1];
		Matrix y = new Matrix(b.getArray());
		Matrix x = new Matrix(b.getRowDimension(), 1);
		for (int i = 0; i < y.getRowDimension(); i++) {
			double xx = b.get(i, 0);
			for (int j = 0; j < i; j++) {
				xx -= l.get(i, j) * y.get(j, 0);
			}
			y.set(i, 0, xx);
		}
		for (int i = x.getRowDimension() - 1; i >= 0; i--) {
			double xx = y.get(i, 0);
			for (int j = m.getColumnDimension() - 1; j > i; j--) {
				xx -= u.get(i, j) * x.get(j, 0);
			}
			x.set(i, 0, xx / u.get(i, i));
		}
		theResult.setX(x);
        theResult.setError(normInfinity(multiply(ans[0], ans[1]).minus(m)));
        return theResult;
	}

    public static TheResult solve_qr_b(Matrix a, Matrix b) {
        TheResult theResult = new TheResult();
        TheResult answersHouseHolder = qr_fact_househ(a);
        TheResult answersGivens = qr_fact_givens(a);

        Matrix qHouse = answersHouseHolder.getQ();
        Matrix rHouse = answersHouseHolder.getR();
        Matrix qGivens = answersGivens.getQ();
        Matrix rGivens = answersGivens.getR();

        theResult.setErrorHouse(answersHouseHolder.getError());
        theResult.setErrorGivens(answersGivens.getError());
        // For HouseHolder
        Matrix y = multiply(qHouse.transpose(), b);
        Matrix x = new Matrix(b.getRowDimension(), 1);
        for (int i = b.getRowDimension() - 1; i >= 0; i--) {
            double temp = y.get(i,0);
            for (int j = i + 1; j < b.getRowDimension(); j++) {
                 temp -= (rHouse.get(i, j) * x.get(j, 0));
            }
            x.set(i, 0, temp/ rHouse.get(i,i));
        }
        theResult.setXHouse(x);

        // For Givens
        y = multiply(qGivens.transpose(), b);
        Matrix x2 = new Matrix(b.getRowDimension(), 1);
            double temp = y.get(i,0);
            for (int j = i + 1; j < b.getRowDimension(); j++) {
                 temp -= (rGivens.get(i, j) * x.get(j, 0));
            }
            x2.set(i, 0, temp/ rGivens.get(i,i));
        }
        theResult.setXGivens(x2);
        return theResult;
    }

    // ------------------------------------------------------------
    // Solving Px = b
    // ------------------------------------------------------------

    public static void pxb() {
        Matrix p = new Matrix(0, 0);
        Matrix b = new Matrix(0, 0);
        Matrix error = new Matrix(0, 0);
        Matrix error2 = new Matrix(0, 0);
        TheResult theResult = new TheResult();
        Matrix[] solutions = new Matrix[2];
        // LU
        System.out.println("-------------------------------------------------");
        System.out.println("USING LU TO SOLVE Px = b");
        System.out.println("-------------------------------------------------");
        for (int i = 2; i <= 12; i++) {
            p = makeP(i);
            Matrix pCopy = new Matrix(p.getArrayCopy());
            b = makeB(i);
            theResult = solve_lu_b(p, b);
            error = multiply(pCopy, theResult.getX());
            error.minusEquals(b);
            System.out.println("------------------");
            System.out.println("n = " + i);
            System.out.println("------------------\n");
            System.out.println("Xsol =");
            theResult.getX().print(2, 3);
            System.out.println("ERROR: ||LU - P||∞ = " + theResult.getError());
            System.out.println("ERROR: ||PXsol - b||∞ = " + normInfinity(error) + "\n");
        }
        System.out.println("-------------------------------------------------");
        System.out.println("USING QR TO SOLVE Px = b");
        System.out.println("-------------------------------------------------");
        // HouseHolder
        for (int i = 2; i <= 12; i++) {
            p = makeP(i);
            Matrix pCopy = new Matrix(p.getArrayCopy());
            b = makeB(i);
            theResult = solve_qr_b(p, b);
            error = multiply(pCopy, theResult.getXHouse());
            error.minusEquals(b);
            error2 = multiply(pCopy, theResult.getXGivens());
            error.minusEquals(b);
            System.out.println("------------------");
            System.out.println("n = " + i);
            System.out.println("------------------\n");
            System.out.println("HOUSEHOLDER");
            System.out.println("-----------\n");
            System.out.println("Xsol =");
            theResult.getXHouse().print(2, 3);
            System.out.println("ERROR: ||QR - P||∞ = " + theResult.getErrorHouse());
            System.out.println("ERROR: ||PXsol - b||∞ = " + normInfinity(error) + "\n");
            System.out.println("GIVENS");
            System.out.println("------\n");
            System.out.println("Xsol =");
            theResult.getXGivens().print(2, 3);
            System.out.println("ERROR: ||QR - P||∞ = " + theResult.getErrorGivens());
            System.out.println("ERROR: ||PXsol - b||∞ = " + normInfinity(error2) + "\n");
        }
    }

    // ------------------------------------------------------------
    // Helper Methods/Classes
    // ------------------------------------------------------------


    public static Matrix makeP(int n) {
    		Matrix ans = new Matrix(n,n);
    		for (int i = 0; i < n; i++) {
    			ans.set(i, 0, 1);
    			ans.set(0, i, 1);
    		}
    		for (int i = 1; i < n; i++) {
    			for (int j = 1; j < n; j++) {
    				double k = ans.get(i-1,j) + ans.get(i, j-1);
    				ans.set(i, j, k);
    			}
    		}
    		return ans;
    }

    public static Matrix makeB(int n) {
        Matrix ans = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double k = (double) 1 / (i + 1);
            ans.set(i, 0, k);
        }
        return ans;
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

    public static Matrix readDat(String message) throws IOException{
        Scanner sc = new Scanner(System.in);
        System.out.print(message);
        String fileName = sc.next();
        File theFile = new File(fileName);
        Scanner sc2 = new Scanner(theFile);
        Scanner sc4 = new Scanner(theFile);
        Scanner sc3;
        int i = 0;
        int j = 0;
        int rows = 0;
        int columns = 0;
        while (sc2.hasNextLine()) {
            rows++;
            String line = sc2.nextLine();
            if (columns == 0) {
                sc3 = new Scanner(line);
                while (sc3.hasNextDouble()) {
                    sc3.nextDouble();
                    columns++;
                }
            }
        }
        double[][] data = new double[rows][columns];
        while (sc4.hasNextLine()) {
            String line = sc4.nextLine();
                sc3 = new Scanner(line);
                while (sc3.hasNextDouble()) {
                    data[i][j] = sc3.nextDouble();
                    j++;
                }
                j = 0;
                i++;
        }
        Matrix a = new Matrix(data);
        return a;
    }


    private static class TheResult {
        Matrix q;
        Matrix r;
        Matrix x;
        double error;
        Matrix xHouse;
        Matrix xGivens;
        double errorHouse;
        double errorGivens;

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
        public Matrix getX() {
            return x;
        }
        public void setX(Matrix x) {
            this.x = x;
        }
        public Matrix getXHouse() {
            return xHouse;
        }
        public void setXHouse(Matrix x) {
            this.xHouse = x;
        }
        public Matrix getXGivens() {
            return xGivens;
        }
        public void setXGivens(Matrix x) {
            this.xGivens = x;
        }
        public double getErrorHouse() {
            return errorHouse;
        }

        public void setErrorHouse(double error) {
            this.errorHouse = error;
        }
        public double getErrorGivens() {
            return errorGivens;
        }

        public void setErrorGivens(double error) {
            this.errorGivens = error;
        }
    }

}
