import jama.*;

public class SymmetricPascalMatrix {
    // LU
    public static Matrix[] lu_fact(Matrix m) {
		System.out.println("A = ");
		m.print(2, 1);
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

		// m.print(2,1);
		Matrix la = new Matrix(l);
		Matrix ua = new Matrix(u);
		Matrix[] ans = { la, ua };
		System.out.println("L = ");
		ans[0].print(2, 20);
		System.out.println("U = ");
		ans[1].print(2, 20);
		System.out.println("L*U = ");
		multiply(ans[0], ans[1]).print(2, 20);
		System.out.println("||LU - A||inf = "
				+ normInfinity(multiply(ans[0], ans[1]).minus(m)));
		return ans;
	}

	public static Matrix solve_lu_b(Matrix m, Matrix b) {
		Matrix[] ans = lu_fact(m);
		Matrix l = ans[0];
		Matrix u = ans[1];
		Matrix y = new Matrix(b.getArray());
		Matrix x = new Matrix(b.getRowDimension(), 1);
		for (int i = 0; i < y.getRowDimension(); i++) {
			double xx = b.get(i, 0);
			for (int j = 0; j < i; j++) {
				// System.out.println(l.get(i, j) + "   " + y.get(j, 0));
				xx -= l.get(i, j) * y.get(j, 0);
			}
			y.set(i, 0, xx);
		}
		for (int i = x.getRowDimension() - 1; i >= 0; i--) {
			// System.out.println("I iteration");
			double xx = y.get(i, 0);
			// System.out.println(xx);
			for (int j = m.getColumnDimension() - 1; j > i; j--) {
				xx -= u.get(i, j) * x.get(j, 0);
				// System.out.println("i and j: " + i + "," + j + " " + u.get(i,
				// j) + "  " + x.get(j, 0) + " " + xx);
			}
			x.set(i, 0, xx / u.get(i, i));
		}
		// System.out.println("Y: ");
		// y.print(2, 3);
		System.out.println("X: ");
		x.print(2, 3);
		return x;
	}

	public static Matrix multiply(Matrix a, Matrix b) {
		if (a.getColumnDimension() == b.getRowDimension()) {
			Matrix ans = new Matrix(a.getRowDimension(), b.getColumnDimension());
			for (int i = 0; i < ans.getRowDimension(); i++) {
				for (int j = 0; j < ans.getColumnDimension(); j++) {
					for (int k = 0; k < a.getColumnDimension(); k++) {
						ans.set(i, j,
								ans.get(i, j) + (a.get(i, k) * b.get(k, j)));
					}
				}
			}
			return ans;
		} else if (a.getColumnDimension() == 2 && a.getRowDimension() == 2) {
			if (b.getRowDimension() == 2 && b.getColumnDimension() == 1) {
				Matrix ans = new Matrix(2, 1);
				ans.set(0, 0, (a.get(0, 0) * b.get(0, 0)) + (a.get(0, 1) * b.get(1, 0)));
				ans.set(1, 0, (a.get(1, 0) * a.get(0, 0)) + (a.get(1, 1) * b.get(1, 0)));
				return ans;
			}
			return null;
		} else {
			return null;
		}
	}

	public static double normInfinity(Matrix m) {
		double[] ans = new double[m.getRowDimension()];
		for (int i = 0; i < ans.length; i++) {
			for (int j = 0; j < m.getColumnDimension(); j++) {
				ans[i] += Math.abs(m.get(i, j));
				// System.out.println(ans[i]);
			}
		}
		double answer = ans[0];
		for (int i = 1; i < ans.length; i++) {
			if (ans[i] > answer) {
				answer = ans[i];
			}
		}
		return answer;
	}

	public static void main(String args[]) {
		double[][] mat3 = { { 1, 1, 1, 1 }, { 1, 2, 3, 4 }, { 1, 3, 6, 10 },
				{ 1, 4, 10, 20 } };
		Matrix m = new Matrix(mat3);
		double[][] b2 = { { 1 }, { .5 }, { 1.0 / 3.0 }, { (1.0 / 4.0) } };
		// System.out.println(b2[2][0]);
		Matrix b = new Matrix(b2);
		solve_lu_b(m, b);
	}
    //--------------------------------------------
    //HouseHolder

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

    private static class TheResult {
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


    // ------------------------------------------------------------
    // Givens

    static Matrix q;
    static Matrix r;
    static double error;

    public static void qr_fact_givens(Matrix a) {
        double[][] copy = a.getArrayCopy();
        q = Matrix.identity(a.getRowDimension(),a.getColumnDimension());
        r = a;
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
        Matrix errorMatrix = multiply(q, r);
        errorMatrix.minusEquals(a);
        error = normInfinity(errorMatrix);
        System.out.print("Q:");
        q.print(7, 7);
        System.out.print("R:");
        r.print(7, 7);
        System.out.println("Error: " + Givens.getError());
    }
    public static Matrix solve_qr_b(Matrix a, Matrix b) {
        qr_fact_givens(a);
        Matrix y = multiply(q.transpose(), b);
        Matrix x = new Matrix(b.getRowDimension(), 1);
//        x.set(b.getRowDimension() - 1, 0, (y.get(b.getRowDimension() - 1, 0) / r.get(b.getRowDimension() - 1, b.getRowDimension() - 1)));
        for (int i = b.getRowDimension() - 1; i >= 0; i--) {
            double temp = y.get(i,0);
            for (int j = i + 1; j < b.getRowDimension(); j++) {
                 temp -= (r.get(i, j) * x.get(j, 0));
            }
            x.set(i, 0, temp/ r.get(i,i));
        }
        System.out.println("Y: ");
        y.print(2, 3);
        System.out.println("X: ");
        x.print(2, 3);
        return x;
    }

    public static Matrix getQ() {
        return q;
    }

    public static Matrix getR() {
        return r;
    }

    public static double getError() {
        return error;
    }

    public static void main(String[] args) {
        Matrix test = new Matrix(new double[4][4]);
        test.set(0, 0, 1);
        test.set(1, 0, 1);
        test.set(2, 0, 1);
        test.set(3, 0, 1);
        test.set(0, 1, 1);
        test.set(1, 1, 2);
        test.set(2, 1, 3);
        test.set(3, 1, 4);
        test.set(0, 2, 1);
        test.set(1, 2, 3);
        test.set(2, 2, 6);
        test.set(3, 2, 10);
        test.set(0, 3, 1);
        test.set(1, 3, 4);
        test.set(2, 3, 10);
        test.set(3, 3, 20);
        Givens.qr_fact_givens(test);
        Matrix b = new Matrix(new double [4][1]);
        b.set(0, 0, 1);
        b.set(1, 0, (double) 1 / 2);
        b.set(2, 0, (double) 1 / 3);
        b.set(3, 0, (double) 1 / 4);
        Givens.solve_qr_b(test, b);
    }




    // ------------------------------------------------------
    // Helper Methods

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
