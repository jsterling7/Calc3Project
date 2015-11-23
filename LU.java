import Jama.Matrix;

public class LU {
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
}
