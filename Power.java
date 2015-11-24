import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import Jama.Matrix;

//http://college.cengage.com/mathematics/larson/elementary_linear/5e/students/ch08-10/chap_10_3.pdf
public class Power {
	static double[] dets = new double[1000];
	static double[] traces = new double[1000];
	static double[] iterations = new double[1000];
	static double[] detsI = new double[1000];
	static double[] tracesI = new double[1000];
	static double[] iterationsI = new double[1000];
	static Matrix[] eigenVectors = new Matrix[1000];
	static double[] eigenValues = new double[1000];
	static Matrix[] eigenVectorsI = new Matrix[1000];
	static double[] eigenValuesI = new double[1000];
	static int index = 0;
	static boolean inverse = false;
	static FileWriter writer = null;
	static int numfail = 0;
	static boolean input = false;

	public static Matrix[] make2x2(int n) {
		Matrix[] ans = new Matrix[n];
		for (int i = 0; i < n; i++) {
			double[][] mat = new double[2][2];
			mat[0][1] = (Math.random() * 4) - 2;
			mat[0][0] = (Math.random() * 4) - 2;
			mat[1][0] = (Math.random() * 4) - 2;
			mat[1][1] = (Math.random() * 4) - 2;
			ans[i] = new Matrix(mat);
		}
		return ans;
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
				ans.set(0,0,(a.get(0, 0) * b.get(0, 0))+ (a.get(0, 1) * b.get(1, 0)));
				ans.set(1,0,(a.get(1, 0) * a.get(0, 0))+ (a.get(1, 1) * b.get(1, 0)));
				return ans;
			}
		}
		return null;
	}

	public static Matrix inverse2x2(Matrix m) {
		Matrix ans = new Matrix(2, 2);
		double det = determinateOf2x2(m);
		ans.set(0, 0, m.get(1, 1) * det);
		ans.set(1, 1, m.get(0, 0) * det);
		ans.set(0, 1, m.get(0, 1) * -1 * det);
		ans.set(1, 0, m.get(1, 0) * -1 * det);
		return ans;
	}

	public static int power_method(Matrix a, Matrix v, double tol, int n) {
		Matrix u = new Matrix(v.getArrayCopy());
		Matrix lastAns = new Matrix(v.getRowDimension(), v.getColumnDimension());
		for (int i = 0; i < n; i++) {
			u = simplifyForPower(multiply(a, u));
			if (equal(u, lastAns, tol)) {
				if (inverse) {
					eigenValuesI[index] = rayleigh(a, u);
					eigenVectorsI[index] = u;
				} else {
					eigenValues[index] = rayleigh(a, u);
					eigenVectors[index] = u;
				}
				if(input) {
					System.out.println("For given matrix");
					a.print(2,5);
					System.out.println();
					System.out.println("It took " + i + " iterations");
					System.out.println("With a margain for error of .00005");
					System.out.println("The largest eigenvector is ");
					u.print(2, 5);
					System.out.println("The corresponding eigenvalue is ");
					System.out.println(rayleigh(a, u) + "\n \n");
				}
				return i;
			}
			lastAns = new Matrix(u.getArrayCopy());
		}
		numfail++;
		if (input) {
			System.out.println("For given matrix");
			a.print(2,5);
			System.out.println(n + " iterations was not enough to draw a meaningful result");
		}
		return n;
	}

	private static double rayleigh(Matrix a, Matrix x) {
		return Math.round(dotProduct(multiply(a, x), x) / dotProduct(x, x));
	}

	private static double dotProduct(Matrix a, Matrix b) {
		if (a.getColumnDimension() > 1 || b.getColumnDimension() > 1
				|| (a.getRowDimension() != b.getRowDimension())) {
			System.out.println("error in dot product");
			return -999999999999.0;
		} else {
			double ans = 0;
			for (int i = 0; i < a.getRowDimension(); i++) {
				ans += a.get(i, 0) * b.get(i, 0);
			}
			return ans;
		}
	}


	private static boolean equal(Matrix a, Matrix b, double tol) {
		if (a.getColumnDimension() == b.getColumnDimension() && a.getRowDimension() == b.getRowDimension()) {
			for (int i = 0; i < a.getRowDimension(); i++) {
				for (int j = 0; j < b.getColumnDimension(); j++) {
					if (Math.abs(a.get(i, j) - b.get(i, j)) >= tol) {
						return false;
					}
				}
			}
			return true;
		}
		return false;
	}

	private static Matrix simplifyForPower(Matrix m) {
		double biggest = Double.MIN_VALUE;
		for (int i = 0; i < m.getRowDimension(); i++) {
			for (int j = 0; j < m.getColumnDimension(); j++) {
				if (Math.abs(m.get(i, j)) >= biggest) {
					biggest = m.get(i, j);
				}
			}
		}
		for (int i = 0; i < m.getRowDimension(); i++) {
			for (int j = 0; j < m.getColumnDimension(); j++) {
				m.set(i, j, m.get(i, j) / biggest);
			}
		}
		return m;
	}

	private static double trace(Matrix m) {
		double ans = 0;
		for (int i = 0; i < m.getRowDimension(); i++) {
			ans += m.get(i, i);
		}
		return ans;
	}

	private static double determinateOf2x2(Matrix m) {
		return (1 / (m.get(0, 0) * m.get(1, 1) - m.get(1, 0) * m.get(0, 1)));
	}

	public static void main(String[] args) {
		File file = new File("answer.txt");
		try {
			writer = new FileWriter(file);
			double[][] mat2 = { { 1 }, { 1 } };
			Matrix a = new Matrix(mat2);

			Matrix[] mat = make2x2(1000);
			Matrix[] matI = new Matrix[1000];
			for (int i = 0; i < 1000; i++) {
				matI[i] = inverse2x2(mat[i]);
				inverse = false;
				index = i;
				writer.write((i + 1) + "th power method" + "\n");
				iterations[i] = power_method(mat[i], a, 0.00005, 100);
				dets[i] = determinateOf2x2(mat[i]);
				traces[i] = trace(mat[i]);

				writer.write("For matrix A = " + "\n");
				// a.print(2,5);
				double[][] arr = mat[i].getArray();
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						writer.write(String.valueOf(arr[j][k]) + " ");
					}
					writer.write("\n");
				}
				if (iterations[i] < 100 ) {
					writer.write("It took " + iterations[i] + " iterations"
							+ " \n");
					writer.write("With a margain for error of 0.00005" + " \n");
					writer.write("The largest eigenvector is " + " \n");
					arr = eigenVectors[i].getArray();
					for (int j = 0; j < arr.length; j++) {
						for (int k = 0; k < arr[0].length; k++) {
							writer.write(String.valueOf(arr[j][k]) + " ");
						}
						writer.write("\n");
					}
					writer.write("The corresponding eigenvalue is " + "\n");
					writer.write(String.valueOf(eigenValues[i]) + "\n");
					writer.write("The trace of A is " + trace(mat[i]) + "\n");
					writer.write("The determinate of A is "
							+ determinateOf2x2(mat[i]) + "\n \n");
				} else {
					writer.write("After 100 iterations, failure \n \n");
				}
				inverse = true;
				writer.write((i + 1) + "th power method inverse \n");
				iterationsI[i] = power_method(matI[i], a, 0.00005, 100);
				detsI[i] = determinateOf2x2(matI[i]);
				tracesI[i] = trace(matI[i]);
				writer.write("For matrix A = " + "\n");
				// a.print(2,5);
				arr = matI[i].getArray();
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						writer.write(String.valueOf(arr[j][k]) + " ");
					}
					writer.write("\n");
				}
				if (iterationsI[i] < 100) {
				writer.write("It took " + iterationsI[i] + " iterations"
						+ " \n");
				writer.write("With a margain for error of 0.00005" + " \n");
				writer.write("The largest eigenvector is " + " \n");
				arr = eigenVectorsI[i].getArray();
				for (int j = 0; j < arr.length; j++) {
					for (int k = 0; k < arr[0].length; k++) {
						writer.write(String.valueOf(arr[j][k]) + " ");
					}
					writer.write("\n");
				}
				writer.write("The corresponding eigenvalue is " + "\n");
				writer.write(String.valueOf(eigenValuesI[i]) + "\n");
				writer.write("The trace of A is " + trace(matI[i]) + "\n");
				writer.write("The determinate of A is "
						+ determinateOf2x2(matI[i]) + "\n \n");
				} else {
					writer.write("After 100 iterations, failure \n \n");
				}
			}
		} catch (IOException e) {
			e.printStackTrace(); // I'd rather declare method with throws
									// IOException and omit this catch.
		} finally {
			if (writer != null)
				try {
					writer.close();
				} catch (IOException ignore) {
				}
		}
		System.out.printf("File holding the information on the generated matricies is located at %s%n", file.getAbsolutePath());
		System.out.println(numfail + "failed");
		input = true;

//		double[][] mat2 = { { 1 }, { 1 } };
//		Matrix a = new Matrix(mat2);
//
//		Matrix[] mat = make2x2(1000);
//		Matrix[] matI = new Matrix[1000];
//		for (int i = 0; i < 1000; i++) {
//			matI[i] = inverse2x2(mat[i]);
//			System.out.println((i + 1) + "th power method");
//			iterations[i] = power_method(mat[i], a, .00005, 100);
//			dets[i] = determinateOf2x2(mat[i]);
//			traces[i] = trace(mat[i]);
//			System.out.println((i + 1) + "th power method inverse");
//			iterationsI[i] = power_method(matI[i], a, .00005, 100);
//			detsI[i] = determinateOf2x2(matI[i]);
//			tracesI[i] = trace(matI[i]);
//		}

		// System.out.println("Iterations");
		// for(int i = 0; i < 1000; i++) {
		// System.out.println(iterations[i]+ " ");
		// }
		// System.out.println();
		// System.out.println("dets");
		// for(int i = 0; i < 1000; i++) {
		// System.out.println(dets[i]+ " ");
		// }
		// System.out.println();
		// System.out.println("Traces");
		// for(int i = 0; i < 1000; i++) {
		// System.out.println(traces[i]+ " ");
		// }
		// System.out.println();
		// System.out.println("Iterations I");
		// for(int i = 0; i < 1000; i++) {
		// System.out.println(iterationsI[i]+ " ");
		// }
		// System.out.println();
		// System.out.println("dets I");
		// for(int i = 0; i < 1000; i++) {
		// System.out.println(detsI[i]+ " ");
		// }
		// System.out.println();
		// System.out.println("Traces I");
		// for(int i = 0; i < 1000; i++) {
		// System.out.println(tracesI[i]+ " ");
		// }
	}
}
