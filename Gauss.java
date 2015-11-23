import Jama.*;

public class Gauss {
    public static void main(String[] args) {
        double[][] test = {{0}, {0}, {0}};
        Matrix a = new Matrix(test);
        TheResult results = gs_iter(a, .00005, 100);
        System.out.println("The X Solution");
        results.getXSolution().print(2, 10);
        System.out.println("Number of Iterations");
        System.out.println(results.getIterations());
    }

    public static class TheResult {
        Matrix xSolution;
        int iterations;

        public Matrix getXSolution() {
            return xSolution;
        }

        public void setXSolution(Matrix x) {
            xSolution = x;
        }

        public int getIterations() {
            return iterations;
        }

        public void setIterations(int iterations) {
            this.iterations = iterations;
        }
    }

    public static TheResult gs_iter(Matrix x, double tolerance, int maxIterations) {
        TheResult theResult = new TheResult();
        double one = 1;
        double three = 3;

        double[][] lowTriMatrixArr = {{1, 0, 0}, {.5, 1, 0}, {one / three, .25, 1}};
        Matrix lowTriMatrix = new Matrix(lowTriMatrixArr);
        lowTriMatrix = inverseLowerTri(lowTriMatrix);

        double[][] bArr = {{.1}, {.1}, {.1}};
        Matrix b = new Matrix(bArr);
        b.print(2, 2);

        double[][] uArr = {{0, .5, one / three} ,  {0, 0, .25}, {0, 0, 0}};
        Matrix u = new Matrix(uArr);

        Matrix uTimesx = multiply(u, x);
        Matrix bMinuesuTimesx = b.minus(uTimesx);
        Matrix sol1 = multiply(lowTriMatrix, bMinuesuTimesx);
        Matrix sol2 = new Matrix(3 , 1);
        sol1.print(2, 2);
        sol2.print(2, 2);


        for (int i = 0; i < maxIterations; i++) {
            uTimesx = multiply(u, sol1);
            bMinuesuTimesx = b.minus(uTimesx);
            sol2 = multiply(lowTriMatrix, bMinuesuTimesx);
            Matrix check = sol2.minus(sol1);
            double checkInfinNorm = normInfinity(check);
            if (checkInfinNorm <= tolerance) {
                theResult.setXSolution(sol2);
                theResult.setIterations(i + 1);
                return theResult;
            } else {
                sol1 = sol2;
            }
        }

        return null;
    }


    public static Matrix inverseLowerTri(Matrix lowTriMatrix) {
        double[][] theMatrix = lowTriMatrix.getArrayCopy();
        double[][] theMatrixcopy = lowTriMatrix.getArrayCopy();
        for (int i = 0; i < theMatrix.length; i++) {
            for (int j = 0; j < theMatrix[i].length; j++) {
                if ((i == 0 && j == 0) || (i == 1 && j == 1) || (i == 2 && j == 2)) {
                    theMatrix[i][j] = (double) 1 / theMatrixcopy[i][j];
                }
                if (i == 1 && j == 0) {
                    theMatrix[i][j] = ((double) (-1) * theMatrixcopy[i][j]) / (theMatrixcopy[0][0] * theMatrixcopy[1][1]);
                }
                if (i == 2 && j == 0) {
                    theMatrix[i][j] = ((double) (-1) * theMatrixcopy[1][1] * theMatrixcopy[i][j] + theMatrixcopy[1][0] * theMatrixcopy[2][1]) /
                        (theMatrixcopy[0][0] * theMatrixcopy[1][1] * theMatrixcopy[2][2]);
                }
                if (i == 2 && j == 1) {
                    theMatrix[i][j] = ((-1) * theMatrixcopy[i][j]) / (theMatrixcopy[1][1] * theMatrixcopy[2][2]);
                }
            }
        }
        return lowTriMatrix.constructWithCopy(theMatrix);
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
}
