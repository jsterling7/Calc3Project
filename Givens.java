public class Givens {
    static Matrix q;
    static Matrix r;
    static double error;

    public static void solve_qr_givens(Matrix a) {
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
}
