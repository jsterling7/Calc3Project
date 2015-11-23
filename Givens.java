public class Givens {
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
}
