public class Givens {
    static Matrix q;
    static Matrix r;

    public static void solve_qr_givens(Matrix a) {
        double[][] copy = a.getArrayCopy();
        r = a;
        int currentColumn = 0;
        int currentRow = 0;
        for (currentColumn = 0; currentColumn < a.getColumnDimension(); currentColumn++) {
            for (currentRow = currentColumn; currentRow < a.getRowDimension() - 1; currentRow++) {
                if (copy[currentColumn][currentColumn] != 0 && copy[currentColumn][currentRow + 1] != 0) {
                    copy[currentColumn][currentColumn] = copy[currentColumn][currentColumn]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2))
                            + Math.pow(copy[currentColumn][currentRow + 1], 2));
                    copy[currentColumn][currentRow + 1] = copy[currentColumn][currentRow + 1]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2))
                            + Math.pow(copy[currentColumn][currentRow + 1], 2));
                    copy[currentRow + 1][currentColumn] = -(copy[currentColumn][currentColumn]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2))
                            + Math.pow(copy[currentColumn][currentRow + 1], 2)));
                    copy[currentRow + 1][currentRow + 1] = copy[currentColumn][currentRow + 1]
                            / (Math.sqrt(Math.pow(copy[currentColumn][currentColumn], 2))
                            + Math.pow(copy[currentColumn][currentRow + 1], 2));
                    for (int x = 0; x < a.getColumnDimension(); x++) {
                        for (int y = 0; y < a.getRowDimension(); y++) {
                            if (x != currentColumn && y != currentRow + 1) {
                                if (x==y) {
                                    copy[x][y] = 1;
                                } else {
                                    copy[x][y] = 0;
                                }
                            }
                        }
                    }
                    Matrix rtemp = Matrix.constructWithCopy(copy);
                    r = rtemp.arrayTimes(r);
                    copy = r.getArrayCopy();
                }
            }
        }
    }

    public static Matrix getQ() {
        return q;
    }

    public static Matrix getR() {
        return r;
    }
}
