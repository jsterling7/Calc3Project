public class Test {
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

        Matrix a = new Matrix(3, 3);
        a.set(0, 0, 1);
        a.set(1, 0, (double) 1 / 2);
        a.set(2, 0, (double) 1 / 3);
        a.set(0, 1, (double) 1 / 2);
        a.set(1, 1, 1);
        a.set(2, 1, (double) 1 / 4);
        a.set(0, 2, (double) 1 / 3);
        a.set(1, 2, (double) 1 / 4);
        a.set(2, 2, 1);
        Matrix c = new Matrix(3, 1);
        c.set(0, 0, 0.1);
        c.set(1, 0, 0.1);
        c.set(2, 0, 0.1);
        Matrix approxTotal = Jacobi.calculateRandom(a, c);
        Matrix actualTotal = new Matrix(3,1);
        actualTotal.set(0, 0, (double) 9 / 190);
        actualTotal.set(1, 0, (double) 28 / 475);
        actualTotal.set(2, 0, (double) 33 / 475);
        System.out.println("");
        Jacobi.iterativeError(approxTotal, actualTotal);
    }
}
