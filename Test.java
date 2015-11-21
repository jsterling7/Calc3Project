public class Test {
    public static void main(String[] args) {
        Matrix test = new Matrix(new double[4][4]);
        test.set(0,0,1);
        test.set(1,0,1);
        test.set(2,0,1);
        test.set(3,0,1);
        test.set(0,1,1);
        test.set(1,1,2);
        test.set(2,1,3);
        test.set(3,1,4);
        test.set(0,2,1);
        test.set(1,2,3);
        test.set(2,2,6);
        test.set(3,2,10);
        test.set(0,3,1);
        test.set(1,3,4);
        test.set(2,3,10);
        test.set(3, 3, 20);
        System.out.print("Original:");
        test.print(4, 4);
        Givens.solve_qr_givens(test);
        System.out.print("Q:");
        Givens.getQ().print(4, 4);
        System.out.print("R:");
        Givens.getR().print(4, 4);
        System.out.println("Error: " + Givens.getError());
    }
}
