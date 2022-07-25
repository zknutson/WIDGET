package widget;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class Output {
    private BufferedWriter o;

    public Output(String fileName) {
        try {
            o = new BufferedWriter(new FileWriter(fileName));
        } catch (Exception e) {
            e.getStackTrace();
        }
    }

    public void writeTitle(String[] str) {
        try {
            for (int i = 0; i < str.length - 1; i++) {
                o.write(str[i] + ",");
            }
            o.write(str[str.length - 1]);
            o.newLine();
        } catch (Exception e) {
            System.out.println(e);
            e.getStackTrace();
        }
    }

    public void writeEntry(double[] d) {
        try {
            for (int i = 0; i < d.length - 1; i++) {
                o.write(d[i] + ",");
            }
            o.write(String.valueOf(d[d.length - 1]));
            o.newLine();
        } catch (Exception e) {
            System.out.println(e);
            e.getStackTrace();
        }
    }

    public void close() {
        try {
            o.close();
        } catch (Exception e) {
            e.getStackTrace();
        }
    }
}