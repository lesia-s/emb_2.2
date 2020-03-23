import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;

public class Main {
    public static void main(String[] args) {
        JFrame frame = new JFrame("Lab4");
        frame.setLayout(new FlowLayout());

        Signal signal = new Signal(64, 10, 120);

        new Chart(frame, "x(t)", signal.getArr(), "Mx = " + signal.getMx() + "\n" + "Dx = " + signal.getDx());
        new Chart(frame, "", signal.getfArr(), "");

        frame.pack();
        frame.setVisible(true);
    }
}

class Signal {
    private final double[] xArr;
    private final double Mx;
    private final double Dx;

    public Signal(int N, int n, int freq) {

        this.xArr = new double[N];

        for (int i = 0; i < N; i++) {
            double x = 0;
            double A = Math.random();
            double fi = Math.random()*2*Math.PI;

            for (int b = 0, c = freq; b < n; b++, c += freq) {
                x += A * Math.sin(c * i + fi);
            }
            this.xArr[i] = x;
        }

        double mSum = 0;
        for (int i = 0; i < N; i++) {
            mSum += xArr[i];
        }

        this.Mx = mSum/N;

        double dSum = 0;
        for (int i = 0; i < N; i++) {
            dSum += Math.pow((xArr[i] - Mx), 2);
        }

        this.Dx = dSum/(N-1);
    }

    public double[] getfArr() {
        int N = xArr.length;
        double[] fArr = new double[N - 1];
        for (int p = 0; p < N - 1; p++) {
            double i1 = 0, i2 = 0, r1 = 0, r2 = 0;
            for (int k = 0; k < (N / 2) - 1; k++) {
                double temp1 = 4 * Math.PI / N * p * k;
                i1 += xArr[2 * k] * Math.sin(temp1);
                r1 += xArr[2 * k] * Math.cos(temp1);
                double temp2 = 2 * Math.PI / N * p * (2 * k + 1);
                i2 += xArr[2 * k + 1] * Math.sin(temp2);
                r2 += xArr[2 * k + 1] * Math.cos(temp2);
            }
            fArr[p] = Math.sqrt(Math.pow(r1 + r2, 2) + Math.pow(i1 + i2, 2));
        }
        return fArr;
    }

    public double[] getArr(){
        return xArr;
    }

    public double getDx(){
        return this.Dx;
    }

    public double getMx(){
        return this.Mx;
    }
}

class Chart {
    public Chart(JFrame frame, String name, double[] arr, String title) {

        XYSeries series = new XYSeries(name);

        for(int i = 0; i < arr.length; i++){
            series.add(i, arr[i]);
        }

        XYDataset xyDataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory
                .createXYLineChart(title, "", "",
                        xyDataset,
                        PlotOrientation.VERTICAL,
                        true, true, true);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(400, 270));

        frame.getContentPane()
                .add(chartPanel);
    }
}
