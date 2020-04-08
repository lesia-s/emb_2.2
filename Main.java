package additional;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

public class Main {
    public static void main(String[] args) {
        JFrame frame = new JFrame("Additional task");
        frame.setLayout(new FlowLayout());

        Signal signal = new Signal(5, 12, 10, 120);
        signal.countDependencies();

        new Chart(frame, "FFT", "DFT", signal.getFFT(), signal.getDFT(), "");

        frame.pack();
        frame.setVisible(true);
    }
}

class Signal {
    private final int n;
    private final int freq;
    private final int fromN;
    private final int toN;
    private HashMap<Integer, Long> FFT;
    private HashMap<Integer, Long> DFT;

    public Signal(int fromN, int toN, int n, int freq) {
        this.n = n;
        this.freq = freq;
        this.fromN = fromN;
        this.toN = toN;
    }

    public void countDependencies() {
        FFT = new HashMap<>();
        DFT = new HashMap<>();
        for (int i = this.fromN; i <= this.toN; i++) {
            int N = (int) Math.pow(2, i);
            double[] xArr = generateSignal(N, this.n, this.freq);
            FFT.put(N, countFFT(xArr));
            DFT.put(N, countDFT(xArr));
        }
    }

    public double[] generateSignal(int N, int n, int freq) {

        double[] xArr = new double[N];

        for (int i = 0; i < N; i++) {
            double x = 0;
            double A = Math.random();
            double fi = Math.random()*2*Math.PI;

            for (int b = 0, c = freq; b < n; b++, c += freq) {
                x += A * Math.sin(c * i + fi);
            }
            xArr[i] = x;
        }

        return xArr;
    }

    public long countFFT(double[] xArr) {
        long start = System.currentTimeMillis();
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
        return System.currentTimeMillis() - start;
    }

    public long countDFT(double[] xArr) {
        long start = System.currentTimeMillis();
        int N = xArr.length;
        double[][] fCos = new double[N-1][N-1];
        double[][] fSin = new double[N-1][N-1];
        for (int p = 0; p < N-1; p++) {
            for (int k = 0; k < N - 1; k++) {
                fCos[p][k] += Math.cos(2 * Math.PI / N * p * k);
                fSin[p][k] += Math.sin(2 * Math.PI / N * p * k);
            }
        }
        double[] re = new double[N - 1];
        double[] im = new double[N - 1];
        double[] fArr = new double[N - 1];
        for (int p = 0; p < N-1; p++) {
            double l = 0;
            double h = 0;
            for (int k = 0; k < N-1; k++) {
                l += xArr[k]*fCos[p][k];
                h += xArr[k]*fSin[p][k];
            }
            re[p] = l;
            im[p] = h;
        }
        for (int p = 0; p < N-1; p++) {
            fArr[p] = Math.sqrt(Math.pow(re[p], 2) + Math.pow(im[p], 2));
        }
        return System.currentTimeMillis() - start;
    }

    public HashMap<Integer, Long> getDFT() {
        return this.DFT;
    }

    public HashMap<Integer, Long> getFFT() {
        return this.FFT;
    }
}

class Chart {
    public Chart(JFrame frame, String nameFFT, String nameDFT, HashMap<Integer, Long> fft, HashMap<Integer, Long> dft, String title) {

        XYSeries seriesFFT = new XYSeries(nameFFT);
        XYSeries seriesDFT = new XYSeries(nameDFT);

        for (Map.Entry<Integer, Long> pair : fft.entrySet()) {
            seriesFFT.add(pair.getKey(), pair.getValue());
        }

        for (Map.Entry<Integer, Long> pair : dft.entrySet()) {
            seriesDFT.add(pair.getKey(), pair.getValue());
        }

        XYSeriesCollection xyDataset = new XYSeriesCollection();
        xyDataset.addSeries(seriesDFT);
        xyDataset.addSeries(seriesFFT);

        JFreeChart chart = ChartFactory
                .createXYLineChart(title, "N", "ms",
                        xyDataset,
                        PlotOrientation.VERTICAL,
                        true, true, true);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 540));

        frame.getContentPane()
                .add(chartPanel);
    }
}
