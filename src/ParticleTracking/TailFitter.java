/*
 * Copyright (C) 2014 David Barry <david.barry at cancer.org.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package ParticleTracking;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.io.File;
import java.util.Random;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.MathArrays;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class TailFitter extends IsoGaussianFitter {

    private final double spatialRes, sigmaEst;
    double sqrt2 = Math.pow(2.0, 0.5);
    private static final int GAUSSIAN = 0, EMG = 1, EMG_PLUS_GAUSSIAN = 2;
    protected static final String eqLabels[] = {"Gaussian", "EMG", "Gaussian Plus EMG"};
    private final int eqChoice;

    public TailFitter(int eqChoice, double spatialRes, double sigmaEst) {
        super();
        this.eqChoice = eqChoice;
        this.spatialRes = spatialRes;
        this.sigmaEst = sigmaEst;
    }

    public void loadData(double[] xVals, double[] yVals, double[][] zVals) {
        this.xData = xVals;
        this.yData = yVals;
        this.zData = new double[xData.length * yData.length];
        for (int x = 0; x < xData.length; x++) {
            for (int y = 0; y < yData.length; y++) {
                this.zData[y * xData.length + x] = zVals[x][y];
            }
        }
        if (xData != null && yData != null) {
            numPoints = xVals.length * yVals.length;
            for (int i = xVals.length - 1; i >= 0; i--) {
                xData[i] -= xData[0];
            }
            for (int j = yVals.length - 1; j >= 0; j--) {
                yData[j] -= yData[0];
            }
        } else {
            numPoints = 0;
        }
    }

    boolean initialize(double sigmaEst) {
        if (sigmaEst <= 0.0 || xData == null || yData == null || zData == null) {
            return false;
        }
        switch (eqChoice) {
            case (TailFitter.GAUSSIAN):
                numParams = 4;
                break;
            case (TailFitter.EMG):
                numParams = 5;
                break;
            case (TailFitter.EMG_PLUS_GAUSSIAN):
                numParams = 5;
                break;
        }
        // Calculate some things that might be useful for predicting parametres
        numVertices = numParams + 1;      // need 1 more vertice than parametres,
        simp = new double[numVertices][numVertices];
        next = new double[numVertices];
        maxIter = IterFactor * numParams * numParams;  // Where does this estimate come from?
        restarts = defaultRestarts;
        nRestarts = 0;
        Random r = new Random();
        double noise = 0.1;
        simp[0][0] = 1.0 + r.nextDouble() * noise; //A, lambda
        simp[0][1] = 0.95 + r.nextDouble() * noise; //mu
        simp[0][2] = 0.1 + r.nextDouble() * noise; //sigma
        simp[0][3] = 0.0 + r.nextDouble() * noise; //noise
        if (eqChoice != TailFitter.GAUSSIAN) {
            simp[0][4] = 0.5 + r.nextDouble() * noise; //nu
        }
        return true;
    }

    public double evaluate1DEMG(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }
        double p22 = p[2] * p[2];
        double a = 0.5 * p[0] * (2.0 * p[1] + p[0] * p22 - 2.0 * x);
        double b = (p[1] + p[0] * p22 - x) / (sqrt2 * p[2]);

        return p[4] * p[0] * Math.exp(a) * Erf.erfc(b) + p[3];
    }

    public double evaluate1DEMGDerivative(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }
        double p22 = p[2] * p[2];
        double a = (p[0] / 2.0) * (2.0 * p[1] + p[0] * p22 - 2.0 * x);
        double b = (p[1] + p[0] * p22 - x) / (sqrt2 * p[2]);
        double b2 = b * b;
        double c = p[0] * p[3] * sqrt2 / (p[2] * Math.sqrt(Math.PI));
        double d = p[0] * p[0] * p[3];

        return Math.exp(a) * (c * Math.exp(-b2) - d * Erf.erfc(b));
    }

    double findPeak(int Nmax, double a, double b, double tol, double[] p) {
        int N = 1;
        while (N < Nmax) {
            double c = (a + b) / 2;
            double fa, fc;
            if (eqChoice == TailFitter.EMG) {
                fa = evaluate1DEMGDerivative(p, a);
                fc = evaluate1DEMGDerivative(p, c);
            } else if (eqChoice == TailFitter.EMG_PLUS_GAUSSIAN) {
                fa = evaluate1DGaussianPlusEMG1stDerivative(p, a);
                fc = evaluate1DGaussianPlusEMG1stDerivative(p, c);
            } else {
                return Double.NaN;
            }
            if (fc == 0.0 || (b - a) / 2 < tol) {
                return c;
            }
            N++;
            if (fa * fc > 0.0) {
                a = c;
            } else {
                b = c;
            }
        }
        return Double.NaN;
    }

    double estimateLength(int Nmax, double a, double b, double tol, double intens, double[] p) {
        int N = 1;
        double I = intens * (1.0 - p[3]) + p[3];
        while (N < Nmax) {
            double c = (a + b) / 2;
            double fa, fc;
            if (eqChoice == TailFitter.EMG) {
                fa = evaluate1DEMG(p, a);
                fc = evaluate1DEMG(p, c);
            } else if (eqChoice == TailFitter.EMG_PLUS_GAUSSIAN) {
                fa = evaluate1DGaussianPlusEMG(p, a);
                fc = evaluate1DGaussianPlusEMG(p, c);
            } else {
                return Double.NaN;
            }
            if (fc == I || (b - a) / 2 < tol) {
                return c;
            }
            N++;
            if ((fa - I) * (fc - I) > 0.0) {
                a = c;
            } else {
                b = c;
            }
        }
        return Double.NaN;
    }

    public double evaluate1DGaussianPlusEMG(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }
        double lambda = 0.7960;
        double mu = 0.3675 + .75;
        double sigma = 0.2185;
        double a = 0.5 * lambda * (2.0 * mu + lambda * sigma * sigma - 2.0 * x);
        double b = (mu + lambda * sigma * sigma - x) / (Math.sqrt(2.0) * sigma);

        return evaluate1DGaussian(p, x) + p[4] * lambda * Math.exp(a) * Erf.erfc(b);
    }

    public double evaluate1DGaussianPlusEMG1stDerivative(double[] p, double x) {
        return evaluate1DEMGDerivative(p, x) + ((evaluate1DGaussianPlusEMG(p, x) - p[3]) * (p[1] - x) / (p[2] * p[2]));
    }

    public double evaluate1DGaussian(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }

        return p[0] * Math.exp(-0.5 * (Math.pow((x - p[1]) / p[2], 2.0))) + p[3];
    }

    public double evaluate2D(double[] p, double xVal, double y) {
        if (p == null) {
            return Double.NaN;
        }
        double v = Math.pow((y - p[5]) / (p[6] * sqrt2), 2.0);

        return xVal * Math.exp(-0.5 * v) + p[4];
    }

    protected boolean sumResiduals(double[] x) {
        if (x == null) {
            return false;
        }
        /*
         * x[numParams] = sumResiduals(x, xData, yData, zData); return true;
         */
        double e;
        x[numParams] = 0.0;
        double tail1d[] = buildTail(x);
        for (int i = 0; i < xData.length; i++) {
            for (int j = 0; j < yData.length; j++) {
                e = tail1d[i + xData.length / 2] - zData[j * xData.length + i];
//                e = tail1d[i] - zData[j * xData.length + i];
                x[numParams] = x[numParams] + (e * e);
            }
        }
        return true;
    }

    double[] buildTail(double[] p) {
        double[] emg = new double[xData.length];
        Gaussian gauss = new Gaussian((emg.length - 1.0) / 2.0, sigmaEst / spatialRes);
        double[] gaussian = new double[emg.length];
        for (int i = 0; i < emg.length; i++) {
            gaussian[i] = gauss.value(i);
            switch (eqChoice) {
                case (TailFitter.EMG_PLUS_GAUSSIAN):
                    emg[i] = evaluate1DGaussianPlusEMG(p, xData[i]);
                    break;
                case (TailFitter.GAUSSIAN):
                    emg[i] = evaluate1DGaussian(p, xData[i]);
                    break;
                case (TailFitter.EMG):
                    emg[i] = evaluate1DEMG(p, xData[i]);
                    break;
            }
        }
        return MathArrays.convolve(emg, gaussian);
//        return emg;
    }

    public void printParams() {
        double params[] = getParams();
        for (int i = 0; i < numParams; i++) {
            System.out.print("p[" + String.valueOf(i) + "]:," + params[i] + ",");
        }
        if (eqChoice != TailFitter.GAUSSIAN) {
            double peak = findPeak(10000, params[1] - 2.0 * params[2],
                    params[1] + 2.0 * params[2], 1.0E-10, params);
            double length = estimateLength(10000, peak, peak + 46.0 * spatialRes,
                    1.0E-10, 0.1, params);
            System.out.print("Peak:,x=," + peak + ",Length:,x=," + length);
        }
        System.out.println();
    }

    void printImage(File directory, String filename) {
        FloatProcessor deconvolved = new FloatProcessor(xData.length, yData.length);
        FloatProcessor convolved = new FloatProcessor(xData.length, yData.length);
        FloatProcessor derivative = new FloatProcessor(xData.length, yData.length);
        double tail1d[] = buildTail(simp[best]);
        for (int y = 0; y < yData.length; y++) {
            for (int x = 0; x < xData.length; x++) {
                convolved.putPixelValue(x, y, tail1d[x + xData.length / 2]);
                switch (eqChoice) {
                    case (TailFitter.EMG_PLUS_GAUSSIAN):
                        deconvolved.putPixelValue(x, y, evaluate1DGaussianPlusEMG(simp[best], xData[x]));
                        derivative.putPixelValue(x, y, evaluate1DGaussianPlusEMG1stDerivative(simp[best], xData[x]));
                        break;
                    case (TailFitter.EMG):
                        deconvolved.putPixelValue(x, y, evaluate1DEMG(simp[best], xData[x]));
                        derivative.putPixelValue(x, y, evaluate1DEMGDerivative(simp[best], xData[x]));
                        break;
                    case (TailFitter.GAUSSIAN):
                        deconvolved.putPixelValue(x, y, evaluate1DGaussian(simp[best], xData[x]));
                        break;
                }
            }
        }
        IJ.saveAs((new ImagePlus("", convolved)), "text image", filename + "Convolved.txt");
        IJ.saveAs((new ImagePlus("", deconvolved)), "text image", filename + "Deconvolved.txt");
        if (eqChoice != TailFitter.GAUSSIAN) {
            IJ.saveAs((new ImagePlus("", derivative)), "text image", filename + "Derivative.txt");
        }
    }
}
