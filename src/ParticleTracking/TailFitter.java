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

import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.ZProjector;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.MathArrays;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class TailFitter extends IsoGaussianFitter {

    private static double spatialRes = 0.133333;
    private static double sigmaEst = 0.120;
    double sqrt2 = Math.pow(2.0, 0.5);
    int minSlices = 20;
    private static final int GAUSSIAN = 0, EMG = 1, EMG_PLUS_GAUSSIAN = 2, C1 = 0, C2 = 1;
    private static String protein = "ALLSTARS";
    private boolean randomize = false, cellByCell = false;
    private final String channelLabels[] = {"C0_Output", "C1_Output"};
    private final String eqLabels[] = {"Gaussian", "EMG", "Gaussian Plus EMG"};
    private final String checkBoxLabels[] = {"Randomize", "Cell by Cell"};
    private boolean checks[] = {randomize, cellByCell};
    private static int chanChoice = 1, eqChoice = 1, iterations = 1;
    double minVersion = 5.017;

    public static void main(String args[]) {
        TailFitter tf = new TailFitter();
        tf.showDialog();
        tf.run();
        System.exit(0);
    }

    public void run() {
        File parentDir = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails"), "Select input folder", true);
        if (parentDir == null) {
            return;
        }
        String channel = channelLabels[chanChoice];
        String chan = channel.substring(0, 2);
        ArrayList<File> subDirs = new ArrayList();
        listSubDirs(parentDir, subDirs, protein);
        ArrayList<File> selectedSubDirs = new ArrayList();
        selectSubDirs(subDirs, selectedSubDirs, channelLabels[chanChoice]);
        ImagePlus temp = IJ.openImage(selectedSubDirs.get(0).listFiles()[0].getAbsolutePath());
        ImageProcessor tempIP = temp.getProcessor();
        int stackWidth = tempIP.getWidth();
        int stackHeight = tempIP.getHeight();
        temp.close();
        for (int i = 0; i < iterations; i++) {
            System.out.print(i + ",");
            ImageProcessor stackAverage;
            if (cellByCell) {
                stackAverage = projectStack(buildStackAverageTailByTail(stackWidth, stackHeight, selectedSubDirs, chan, randomize));
            } else {
                stackAverage = projectStack(buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, randomize));
            }
            Rectangle cropRoi = new Rectangle(0, 15, stackAverage.getWidth() - 16, 1);
            stackAverage.setRoi(cropRoi);
            stackAverage = stackAverage.crop();
            ImageStatistics stats = stackAverage.getStatistics();
            double max = stats.max;
            double min = stats.min;
            stackAverage.subtract(min);
            stackAverage.multiply(1.0 / (max - min));
            int width = stackAverage.getWidth();
            int height = stackAverage.getHeight();
            double xVals[] = new double[width];
            double yVals[] = new double[height];
            double zVals[][] = new double[width][height];
            for (int y = 0; y < height; y++) {
                yVals[y] = y * spatialRes;
                for (int x = 0; x < width; x++) {
                    xVals[x] = x * spatialRes;
                    zVals[x][y] = stackAverage.getPixelValue(x, y);
                }
            }
            loadData(xVals, yVals, zVals);
            doFit(TailFitter.sigmaEst);
            printParams();
        }
        generateImages(selectedSubDirs, chan, parentDir, stackWidth, stackHeight);
    }

    void listSubDirs(File directory, ArrayList<File> subDirs, String key) {
        File files[] = directory.listFiles();
        for (File file : files) {
            if (file.isDirectory()) {
                if (file.getName().contains(key)) {
                    subDirs.add(file);
                }
                listSubDirs(file, subDirs, key);
            }
        }
    }

    ArrayList<File> selectSubDirs(ArrayList<File> subDirs, ArrayList<File> selectedSubDirs, String key) {
        GenericDialog gd = new GenericDialog("Select subdirectories");
        for (File subDir : subDirs) {
            gd.addCheckbox(subDir.getAbsolutePath(), false);
        }
        gd.showDialog();
        if (!gd.wasOKed()) {
            return null;
        }
        for (int i = 0; i < subDirs.size(); i++) {
            if (gd.getNextBoolean()) {
                File[] subSubDirs = subDirs.get(i).listFiles();
                for (File subSubDir : subSubDirs) {
                    if (subSubDir.isDirectory() && subSubDir.getName().contains("Capture")) {
                        File[] resultsDirs = subSubDir.listFiles();
                        for (File resultsDir : resultsDirs) {
                            if (resultsDir.isDirectory()
                                    && (resultsDir.getName().contains("Particle Tracker"))) {
                                Scanner scanner = (new Scanner(resultsDir.getName())).useDelimiter("[_v]");
                                double version = 0;
                                while (scanner.hasNext()) {
                                    if (scanner.hasNextDouble()) {
                                        version = scanner.nextDouble();
                                    } else {
                                        scanner.next();
                                    }
                                }
                                if (version >= minVersion) {
                                    selectedSubDirs.add(new File(resultsDir.getAbsolutePath() + "/" + key));
                                    System.out.println(selectedSubDirs.get(selectedSubDirs.size() - 1).getAbsolutePath());
                                }
                            }
                        }
                    }
                }
            }
        }
        return selectedSubDirs;
    }

    ImageStack buildStackAverageTailByTail(int stackwidth, int stackheight, ArrayList<File> subDirs, String channel, boolean randomize) {
        ImageStack overallStack = new ImageStack(stackwidth, stackheight);
        int nDirs = subDirs.size();
        Random r = new Random();
        int count = 0;
        for (int j = 0; j < nDirs; j++) {
            File files[] = subDirs.get(j).listFiles();
            ArrayList<ArrayList> sortedFiles = sortFiles(files, channel);
            int n = sortedFiles.size();
            ImageStack cellStack = new ImageStack(stackwidth, stackheight);
            for (int i = 0; i < n; i++) {
                ArrayList<File> theseFiles = sortedFiles.get(i);
                int m = theseFiles.size();
                count += m;
                ImageStack tailStack = new ImageStack(stackwidth, stackheight);
                for (int l = 0; l < m; l++) {
                    int fileindex = randomize ? r.nextInt(m) : l;
                    ImagePlus imp = IJ.openImage(theseFiles.get(fileindex).getAbsolutePath());
                    tailStack.addSlice(imp.getProcessor());
                    imp.close();
                }
                if (tailStack.getSize() > 0) {
                    cellStack.addSlice(projectStack(tailStack));
                }
            }
            if (cellStack.getSize() > 0) {
                overallStack.addSlice(projectStack(cellStack));
            }
        }
        System.out.println("Images: " + String.valueOf(count));
        return overallStack;
//        if (overallStack.getSize() > 0) {
//            return projectStack(overallStack);
//        } else {
//            return null;
//        }
    }

    ImageStack buildStackAverageOverall(int stackwidth, int stackheight, ArrayList<File> subDirs, String channel, boolean randomize) {
        ImageStack output = new ImageStack(stackwidth, stackheight);
        int nDirs = subDirs.size();
        Random r = new Random();
        ArrayList<File> allFiles = new ArrayList();
        int nTails = 0;
        for (int j = 0; j < nDirs; j++) {
            File files[] = subDirs.get(j).listFiles();
            ArrayList<ArrayList> sortedFiles = sortFiles(files, channel);
            int n = sortedFiles.size();
            nTails += n;
            for (int i = 0; i < n; i++) {
                ArrayList<File> theseFiles = sortedFiles.get(i);
                int m = theseFiles.size();
                for (int l = 0; l < m; l++) {
                    allFiles.add(theseFiles.get(l));
                }
            }
        }
        int n = allFiles.size();
        for (int k = 0; k < n; k++) {
            int fileIndex = randomize ? r.nextInt(n) : k;
            ImagePlus imp = IJ.openImage(allFiles.get(fileIndex).getAbsolutePath());
            ImageProcessor ip = imp.getProcessor();
            ImageStatistics stats = ip.getStatistics();
            double max = stats.max;
            double min = stats.min;
            ip.subtract(min);
            ip.multiply(1.0 / (max - min));
            output.addSlice(ip);
            imp.close();
        }
        System.out.println("Cells: " + nDirs + " Tails: " + nTails + " Images: " + String.valueOf(output.getSize()));
        return output;
    }

    ImageProcessor projectStack(ImageStack stack) {
        ZProjector zproj = new ZProjector(new ImagePlus("", stack));
        zproj.setMethod(ZProjector.AVG_METHOD);
        zproj.doProjection();
        return zproj.getProjection().getProcessor();
    }

    ArrayList<ArrayList> sortFiles(File[] unsorted, String channel) {
        ArrayList<ArrayList> files = new ArrayList();
        int n = unsorted.length;
        for (int i = 0; i < n; i++) {
            String name = unsorted[i].getName();
            Scanner scanner = new Scanner(name).useDelimiter("-");
            String thisChannel = scanner.next();
            if (channel.equalsIgnoreCase(thisChannel)) {
                int index = Integer.parseInt(scanner.next());
                while (files.size() < index) {
                    files.add(new ArrayList<File>());
                }
                files.get(index - 1).add(unsorted[i]);
            }
        }
        return files;
    }

    public void generateImages(ArrayList<File> selectedSubDirs, String chan, File directory, int stackWidth, int stackHeight) {
        ImageStack stack;
        if (cellByCell) {
            stack = buildStackAverageTailByTail(stackWidth, stackHeight, selectedSubDirs, chan, false);
        } else {
            stack = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, false);
        }
//        ImageStack stack = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, false);
        ZProjector zproj = new ZProjector(new ImagePlus("", stack));
        zproj.setMethod(ZProjector.AVG_METHOD);
        zproj.doProjection();
        ImageProcessor stackAverage = zproj.getProjection().getProcessor();
        String fileBaseName = directory + "/" + protein + "_" + eqLabels[eqChoice] + "_" + "Randomize-" + randomize + "_" + "CellByCell-" + cellByCell + "_" + channelLabels[chanChoice] + "_";
        IJ.saveAs((new ImagePlus("", stackAverage)), "text image", fileBaseName + "Average.txt");
        zproj.setMethod(ZProjector.SD_METHOD);
        zproj.doProjection();
        IJ.saveAs((new ImagePlus("", zproj.getProjection().getProcessor())), "text image", fileBaseName + "SD.txt");
        Rectangle cropRoi = new Rectangle(0, 15, stackAverage.getWidth() - 2, 1);
        stackAverage.setRoi(cropRoi);
        stackAverage = stackAverage.crop();
        ImageStatistics stats = stackAverage.getStatistics();
        double max = stats.max;
        double min = stats.min;
        stackAverage.subtract(min);
        stackAverage.multiply(1.0 / (max - min));
        printImage(directory, fileBaseName);
        IJ.saveAs((new ImagePlus("", stackAverage)), "text image", fileBaseName + "NormAverage.txt");
    }

    public TailFitter() {
        super();
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
            System.out.print("Peak:,x=," + findPeak(10000, params[1] - 2.0 * params[2],
                    params[1] + 2.0 * params[2], 1.0E-10, params) + ",");
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

    public boolean showDialog() {
        GenericDialog gd = new GenericDialog("Tail Fitter");
        gd.addStringField("Protein Label:", protein, 5);
        gd.addRadioButtonGroup("Select Channel:", channelLabels, 1, 2, channelLabels[chanChoice]);
        gd.addChoice("Select Equation to Fit:", eqLabels, eqLabels[eqChoice]);
        gd.addCheckboxGroup(1, 2, checkBoxLabels, checks);
        gd.addNumericField("Iterations: ", iterations, 0);
        gd.addNumericField("Min Version: ", minVersion, 3);
        gd.showDialog();
        if (!gd.wasOKed()) {
            return false;
        }
        protein = gd.getNextString();
        chanChoice = gd.getNextRadioButton().equalsIgnoreCase(channelLabels[0]) ? TailFitter.C1 : TailFitter.C2;
        eqChoice = gd.getNextChoiceIndex();
        randomize = gd.getNextBoolean();
        cellByCell = gd.getNextBoolean();
        iterations = (int) Math.round(gd.getNextNumber());
        minVersion = gd.getNextNumber();
        return true;
    }
}
