/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import IAClasses.CrossCorrelation;
import IAClasses.ProgressDialog;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.ZProjector;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.FloatStatistics;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Rectangle;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

/**
 *
 * @author David Barry <david.barry at crick.ac.uk>
 */
public class TailAnalyser {

    protected DecimalFormat numFormat = new DecimalFormat("0.000");
//    private final int VEL_BINS = 5;
//    private final double VEL_BINS[] = {0.025, 0.05, 0.1, 0.15, 0.2};
    private final double VEL_BINS[] = {100.0};
//    private double sigSums[] = new double[100];
//    private int sigCounts[] = new int[100];
    ArrayList<Double> sigSums[] = new ArrayList[100];
    ImageStack sigStacks[] = new ImageStack[100];
    private static double spatialRes = 0.133333;
    int minSlices = 20;
    private static final int C1 = 0, C2 = 1;
    private static String protein = "Lifeact";
    private boolean randomize = false, tailByTail = false;
    private final String channelLabels[] = {"C0_Output", "C1_Output"};
    private final String checkBoxLabels[] = {"Randomize", "Tail by Tail"};
    private boolean checks[] = {randomize, tailByTail};
    private static int chanChoice = 1, eqChoice = 1, iterations = 1;
    double minVersion = 5.028;
    private boolean correlate = true;

//    public static void main(String args[]) {
//        TailAnalyser ta = new TailAnalyser();
//        ta.showDialog();
//        ta.run();
//        System.exit(0);
//    }

    public boolean showDialog() {
        GenericDialog gd = new GenericDialog("Tail Fitter");
        gd.addStringField("Protein Label:", protein, 5);
        gd.addRadioButtonGroup("Select Channel:", channelLabels, 1, 2, channelLabels[chanChoice]);
        gd.addChoice("Select Equation to Fit:", TailFitter.eqLabels, TailFitter.eqLabels[eqChoice]);
        gd.addCheckboxGroup(1, 2, checkBoxLabels, checks);
        gd.addNumericField("Iterations: ", iterations, 0);
        gd.addNumericField("Min Version: ", minVersion, 3);
        gd.showDialog();
        if (!gd.wasOKed()) {
            return false;
        }
        protein = gd.getNextString();
        chanChoice = gd.getNextRadioButton().equalsIgnoreCase(channelLabels[0]) ? TailAnalyser.C1 : TailAnalyser.C2;
        eqChoice = gd.getNextChoiceIndex();
        randomize = gd.getNextBoolean();
        tailByTail = gd.getNextBoolean();
        iterations = (int) Math.round(gd.getNextNumber());
        minVersion = gd.getNextNumber();
        return true;
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
        File resultsDir = new File(GenUtils.openResultsDirectory(parentDir + "/" + protein, "/"));
//        Arrays.fill(sigSums, 0.0);
//        Arrays.fill(sigCounts, 0);
        for (int i = 0; i < iterations; i++) {
            System.out.print(i + ",");
            if (tailByTail) {
                for (int j = 0; j < selectedSubDirs.size(); j++) {
                    File files[] = selectedSubDirs.get(j).listFiles();
                    ArrayList<ArrayList[]> sortedFiles = sortFiles(files, chan, j);
                    for (int n = 0; n < sortedFiles.size(); n++) {
                        ArrayList<FileWrapper> tailFiles[] = sortedFiles.get(n);
                        for (int m = 0; m < tailFiles.length; m++) {
                            System.out.print(j + ", " + n + "," + numFormat.format(VEL_BINS[m]) + " ");
                            ImageProcessor stackAverage = projectStack(buildTailAverage(stackWidth, stackHeight, tailFiles[m], chan, randomize));
                            if (stackAverage != null) {
                                processStackAverage(stackAverage);
                            }
                        }
                    }
                }
//            } else {
//                ImageStack stacks[] = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, randomize, resultsDir);
//                for (int j = 0; j < stacks.length; j++) {
//                    System.out.print("Vel:," + numFormat.format(VEL_BINS[j]) + ",");
//                    ImageProcessor stackAverage = projectStack(stacks[j]);
//                    processStackAverage(stackAverage);
//                }
            }
        }
        generateImages(selectedSubDirs, chan, resultsDir, stackWidth, stackHeight);
        if (correlate) {
            for (int i = 0; i < sigSums.length; i++) {
//                System.out.print((i / 100.0) + ",");
                System.out.print(i + ",");
                if (sigSums[i] != null) {
                    IJ.saveAs(new ImagePlus("", sigStacks[i]), "TIF", resultsDir + "/" + i);
                    int s = sigSums[i].size();
                    double data[] = new double[s];
                    for (int j = 0; j < s; j++) {
                        data[j] = sigSums[i].get(j);
                    }
                    Mean mean = new Mean();
                    StandardDeviation sd = new StandardDeviation();
                    System.out.print(mean.evaluate(data, 0, s) + "," + sd.evaluate(data, 0, s) + "," + s);
                }
                System.out.println();
            }
        }
    }

    TailFitter processStackAverage(ImageProcessor stackAverage) {
        Rectangle cropRoi = new Rectangle(0, 15, stackAverage.getWidth() - 1, 1);
        stackAverage.setRoi(cropRoi);
        stackAverage = stackAverage.crop();
        FloatStatistics stats = new FloatStatistics(stackAverage, ImageStatistics.MEDIAN + ImageStatistics.MIN_MAX, null);
        double median = stats.median;
        double min = stats.min;
        stackAverage.subtract(min);
        stackAverage.multiply(1.0 / (median - min));
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
        TailFitter tf = new TailFitter(eqChoice, spatialRes, Analyse_.SIG_EST_GREEN);
        tf.loadData(xVals, yVals, zVals);
        tf.doFit(Analyse_.SIG_EST_GREEN);
        return tf;
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

    ImageStack[] buildStackAverageOverall(int width, int height, ArrayList<File> subDirs, String channel, boolean randomize, File parentDir) {
        ImageStack output[] = new ImageStack[VEL_BINS.length];
        int nDirs = subDirs.size();
        Random r = new Random();
        ArrayList<FileWrapper> allFiles[] = new ArrayList[VEL_BINS.length];
        int nTails = 0;
        int imageCount = 0;
        for (int j = 0; j < nDirs; j++) {
            File files[] = subDirs.get(j).listFiles();
            ArrayList<ArrayList[]> sortedFiles = sortFiles(files, channel, j);
            int n = sortedFiles.size();
            nTails += n;
            ProgressDialog dpd = new ProgressDialog(null, "Correlating " + (j + 1) + " of " + nDirs, false, null, false);
            dpd.setVisible(true);
            for (int i = 0; i < n; i++) {
                dpd.updateProgress(i, n);
                ArrayList<FileWrapper> theseFiles[] = sortedFiles.get(i);
                for (int k = 0; k < VEL_BINS.length; k++) {
//                    if (correlate && parentDir != null) {
                    if (correlate) {
                        correlate(theseFiles[k]);
                    }
                    if (allFiles[k] == null) {
                        allFiles[k] = new ArrayList();
                    }
                    int m = theseFiles[k].size();
                    for (int l = 0; l < m; l++) {
                        allFiles[k].add(theseFiles[k].get(l));
                    }
                }
            }
            dpd.dispose();
        }
        for (int j = 0; j < VEL_BINS.length; j++) {
            output[j] = new ImageStack(width, height);
            int n = allFiles[j].size();
            ProgressDialog ipd = new ProgressDialog(null, "Building Stack " + (j + 1) + " of " + VEL_BINS.length, false, null, false);
            ipd.setVisible(true);
            for (int k = 0; k < n; k++) {
                ipd.updateProgress(k, n);
                int fileIndex = randomize ? r.nextInt(n) : k;
                ImagePlus imp = IJ.openImage(allFiles[j].get(fileIndex).getFile().getAbsolutePath());
                ImageProcessor ip = imp.getProcessor();
                FloatStatistics stats = new FloatStatistics(ip, ImageStatistics.MEDIAN, null);
                ip.multiply(1.0 / stats.median);
                output[j].addSlice(ip);
                imp.close();
            }
            ipd.dispose();
            imageCount += output[j].getSize();
            System.out.println("Images: " + output[j].getSize());
        }
        System.out.println("Cells: " + nDirs + " Tails: " + nTails + " Images: " + String.valueOf(imageCount));
        return output;
    }

    void buildCorrelationMaps(ArrayList<FileWrapper> files, int width, int height, File parentDir, int step, String label) {
        int s = files.size();
        if (s < minSlices) {
            return;
        }
        ImageProcessor originalImages[] = new ImageProcessor[s];
        for (int i = 0; i < s; i++) {
            originalImages[i] = IJ.openImage(files.get(i).getFile().getAbsolutePath()).getProcessor();
            (new GaussianBlur()).blurGaussian(originalImages[i], Analyse_.SIG_EST_GREEN, Analyse_.SIG_EST_GREEN, 0.01);
        }
        FloatProcessor accel = new FloatProcessor(s - 2 * step, step);
        FloatProcessor crossCorrelation = new FloatProcessor(width, height);
        for (int i = step; i < s - step; i++) {
            double vdiff = 0.0, tdiff = 0.0;
            for (int j = i - step; j < i + step; j++) {
                vdiff += files.get(j + 1).getVelocity() - files.get(j).getVelocity();
                tdiff += files.get(j + 1).getTimepoint() - files.get(j).getTimepoint();
            }
            accel.putPixelValue(i - step, 0, vdiff / tdiff);
        }
        double meanAccel = (new FloatStatistics(accel)).mean;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                FloatProcessor sigChange = new FloatProcessor(s - 2 * step, step);
                for (int i = step; i < s - step; i++) {
                    double sdiff = 0.0, tdiff = 0.0;
                    for (int j = i - step; j < i + step; j++) {
                        sdiff += originalImages[j + 1].getPixelValue(x, y) - originalImages[j].getPixelValue(x, y);
                        tdiff += files.get(j + 1).getTimepoint() - files.get(j).getTimepoint();
                    }
                    sigChange.putPixelValue(i - step, 0, sdiff / tdiff);
                }
                double meanSig = (new FloatStatistics(sigChange)).mean;
                crossCorrelation.putPixelValue(x, y, CrossCorrelation.crossCorrelation(accel, sigChange, 0, 0, meanAccel, meanSig));
            }
        }
        IJ.saveAs(new ImagePlus("", crossCorrelation), "TIF", parentDir + "/"
                + files.get(0).getDirIndex() + "_" + files.get(0).getIndex() + "_"
                + label + "_CC.tif");
    }

    void correlate(ArrayList<FileWrapper> files) {
        int s = files.size();
        if (s < minSlices) {
            return;
        }
        ImageProcessor originalImages[] = new ImageProcessor[s];
        double sigma = Analyse_.SIG_EST_GREEN / UserVariables.getSpatialRes();
        for (int i = 0; i < s; i++) {
            originalImages[i] = IJ.openImage(files.get(i).getFile().getAbsolutePath()).getProcessor();
            (new GaussianBlur()).blurGaussian(originalImages[i], sigma, sigma, 0.01);
            FloatStatistics stats = new FloatStatistics(originalImages[i], ImageStatistics.MEDIAN + ImageStatistics.MIN_MAX, null);
            double median = stats.median;
            double min = stats.min;
            originalImages[i].subtract(min);
            originalImages[i].multiply(1.0 / (median - min));
        }
        for (int i = 0; i < s; i++) {
            FloatProcessor fp = new FloatProcessor(originalImages[0].getWidth() - 1, 1);
//            TailFitter tf = processStackAverage(originalImages[i]);
//            double p[] = tf.getParams();
            double sum = 0.0;
            for (int x = 0; x < originalImages[i].getWidth() - 1; x++) {
//                double e = tf.evaluate(p, x);
                double e = originalImages[i].getPixelValue(x, 15);
                fp.putPixelValue(x, 0, e);
                sum += e;
            }
//            int index = 0;
//            while (files.get(i).getVelocity() > VEL_BINS[index] && index < VEL_BINS.length - 1) {
//                index++;
//            }
            int index = (int) Math.round(files.get(i).getVelocity() * 100.0);
            if (index < sigSums.length) {
                if (sigSums[index] == null) {
                    sigSums[index] = new ArrayList();
                    sigStacks[index] = new ImageStack(originalImages[i].getWidth() - 1, 1);
                }
                sigSums[index].add(sum);
                sigStacks[index].addSlice(fp);
//                sigSums[index] += sum;
//                sigCounts[index]++;
            }
//            System.out.println(files.get(i).getVelocity() + "," + sum);
        }
        return;
    }

//    void stackCorrelate(ImageStack originalImages[]) {
//        int s = originalImages.length;
//        for (int i = 0; i < s; i++) {
//            (new GaussianBlur()).blurGaussian(originalImages[i].getProcessor(i + 1), Analyse_.SIG_EST_GREEN, Analyse_.SIG_EST_GREEN, 0.01);
//        }
//        for (int i = 0; i < s; i++) {
//            FloatProcessor fp = new FloatProcessor(originalImages[0].getWidth() - 1, 1);
//            TailFitter tf = processStackAverage(originalImages[i].getProcessor(i + 1));
//            double p[] = tf.getParams();
//            double sum = 0.0;
//            for (int x = 0; x < originalImages[i].getWidth() - 1; x++) {
//                double e = tf.evaluate(p, x);
//                fp.putPixelValue(x, 0, e);
//                sum += e;
//            }
//            if (sigSums[i] == null) {
//                sigSums[i] = new ArrayList();
//                sigStacks[i] = new ImageStack(originalImages[i].getWidth() - 1, 1);
//            }
//            sigSums[i].add(sum);
//            sigStacks[i].addSlice(fp);
//        }
//    }

    ImageStack buildTailAverage(int stackwidth, int stackheight, ArrayList<FileWrapper> subDirs, String channel, boolean randomize) {
        ImageStack output = new ImageStack(stackwidth, stackheight);
        Random r = new Random();
        int n = subDirs.size();
        for (int k = 0; k < n; k++) {
            int fileIndex = randomize ? r.nextInt(n) : k;
            ImagePlus imp = IJ.openImage(subDirs.get(fileIndex).getFile().getAbsolutePath());
            ImageProcessor ip = imp.getProcessor();
            ImageStatistics stats = ip.getStatistics();
            double max = stats.max;
            double min = stats.min;
            ip.subtract(min);
            ip.multiply(1.0 / (max - min));
            output.addSlice(ip);
            imp.close();
        }
        return output;
    }

    ImageProcessor projectStack(ImageStack stack) {
        if (stack.getSize() > 0) {
            ZProjector zproj = new ZProjector(new ImagePlus("", stack));
            zproj.setMethod(ZProjector.AVG_METHOD);
            zproj.doProjection();
            return zproj.getProjection().getProcessor();
        } else {
            return null;
        }
    }

    ArrayList<ArrayList[]> sortFiles(File[] unsorted, String channel, int dirIndex) {
        ArrayList<ArrayList[]> files = new ArrayList();
        int n = unsorted.length;
        for (int i = 0; i < n; i++) {
            String name = FilenameUtils.getBaseName(unsorted[i].getName());
            Scanner scanner = new Scanner(name).useDelimiter("-");
            String thisChannel = scanner.next();
            if (channel.equalsIgnoreCase(thisChannel)) {
                int index = Integer.parseInt(scanner.next());
                double timepoint = Double.parseDouble(scanner.next());
                double velocity = Double.parseDouble(scanner.next());
                while (files.size() < index) {
                    files.add(new ArrayList[VEL_BINS.length]);
                    for (int j = 0; j < VEL_BINS.length; j++) {
                        files.get(files.size() - 1)[j] = new ArrayList();
                    }
                }
                int velBinIndex = 0;
                while (velocity > VEL_BINS[velBinIndex] && velBinIndex < VEL_BINS.length - 1) {
                    velBinIndex++;
                }
                FileWrapper sorted = new FileWrapper(unsorted[i], velocity, timepoint, index, dirIndex);
                files.get(index - 1)[velBinIndex].add(sorted);
            }
        }
        return files;
    }

    public void generateImages(ArrayList<File> selectedSubDirs, String chan, File directory, int stackWidth, int stackHeight) {
//        ImageStack stack;
//        if (tailByTail) {
//            stack = buildStackAverageTailByTail(stackWidth, stackHeight, selectedSubDirs, chan, false);
//        } else {
        ImageStack stack[] = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, false, null);
//        stackCorrelate(stack);
//        }
//        ImageStack stack = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, false);
        for (int i = 0; i < stack.length; i++) {
            ZProjector zproj = new ZProjector(new ImagePlus("", stack[i]));
            zproj.setMethod(ZProjector.AVG_METHOD);
            zproj.doProjection();
            ImageProcessor stackAverage = zproj.getProjection().getProcessor();
            String fileBaseName = directory + "/" + protein + "_" + TailFitter.eqLabels[eqChoice] + "_" + "Randomize-" + randomize + "_" + numFormat.format(VEL_BINS[i]) + "_" + channelLabels[chanChoice] + "_";
            IJ.saveAs((new ImagePlus("", stackAverage)), "text image", fileBaseName + "Average.txt");
            zproj.setMethod(ZProjector.SD_METHOD);
            zproj.doProjection();
            IJ.saveAs((new ImagePlus("", zproj.getProjection().getProcessor())), "text image", fileBaseName + "SD.txt");
//            Rectangle cropRoi = new Rectangle(0, 15, stackAverage.getWidth() - 2, 1);
//            stackAverage.setRoi(cropRoi);
//            stackAverage = stackAverage.crop();
//            ImageStatistics stats = stackAverage.getStatistics();
//            double max = stats.max;
//            double min = stats.min;
//            stackAverage.subtract(min);
//            stackAverage.multiply(1.0 / (max - min));
//            printImage(directory, fileBaseName);
//            IJ.saveAs((new ImagePlus("", stackAverage)), "text image", fileBaseName + "NormAverage.txt");
        }
    }

    private class FileWrapper {

        private final File file;
        private final double timepoint;
        private final double velocity;
        private final int index;
        private final int dirIndex;

        FileWrapper(File file, double velocity, double timepoint, int index, int dirIndex) {
            this.file = file;
            this.timepoint = timepoint;
            this.velocity = velocity;
            this.index = index;
            this.dirIndex = dirIndex;
        }

        public File getFile() {
            return file;
        }

        public double getTimepoint() {
            return timepoint;
        }

        public double getVelocity() {
            return velocity;
        }

        public int getIndex() {
            return index;
        }

        public int getDirIndex() {
            return dirIndex;
        }

    }

}
