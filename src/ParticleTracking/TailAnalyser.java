/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.ZProjector;
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

/**
 *
 * @author David Barry <david.barry at crick.ac.uk>
 */
public class TailAnalyser {

    protected DecimalFormat numFormat = new DecimalFormat("0.000");
    private final int VEL_BINS = 5;
    private final double VEL_BIN_LIMITS[] = {0.025, 0.05, 0.1, 0.15, 0.2};
    private static double spatialRes = 0.133333;
    int minSlices = 20;
    private static final int C1 = 0, C2 = 1;
    private static String protein = "Coronin";
    private boolean randomize = false, tailByTail = false;
    private final String channelLabels[] = {"C0_Output", "C1_Output"};
    private final String checkBoxLabels[] = {"Randomize", "Tail by Tail"};
    private boolean checks[] = {randomize, tailByTail};
    private static int chanChoice = 1, eqChoice = 1, iterations = 1;
    double minVersion = 5.022;

    public static void main(String args[]) {
        TailAnalyser ta = new TailAnalyser();
        ta.showDialog();
        ta.run();
        System.exit(0);
    }

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
        for (int i = 0; i < iterations; i++) {
            System.out.print(i + ",");
            if (tailByTail) {
                for (int j = 0; j < selectedSubDirs.size(); j++) {
                    File files[] = selectedSubDirs.get(j).listFiles();
                    ArrayList<ArrayList[]> sortedFiles = sortFiles(files, chan);
                    for (int n = 0; n < sortedFiles.size(); n++) {
                        ArrayList<File> tailFiles[] = sortedFiles.get(n);
                        for (int m = 0; m < tailFiles.length; m++) {
                            System.out.print(j + ", " + n + "," + numFormat.format(VEL_BIN_LIMITS[m]) + " ");
                            ImageProcessor stackAverage = projectStack(buildTailAverage(stackWidth, stackHeight, tailFiles[m], chan, randomize));
                            if (stackAverage != null) {
                                processStackAverage(stackAverage);
                            }
                        }
                    }
                }
            } else {
                ImageStack stacks[] = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, randomize);
                for (int j = 0; j < stacks.length; j++) {
                    System.out.print("Vel:," + numFormat.format(VEL_BIN_LIMITS[j]) + ",");
                    ImageProcessor stackAverage = projectStack(stacks[j]);
                    processStackAverage(stackAverage);
                }
            }
        }
        generateImages(selectedSubDirs, chan, parentDir, stackWidth, stackHeight);
    }

    void processStackAverage(ImageProcessor stackAverage) {
        Rectangle cropRoi = new Rectangle(0, 15, stackAverage.getWidth() - 1, 1);
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
        TailFitter tf = new TailFitter(eqChoice, spatialRes, Analyse_.SIG_EST_GREEN);
        tf.loadData(xVals, yVals, zVals);
        tf.doFit(Analyse_.SIG_EST_GREEN);
        tf.printParams();
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

    ImageStack[] buildStackAverageOverall(int stackwidth, int stackheight, ArrayList<File> subDirs, String channel, boolean randomize) {
        ImageStack output[] = new ImageStack[VEL_BINS];
        int nDirs = subDirs.size();
        Random r = new Random();
        ArrayList<File> allFiles[] = new ArrayList[VEL_BINS];
        int nTails = 0;
        int imageCount = 0;
        for (int j = 0; j < nDirs; j++) {
            File files[] = subDirs.get(j).listFiles();
            ArrayList<ArrayList[]> sortedFiles = sortFiles(files, channel);
            int n = sortedFiles.size();
            nTails += n;
            for (int i = 0; i < n; i++) {
                ArrayList<File> theseFiles[] = sortedFiles.get(i);
                for (int k = 0; k < VEL_BIN_LIMITS.length; k++) {
                    if (allFiles[k] == null) {
                        allFiles[k] = new ArrayList();
                    }
                    int m = theseFiles[k].size();
                    for (int l = 0; l < m; l++) {
                        allFiles[k].add(theseFiles[k].get(l));
                    }
                }
            }
        }
        for (int j = 0; j < VEL_BINS; j++) {
            output[j] = new ImageStack(stackwidth, stackheight);
            int n = allFiles[j].size();
            for (int k = 0; k < n; k++) {
                int fileIndex = randomize ? r.nextInt(n) : k;
                ImagePlus imp = IJ.openImage(allFiles[j].get(fileIndex).getAbsolutePath());
                ImageProcessor ip = imp.getProcessor();
                FloatStatistics stats = new FloatStatistics(ip, ImageStatistics.MEDIAN, null);
                ip.multiply(1.0 / stats.median);
                output[j].addSlice(ip);
                imp.close();
            }
            imageCount += output[j].getSize();
            System.out.println("Images: " + output[j].getSize());
        }
        System.out.println("Cells: " + nDirs + " Tails: " + nTails + " Images: " + String.valueOf(imageCount));
        return output;
    }

    ImageStack buildTailAverage(int stackwidth, int stackheight, ArrayList<File> subDirs, String channel, boolean randomize) {
        ImageStack output = new ImageStack(stackwidth, stackheight);
        Random r = new Random();
        int n = subDirs.size();
        for (int k = 0; k < n; k++) {
            int fileIndex = randomize ? r.nextInt(n) : k;
            ImagePlus imp = IJ.openImage(subDirs.get(fileIndex).getAbsolutePath());
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

    ArrayList<ArrayList[]> sortFiles(File[] unsorted, String channel) {
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
                    files.add(new ArrayList[VEL_BINS]);
                    for (int j = 0; j < VEL_BINS; j++) {
                        files.get(files.size() - 1)[j] = new ArrayList();
                    }
                }
                int velBinIndex = 0;
                while (velocity > VEL_BIN_LIMITS[velBinIndex] && velBinIndex < VEL_BINS - 1) {
                    velBinIndex++;
                }
                files.get(index - 1)[velBinIndex].add(unsorted[i]);
            }
        }
        return files;
    }

    public void generateImages(ArrayList<File> selectedSubDirs, String chan, File directory, int stackWidth, int stackHeight) {
//        ImageStack stack;
//        if (tailByTail) {
//            stack = buildStackAverageTailByTail(stackWidth, stackHeight, selectedSubDirs, chan, false);
//        } else {
        ImageStack stack[] = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, false);
//        }
//        ImageStack stack = buildStackAverageOverall(stackWidth, stackHeight, selectedSubDirs, chan, false);
        for (int i = 0; i < stack.length; i++) {
            ZProjector zproj = new ZProjector(new ImagePlus("", stack[i]));
            zproj.setMethod(ZProjector.AVG_METHOD);
            zproj.doProjection();
            ImageProcessor stackAverage = zproj.getProjection().getProcessor();
            String fileBaseName = directory + "/" + protein + "_" + TailFitter.eqLabels[eqChoice] + "_" + "Randomize-" + randomize + "_" + numFormat.format(VEL_BIN_LIMITS[i]) + "_" + channelLabels[chanChoice] + "_";
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

}
