/*
 * Copyright (C) 2017 Dave Barry <david.barry at crick.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package Particle_Analysis;

import Cell.Cell;
import Cell.CellRegion;
import Cell.Cytoplasm;
import Cell.Nucleus;
import IAClasses.Utils;
import static IJUtilities.IJUtils.hideAllImages;
import static IJUtilities.IJUtils.resetRoiManager;
import static IO.DataWriter.getAverageValues;
import static IO.DataWriter.convertArrayToString;
import static IO.DataWriter.saveTextWindow;
import static IO.DataWriter.saveValues;
import Math.Histogram;
import Particle.IsoGaussian;
import Particle.Particle;
import Particle.ParticleArray;
import ParticleTracking.UserVariables;
import Revision.Revision;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.OpenDialog;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.EDM;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.frame.RoiManager;
import ij.process.Blitter;
import ij.process.ByteBlitter;
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.ShortProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import ui.DetectionGUI;

public class Particle_Mapper extends Particle_Tracker {

    private static double histMin = -5.0, histMax = 20.0, threshLevel = 50.0;
    private static boolean useThresh = true, isolateFoci = true, analyseFluorescence = true,
            averageImage = false;
    private static int histNBins = 40;
    String resultsDir;
    private Cell[] cells;
    private final int NUCLEI = 2, CYTO = 1, FOCI = 0;
    private final String NUCLEI_MASK = "Nuclei Mask", FLUO_DIST = "fluorescence_distribution_data.csv",
            FOCI_DIST = "foci_distance_data.csv", FOCI_DETECTIONS = "Foci Detections", DIST_MAP = "Distance Map",
            FOCI_DIST_HIST = "foci_distance_histogram.csv", FOCI_NUC_ASS = "Foci-Nuclei Associations",
            CELL_BOUNDS = "Cell Boundaries", FOCI_DIST_DATA = "Mean Foci Data";
    private final String FLUO_HEADINGS[] = new String[]{"Cell ID", "Cell Mean", "Cell Std Dev",
        "Nuclear Mean", "Nuclear Std Dev", "Cytosolic Mean",
        "Cytosolic Std Dev", "Nuclear Mean / Cytosolic Mean",
        "Nuclear Std Dev / Cytosolic Std Dev"};
    /**
     * Title of the application
     */
    protected String title = "Particle Mapper";
    private final String resultsHeadings = "Cell ID\t X\t Y\t Number of Foci\t Mean Intensity of Foci\t Mean Foci Distance To Nuclear Boundary (" + IJ.micronSymbol + "m)";

    /**
     * Default constructor.
     */
    public Particle_Mapper() {

    }

    public void run(String arg) {
        Prefs.blackBackground = false;
        title = title + "_v" + Revision.VERSION + "." + intFormat.format(Revision.revisionNumber);
        inputs = new ImagePlus[3];
        if (IJ.getInstance() == null) {
            inputs[NUCLEI] = IJ.openImage((new OpenDialog("Specify Nuclei Image", null)).getPath());
            inputs[FOCI] = IJ.openImage((new OpenDialog("Specify Foci Image", null)).getPath());
            inputs[CYTO] = IJ.openImage((new OpenDialog("Specify Image For Thresholding", null)).getPath());
        } else {
            int[] idList = WindowManager.getIDList();
            if (idList == null) {
                GenUtils.error("No Images Open.");
                return;
            }
            for (int i = 0; i < idList.length; i++) {
                inputs[i] = WindowManager.getImage(idList[i]);
            }
        }
        readParamsFromImage();
        if (!showDialog()) {
            return;
        }
        ImageStack[] stacks = getStacks();
        for (ImageStack a : stacks) {
            if (a != null) {
                for (ImageStack b : stacks) {
                    if (b != null) {
                        if (a.size() != b.size()) {
                            GenUtils.error("All stacks must have same number of slices.");
                            return;
                        }
                    }
                }
            }
        }
        hideAllImages();
        File inputDir = buildStacks();
        if (inputDir == null) {
            return;
        }
        if (isolateFoci && !showDetectionGui()) {
            return;
        }
        TextWindow tw = null;
        boolean hideOutputs = stacks[0].size() > 1;
        if (analyseFluorescence && averageImage) {
            tw = new TextWindow("Average Fluorescence Distributions", convertArrayToString("N\t", FLUO_HEADINGS, "\t"), new String(), 640, 480);
        }
        try {
            resultsDir = GenUtils.openResultsDirectory(Utilities.getFolder(inputDir,
                    "Specify directory for output files...", true) + delimiter + title);
            if (resultsDir == null) {
                GenUtils.error("Failed to create output directory.");
                return;
            }
            for (int i = 1; i <= stacks[0].size(); i++) {
                File thisDir = GenUtils.createDirectory(String.format("%s%sSlice_%d", resultsDir, File.separator, i), true);
                ByteProcessor binaryNuclei = checkBinaryImage(stacks[NUCLEI].getProcessor(i));
                IJ.saveAs(new ImagePlus("", binaryNuclei), "PNG", String.format("%s%s%s", thisDir.getAbsolutePath(), File.separator, NUCLEI_MASK));
                if (!findCells(binaryNuclei.duplicate())) {
                    IJ.log(String.format("No cells found in image %d.", i));
                } else {
                    ImageProcessor cellMap = buildTerritories(binaryNuclei.duplicate(), thisDir.getAbsolutePath()).getProcessor();
                    Arrays.sort(cells);
                    if (useThresh) {
                        filterCells(stacks[CYTO].getProcessor(i), new Cytoplasm(), threshLevel, Measurements.MEAN);
                    }
                    if (isolateFoci) {
                        ParticleArray pa = findParticles();
                        assignParticlesToCells(pa, cellMap, thisDir.getAbsolutePath());
                        drawDetections(pa, stacks[FOCI].getWidth(), stacks[FOCI].getHeight(), thisDir.getAbsolutePath());
                        double[][] distances = calcDistances(buildDistanceMap(binaryNuclei, thisDir.getAbsolutePath()));
                        buildHistograms(distances, histNBins, histMax, histMin, thisDir.getAbsoluteFile(), hideOutputs);
                        outputFociDistanceData(distances, thisDir.getAbsolutePath(), resultsHeadings, hideOutputs);
                    }
                    if (analyseFluorescence) {
                        double[][] vals = analyseCellFluorescenceDistribution(stacks[FOCI].getProcessor(i),
                                Measurements.MEAN + Measurements.STD_DEV);
                        String outputFileName = String.format("%s%s%s", thisDir.getAbsolutePath(), File.separator, FLUO_DIST);
                        saveValues(vals, new File(outputFileName), FLUO_HEADINGS);
                        if (averageImage) {
                            tw.append(convertArrayToString(null, getAverageValues(vals, FLUO_HEADINGS.length), "\t"));
                        }
                    }
                }
            }
            if (tw != null) {
                saveTextWindow(tw, new File(String.format("%s%s%s", resultsDir, File.separator, "Mean Image Data.csv")), convertArrayToString("N\t", FLUO_HEADINGS, "\t"));
            }
            cleanUp();
        } catch (Exception e) {
            GenUtils.error(e.getMessage());
        }
        for (ImagePlus imp : inputs) {
            if (imp != null) {
                imp.close();
            }
        }
        IJ.showStatus(String.format("%s done.", title));
    }

    /**
     * Uses {@link ij.plugin.filter.ParticleAnalyzer ParticleAnalyzer} to find
     * centroids of objects in <code>image</code>.
     *
     * @param image a binary image
     * @return a <code>2 x n</code> array, specifying the centroids of any
     * objects found in <code>image</code>. The n<sup>th</sup> centroid, in the
     * form (x, y), is accessed as (centroids[0][n], centroids[1][n]).
     */
    public boolean findCells(ImageProcessor image) {
        ImagePlus imp = new ImagePlus("", image);
        ResultsTable rt = Analyzer.getResultsTable();
        resetRoiManager();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE + ParticleAnalyzer.CLEAR_WORKSHEET + ParticleAnalyzer.ADD_TO_MANAGER, Measurements.CENTROID, rt, 0, Double.MAX_VALUE);
        pa.setHideOutputImage(true);
        pa.analyze(imp);
        int n = rt.size();
        if (!(n > 0)) {
            return false;
        }
        cells = new Cell[n];
        RoiManager roimanager = RoiManager.getInstance();
        roimanager.setVisible(false);
        Roi[] rois = roimanager.getRoisAsArray();
        int x = rt.getColumnIndex("X");
        int y = rt.getColumnIndex("Y");
        for (int i = 0; i < n; i++) {
            double[] centroid = new double[]{rt.getValueAsDouble(x, i), rt.getValueAsDouble(y, i)};
            cells[i] = (new Cell(new Nucleus(rois[i], centroid)));
        }
        return true;
    }

    /**
     * Creates an indexed list of regions, based on a voronoi segmentation of
     * objects in <code>image</code>.
     *
     * @param image greyscale image in which the constructed regions each have a
     * unique label
     * @return
     */
    public ImagePlus buildTerritories(ImageProcessor image, String resultsDir) {
        ImagePlus imp = new ImagePlus("", image);
        ResultsTable rt = Analyzer.getResultsTable();
        EDM edm = new EDM();
        edm.setup("voronoi", imp);
        edm.run(image);
        image.threshold(1);
        image.setColor(Color.white);
        int fontsize = (int) Math.round(0.05 * Math.min(image.getWidth(), image.getHeight()));
        Font font = new Font("Times", Font.BOLD, fontsize);
        image.setFont(font);
        rt.reset();
        resetRoiManager();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_ROI_MASKS + ParticleAnalyzer.ADD_TO_MANAGER, 0, rt, 0, Double.MAX_VALUE);
        pa.setHideOutputImage(true);
        pa.analyze(imp);
        int n = cells.length;
        RoiManager roimanager = RoiManager.getInstance();
        roimanager.setVisible(false);
        Roi[] rois = roimanager.getRoisAsArray();
        ImageProcessor cellMap = pa.getOutputImage().getProcessor();
        int duds = -1;
        for (int i = 0; i < n; i++) {
            Cell c = cells[i];
            double[] centroid = c.getNucleus().getCentroid();
            int xc = (int) Math.round(centroid[0]);
            int yc = (int) Math.round(centroid[1]);
            int id = cellMap.getPixel(xc, yc);
            if (id < 1 || cellExists(new Cell(id))) {
                id = duds--;
                IJ.log(String.format("The object detected at (%d, %d) could not be"
                        + " assigned to a unique cell. This could indicate a multi-nucleate"
                        + " cell, but could also suggest an excessively noisy image."
                        + " Consider preprocessing your nuclei image to reduce noise.", xc, yc));
            }
            c.setID(id);
            if (id > 0) {
                c.addCellRegion(new Cytoplasm(rois[id - 1]));
            }
            image.drawString(String.valueOf(id), xc, yc);
        }
        IJ.saveAs(new ImagePlus("", image), "PNG", String.format("%s%s%s", resultsDir, File.separator, CELL_BOUNDS));
        return pa.getOutputImage();
    }

    boolean cellExists(Cell c) {
        Cell[] cellsCopy = new Cell[cells.length];
        System.arraycopy(cells, 0, cellsCopy, 0, cells.length);
        Arrays.sort(cellsCopy);
        return (Arrays.binarySearch(cellsCopy, c) > 0);
    }

    /**
     *
     * @param pa
     * @param cellMap
     * @param cellCentroids
     * @return
     */
    void assignParticlesToCells(ParticleArray pa, ImageProcessor cellMap, String resultsDir) {
        ByteProcessor map = new ByteProcessor(cellMap.getWidth(), cellMap.getHeight());
        map.setValue(0);
        map.fill();
        map.setValue(255);
        ArrayList<Particle> detections = pa.getLevel(0);
        int N = detections.size();
        LinkedHashMap<Integer, Integer> idToIndexMap = new LinkedHashMap();
        for (int i = 0; i < cells.length; i++) {
            idToIndexMap.put(cells[i].getID(), i);
        }
        for (int i = 0; i < N; i++) {
            Particle p = detections.get(i);
            int xp = (int) Math.round(p.getX() / UserVariables.getSpatialRes());
            int yp = (int) Math.round(p.getY() / UserVariables.getSpatialRes());
            Integer pIndex = idToIndexMap.get(cellMap.getPixel(xp, yp));
            if (pIndex != null) {
                Cell c = cells[pIndex];
                double[] centroid = c.getNucleus().getCentroid();
                map.drawLine(xp, yp,
                        (int) Math.round(centroid[0]), (int) Math.round(centroid[1]));
                c.addParticle(p);
            }
        }
        IJ.saveAs(new ImagePlus("", map), "PNG", String.format("%s%s%s", resultsDir, File.separator, FOCI_NUC_ASS));
    }

    /**
     *
     * @param assignments
     * @param distanceMap
     * @return
     */
    double[][] calcDistances(ImageProcessor distanceMap) {
        int N = cells.length;
        double[][] distances = new double[N][];
        double res = UserVariables.getSpatialRes();
        for (int i = 0; i < N; i++) {
            Cell c = cells[i];
            ArrayList<Particle> particles = c.getParticles();
            if (particles != null) {
                int M = particles.size();
                distances[i] = new double[M];
                for (int j = 0; j < M; j++) {
                    Particle p = particles.get(j);
                    distances[i][j] = res * distanceMap.getInterpolatedValue(p.getX() / res, p.getY() / res);
                }
            }
        }
        return distances;
    }

    /**
     *
     * @param distances
     * @param nBins
     * @param max
     * @param min
     */
    void buildHistograms(double[][] distances, int nBins, double max, double min, File resultsDir, boolean hidePlot) {
        int N = distances.length;
        int[][] histograms = new int[N][nBins];
        float[] meanHistogram = new float[nBins];
        float[] bins = new float[nBins];
        float binSize = (float) (max - min) / nBins;
        for (int i = 0; i < N; i++) {
            if (distances[i] != null) {
                histograms[i] = Histogram.calcHistogram(distances[i], min, max, nBins);
            }
        }
        for (int i = 0; i < nBins; i++) {
            DescriptiveStatistics mean = new DescriptiveStatistics();
            for (int j = 0; j < N; j++) {
                if (distances[j] != null) {
                    mean.addValue((double) 100.0 * histograms[j][i] / distances[j].length);
                }
            }
            meanHistogram[i] = (float) mean.getMean();
            bins[i] = (float) (i * binSize + min);
        }
        Plot histPlot = new Plot("Mean Foci Distance Distribution",
                "Distance (" + IJ.micronSymbol + "m)", "% Frequency");
        histPlot.setLineWidth(3);
        histPlot.setColor(Color.blue);
        histPlot.addPoints(bins, meanHistogram, Plot.CONNECTED_CIRCLES);
        histPlot.draw();
        if (!hidePlot) {
            histPlot.show();
        }
        ResultsTable rt = histPlot.getResultsTable();
        rt.save(String.format("%s%s%s", resultsDir, File.separator, FOCI_DIST_HIST));
    }

    /**
     *
     * @param image
     * @return
     */
    ImageProcessor buildDistanceMap(ImageProcessor image, String resultsDir) {
        image.invert();
        ImageProcessor invertedImage = image.duplicate();
        invertedImage.invert();
        EDM edm = new EDM();
        FloatProcessor fgDistanceMap = edm.makeFloatEDM(image, 0, false);
        FloatProcessor distanceMap = edm.makeFloatEDM(invertedImage, 0, false);
        (new FloatBlitter(distanceMap)).copyBits(fgDistanceMap, 0, 0, Blitter.SUBTRACT);
        IJ.saveAs(new ImagePlus("", distanceMap), "TIF", String.format("%s%s%s", resultsDir, File.separator, DIST_MAP));
        return distanceMap;
    }

    /**
     * Displays the user interface
     *
     * @return true if the user click OK, false otherwise
     */
    public boolean showDialog() {
        String[] imageTitles;
        if (IJ.getInstance() != null) {
            imageTitles = WindowManager.getImageTitles();
        } else {
            imageTitles = new String[3];
            imageTitles[NUCLEI] = inputs[NUCLEI] != null ? inputs[NUCLEI].getTitle() : " ";
            imageTitles[FOCI] = inputs[FOCI] != null ? inputs[FOCI].getTitle() : " ";
            imageTitles[CYTO] = inputs[CYTO] != null ? inputs[CYTO].getTitle() : " ";
        }
        int N = imageTitles.length;
        if (N < 2) {
            IJ.error("Minimum of two images required.");
            return false;
        }
        GenericDialog gd = new GenericDialog(title);
        Font bFont = gd.getFont();
        if (bFont == null) {
            bFont = new Font("Times", Font.BOLD, 12);
        } else {
            bFont = bFont.deriveFont(Font.BOLD);
        }
        gd.centerDialog(true);
        gd.addChoice("Nuclei: ", imageTitles, imageTitles[NUCLEI < N ? NUCLEI : 0]);
        gd.addChoice("Protein distribution to be quantified: ", imageTitles, imageTitles[FOCI < N ? FOCI : 0]);
        gd.addMessage("Do you want to use a third image to select cells based on intensity threshold?", bFont);
        gd.addCheckbox("Use threshold image", useThresh);
        gd.addChoice("Select threshold image: ", imageTitles, imageTitles[CYTO < N ? CYTO : 0]);
        gd.addSlider("Specify threshold level %", 0.0, 100.0, threshLevel);
        gd.addMessage("How do you want to analyse the protein distribution?", bFont);
        String[] checkBoxLabels = new String[]{"Attempt to isolate foci", "Quantify entire distribution"};
        gd.addCheckboxGroup(1, 2, checkBoxLabels, new boolean[]{isolateFoci, analyseFluorescence});
        gd.addMessage("Specify ranges and bin size for distance histogram", bFont);
        gd.addNumericField("Minimum Value:", histMin, 1);
        gd.addNumericField("Maximum Value:", histMax, 1);
        gd.addNumericField("Number of Bins:", histNBins, 0);
        gd.addMessage("How do you want results to be output?", bFont);
        String[] radioButtonLabels = new String[]{"Show me data for each cell", "Summarise data for each image"};
        gd.addRadioButtonGroup(null, radioButtonLabels, 1, 2, radioButtonLabels[averageImage ? 1 : 0]);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        boolean c1, c2, c3;
        int choice1 = gd.getNextChoiceIndex(), choice2 = gd.getNextChoiceIndex(), choice3 = gd.getNextChoiceIndex();
        if (IJ.getInstance() == null) {
            c1 = inputs[choice1] != null;
            c2 = inputs[choice2] != null;
            c3 = inputs[choice3] != null;
        } else {
            c1 = WindowManager.getImage(imageTitles[choice1]) != null;
            c2 = WindowManager.getImage(imageTitles[choice2]) != null;
            c3 = WindowManager.getImage(imageTitles[choice3]) != null;
        }
        useThresh = gd.getNextBoolean();
        threshLevel = gd.getNextNumber();
        isolateFoci = gd.getNextBoolean();
        analyseFluorescence = gd.getNextBoolean();
        histMin = gd.getNextNumber();
        histMax = gd.getNextNumber();
        histNBins = (int) Math.round(gd.getNextNumber());
        averageImage = gd.getNextRadioButton().equals(radioButtonLabels[1]);
        if (!c1) {
            GenUtils.error("You have not specified a nuclei image.");
            return showDialog();
        } else if (!c2) {
            GenUtils.error("You have not specified a protein distribution image.");
            return showDialog();
        } else if (!c3 && useThresh) {
            GenUtils.error("You have not specified an image for thresholding.");
            return showDialog();
        } else {
            if (IJ.getInstance() == null) {
                ImagePlus[] inputsCopy = new ImagePlus[3];
                inputsCopy[NUCLEI] = inputs[NUCLEI] != null ? inputs[NUCLEI].duplicate() : null;
                inputsCopy[FOCI] = inputs[FOCI] != null ? inputs[FOCI].duplicate() : null;
                inputsCopy[CYTO] = inputs[CYTO] != null ? inputs[CYTO].duplicate() : null;
                inputs[NUCLEI] = inputs[choice1] != null ? inputsCopy[choice1].duplicate() : null;
                inputs[FOCI] = inputs[choice2] != null ? inputsCopy[choice2].duplicate() : null;
                inputs[CYTO] = inputs[choice3] != null ? inputsCopy[choice3].duplicate() : null;
            } else {
                inputs[NUCLEI] = WindowManager.getImage(imageTitles[choice1]);
                inputs[FOCI] = WindowManager.getImage(imageTitles[choice2]);
                inputs[CYTO] = WindowManager.getImage(imageTitles[choice3]);
            }
            return true;
        }
    }

    boolean showDetectionGui() {
        DetectionGUI ui = new DetectionGUI(null, true, title, this);
        ui.setVisible(true);
        return ui.isWasOKed();
    }

    /**
     * Generates an image of the particles contained in <code>pa</code>.
     *
     * @param pa array of particles to be drawn
     * @param width width of output image
     * @param height height of output image
     */
    public void drawDetections(ParticleArray pa, int width, int height, String resultsDir) {
        int depth = pa.getDepth();
        for (int d = 0; d < depth; d++) {
            ArrayList<Particle> level = pa.getLevel(d);
            ShortProcessor output = new ShortProcessor(width, height);
            for (Particle p : level) {
                if (p instanceof IsoGaussian) {
                    Utils.draw2DGaussian(output, (IsoGaussian) p, 0.0, UserVariables.getSpatialRes(), false);
                }
            }
            output.multiply(1.0 / normFactor);
            IJ.saveAs(new ImagePlus("", output), "TIF", String.format("%s%s%s", resultsDir, File.separator, FOCI_DETECTIONS));
        }
    }

    /**
     * Duplicates the image containing foci, saves a normalised copy and returns
     * a reference to the saved location
     *
     * @return absolute path to the normalised images
     */
    protected String prepareInputs() {
        if (inputs[FOCI] == null) {
            return null;
        }
        return normaliseStacks(inputs[FOCI].getImageStack(), null);
    }

    /**
     *
     * @param distances
     * @param centroids
     * @param assignments
     * @throws IOException
     */
    void outputFociDistanceData(double[][] distances, String resultsDir, String resultsHeadings, boolean hideTextWindow) throws IOException {
        int N = distances.length;
        TextWindow tw = new TextWindow("Mean Foci Data", resultsHeadings, new String(), 640, 480);
        DecimalFormat df1 = new DecimalFormat("00");
        DecimalFormat df2 = new DecimalFormat("0.000");
        for (int i = 0; i < N; i++) {
            Cell c = cells[i];
            double[] centroid = c.getNucleus().getCentroid();
            ArrayList<Particle> particles = c.getParticles();
            String result = c.getID() + "\t" + df1.format(Math.round(centroid[0])) + "\t" + df1.format(Math.round(centroid[1]));
            if (particles != null) {
                int L = particles.size();
                DescriptiveStatistics intensitites = new DescriptiveStatistics();
                DescriptiveStatistics dist = new DescriptiveStatistics();
                for (int j = 0; j < L; j++) {
                    Particle p = particles.get(j);
                    intensitites.addValue(p.getMagnitude());
                    dist.addValue(distances[i][j]);
                }
                result = result.concat("\t" + L + "\t" + df2.format(intensitites.getMean()) + "\t" + df2.format(dist.getMean()));
            } else {
                result = result.concat("\t 0\t " + df2.format(0) + "\t " + df2.format(0));
            }
            tw.append(result);
        }
        tw.setVisible(!hideTextWindow);
        saveTextWindow(tw, new File(String.format("%s%s%s", resultsDir, File.separator, FOCI_DIST)), resultsHeadings);
    }

    /**
     *
     * @param binaryImage
     */
    ByteProcessor checkBinaryImage(ImageProcessor image) {
        image.resetMinAndMax();
        ByteProcessor binaryImage = (ByteProcessor) (new TypeConverter(image, true)).convertToByte();
        binaryImage.autoThreshold();
        if (binaryImage.isInvertedLut()) {
            binaryImage.invertLut();
        }
        ImageStatistics stats = binaryImage.getStatistics();
        if (stats.histogram[0] > stats.histogram[255]) {
            binaryImage.invert();
        }
        return binaryImage;
    }

    void filterCells(ImageProcessor image, CellRegion regionType, double threshold, int measurement) {
        DescriptiveStatistics ds = new DescriptiveStatistics();
        boolean[] selected = new boolean[cells.length];
        Arrays.fill(selected, false);
        int b = 0;
        for (Cell cell : cells) {
            CellRegion cr = cell.getRegion(regionType);
            if (cr != null) {
                selected[b] = true;
                image.setRoi(cr.getRoi());
                ImageStatistics stats = ImageStatistics.getStatistics(image, measurement, null);
                switch (measurement) {
                    case Measurements.MEAN:
                        ds.addValue(stats.mean);
                        break;
                    case Measurements.STD_DEV:
                        ds.addValue(stats.stdDev);
                        break;
                    default:
                        ds.addValue(0.0);
                }
            }
            b++;
        }
        double percentile = threshold > 0.0 ? ds.getPercentile(threshold) : 0.0;
        double[] measures = ds.getValues();
        ArrayList<Cell> cells2 = new ArrayList();
        for (int i = 0, j = 0; i < cells.length; i++) {
            if (selected[i] && measures[j++] > percentile) {
                cells2.add(cells[i]);
            }
        }
        cells = cells2.toArray(new Cell[]{});
    }

    double[][] analyseCellFluorescenceDistribution(ImageProcessor image, int measurements) {
        int N = cells.length;
        double[][] vals = new double[N][];
        for (int i = 0; i < N; i++) {
            Cell cell = cells[i];
            if (cell.getID() > 0) {
                Nucleus nucleus = (Nucleus) cell.getRegion(new Nucleus());
                Cytoplasm cyto = (Cytoplasm) cell.getRegion(new Cytoplasm());
                Roi nucRoi = nucleus.getRoi();
                ByteProcessor nucMask = (ByteProcessor) nucRoi.getMask();
                image.setRoi(nucRoi);
                image.setMask(nucMask);
                ImageStatistics nucstats = ImageStatistics.getStatistics(image, measurements, null);

                Roi cytoRoi = cyto.getRoi();
                ByteProcessor cytoMask = (ByteProcessor) cytoRoi.getMask();
                image.setRoi(cytoRoi);
                image.setMask(cytoMask);
                ImageStatistics cellstats = ImageStatistics.getStatistics(image, measurements, null);

                int xc = nucRoi.getBounds().x - cytoRoi.getBounds().x;
                int yc = nucRoi.getBounds().y - cytoRoi.getBounds().y;
                (new ByteBlitter(cytoMask)).copyBits(nucMask, xc, yc, Blitter.SUBTRACT);
                image.setRoi(cytoRoi);
                image.setMask(cytoMask);
                ImageStatistics cytostats = ImageStatistics.getStatistics(image, measurements, null);

                vals[i] = new double[]{cell.getID(), cellstats.mean, cellstats.stdDev,
                    nucstats.mean, nucstats.stdDev,
                    cytostats.mean, cytostats.stdDev,
                    nucstats.mean / cytostats.mean,
                    nucstats.stdDev / cytostats.stdDev};
            }
        }
        return vals;
    }
}
