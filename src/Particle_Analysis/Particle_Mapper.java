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
import static IJUtilities.IJUtils.closeAllOpenImages;
import static IJUtilities.IJUtils.hideAllImages;
import static IJUtilities.IJUtils.resetRoiManager;
import static IO.DataWriter.saveTextWindow;
import static IO.DataWriter.saveValues;
import Math.Histogram;
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

    private double histMin = -5.0, histMax = 20.0, threshLevel = 50.0;
    private boolean useThresh = true, isolateFoci = true, analyseFluorescence = true;
    private int histNBins = 40;
    String outputDirName;
    private Cell[] cells;
    ImagePlus[] inputs = new ImagePlus[3];
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
        stacks = new ImageStack[2];
        if (IJ.getInstance() == null) {
            inputs[0] = IJ.openImage((new OpenDialog("Specify Nuclei Image", null)).getPath());
            inputs[1] = IJ.openImage((new OpenDialog("Specify Foci Image", null)).getPath());
            inputs[2] = IJ.openImage((new OpenDialog("Specify Image For Thresholding", null)).getPath());
        }
        if (!showDialog()) {
            return;
        }
        hideAllImages();
        File inputDir = buildStacks();
        if (inputDir == null) {
            return;
        }
        if (isolateFoci && !showDetectionGui()) {
            return;
        }
        outputDirName = GenUtils.openResultsDirectory(Utilities.getFolder(inputDir,
                "Specify directory for output files...", true) + delimiter + title);
        try {
            ByteProcessor binaryNuclei = checkBinaryImage(inputs[0].getProcessor());
            IJ.saveAs(new ImagePlus("", binaryNuclei), "PNG", outputDirName + "/Nuclei Mask");
            findCells(binaryNuclei.duplicate());
            ImageProcessor cellMap = buildTerritories(binaryNuclei.duplicate()).getProcessor();
            Arrays.sort(cells);
            if (useThresh) {
                filterCells(inputs[2].getProcessor(), new Cytoplasm(), threshLevel, Measurements.MEAN);
            }
            if (isolateFoci) {
                ParticleArray pa = findParticles();
                assignParticlesToCells(pa, cellMap);
                drawDetections(pa, stacks[0].getWidth(), stacks[0].getHeight());
                double[][] distances = calcDistances(buildDistanceMap(binaryNuclei));
                buildHistograms(distances, histNBins, histMax, histMin);
                outputFociDistanceData(distances);
            }
            if (analyseFluorescence) {
                saveValues(analyseCellFluorescenceDistribution(stacks[0].getProcessor(1),
                        Measurements.MEAN + Measurements.STD_DEV),
                        new File(outputDirName + "/fluorescence_distribution_data.csv"),
                        new String[]{"Cell ID", "Nuclear Mean", "Nuclear Std Dev", "Cytosolic Mean",
                            "Cytosolic Std Dev", "Nuclear Mean / Cytosolic Mean",
                            "Nuclear Std Dev / Cytosolic Std Dev"});
            }
            cleanUp();
        } catch (IOException e) {
            IJ.error(e.getMessage());
        }
        closeAllOpenImages();
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
    public void findCells(ImageProcessor image) {
        ImagePlus imp = new ImagePlus("", image);
        ResultsTable rt = Analyzer.getResultsTable();
        resetRoiManager();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE + ParticleAnalyzer.CLEAR_WORKSHEET + ParticleAnalyzer.ADD_TO_MANAGER, Measurements.CENTROID, rt, 0, Double.MAX_VALUE);
        pa.setHideOutputImage(true);
        pa.analyze(imp);
        int n = rt.size();
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
    }

    /**
     * Creates an indexed list of regions, based on a voronoi segmentation of
     * objects in <code>image</code>.
     *
     * @param image greyscale image in which the constructed regions each have a
     * unique label
     * @return
     */
    public ImagePlus buildTerritories(ImageProcessor image) {
        ImagePlus imp = new ImagePlus("", image);
        ResultsTable rt = Analyzer.getResultsTable();
        EDM edm = new EDM();
        edm.setup("voronoi", imp);
        edm.run(image);
        image.threshold(1);
        image.setColor(Color.white);
        int fontsize = (int)Math.round(0.05 * Math.min(image.getWidth(), image.getHeight()));
        Font font = new Font("Times",Font.BOLD,fontsize);
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
                IJ.log("The partition of the field of view into distinct cell bodies "
                        + "was not completely successful, possibly due to the nuclei "
                        + "image being excessively noisy.");
            }
            c.setID(id);
            if (id > 0) {
                c.addCellRegion(new Cytoplasm(rois[id - 1]));
            }
            image.drawString(String.valueOf(id), xc, yc);
        }
        IJ.saveAs(new ImagePlus("", image), "PNG", outputDirName + "/Cell Boundaries");
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
    void assignParticlesToCells(ParticleArray pa, ImageProcessor cellMap) {
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
        IJ.saveAs(new ImagePlus("", map), "PNG", outputDirName + "/Foci-Nuclei Associations");
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
    void buildHistograms(double[][] distances, int nBins, double max, double min) {
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
        histPlot.show();
        ResultsTable rt = histPlot.getResultsTable();
        rt.save(outputDirName + "/foci_distance_histogram.csv");
    }

    /**
     *
     * @param image
     * @return
     */
    ImageProcessor buildDistanceMap(ImageProcessor image) {
        image.invert();
        ImageProcessor invertedImage = image.duplicate();
        invertedImage.invert();
        EDM edm = new EDM();
        FloatProcessor fgDistanceMap = edm.makeFloatEDM(image, 0, false);
        FloatProcessor distanceMap = edm.makeFloatEDM(invertedImage, 0, false);
        (new FloatBlitter(distanceMap)).copyBits(fgDistanceMap, 0, 0, Blitter.SUBTRACT);
        IJ.saveAs(new ImagePlus("", distanceMap), "TIF", outputDirName + "/Distance Map");
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
            imageTitles = new String[]{inputs[0].getTitle(), inputs[1].getTitle(), inputs[2].getTitle()};
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
        gd.addChoice("Nuclei: ", imageTitles, imageTitles[0]);
        gd.addChoice("Protein distribution to be quantified: ", imageTitles, imageTitles[1]);
        gd.addMessage("Do you want to use a third image to select cells based on intensity threshold?", bFont);
        gd.addCheckbox("Use threshold image", useThresh);
        gd.addChoice("Select threshold image: ", imageTitles, imageTitles[N > 2 ? 2 : 1]);
        gd.addSlider("Specify threshold level %", 0.0, 100.0, threshLevel);
        gd.addMessage("How do you want to analyse the protein distribution?", bFont);
        String[] checkBoxLabels = new String[]{"Attempt to isolate foci", "Just quantify the distribution"};
        gd.addCheckboxGroup(1, 2, checkBoxLabels, new boolean[]{isolateFoci, analyseFluorescence});
        gd.addMessage("Specify ranges and bin size for distance histogram", bFont);
        gd.addNumericField("Minimum Value:", histMin, 1);
        gd.addNumericField("Maximum Value:", histMax, 1);
        gd.addNumericField("Number of Bins:", histNBins, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        if (IJ.getInstance() == null) {
            ImagePlus[] inputsCopy = new ImagePlus[]{inputs[0].duplicate(), inputs[1].duplicate(), inputs[2].duplicate()};
            inputs[0] = inputsCopy[gd.getNextChoiceIndex()].duplicate();
            inputs[1] = inputsCopy[gd.getNextChoiceIndex()].duplicate();
            inputs[2] = inputsCopy[gd.getNextChoiceIndex()].duplicate();
        } else {
            inputs[0] = WindowManager.getImage(gd.getNextChoice());
            inputs[1] = WindowManager.getImage(gd.getNextChoice());
            inputs[2] = WindowManager.getImage(gd.getNextChoice());
        }
        useThresh = gd.getNextBoolean();
        threshLevel = gd.getNextNumber();
        isolateFoci = gd.getNextBoolean();
        analyseFluorescence = gd.getNextBoolean();
        histMin = gd.getNextNumber();
        histMax = gd.getNextNumber();
        histNBins = (int) Math.round(gd.getNextNumber());
        return true;
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
    public void drawDetections(ParticleArray pa, int width, int height) {
        int depth = pa.getDepth();
        for (int d = 0; d < depth; d++) {
            ArrayList<Particle> level = pa.getLevel(d);
            ShortProcessor output = new ShortProcessor(width, height);
            for (Particle p : level) {
                Utils.draw2DGaussian(output, p.getC1Gaussian(), 0.0, UserVariables.getSpatialRes(), false);
            }
            output.multiply(1.0 / normFactor);
            IJ.saveAs(new ImagePlus("", output), "TIF", outputDirName + "/Foci Detections " + d);
        }
    }

    /**
     * Duplicates the image containing foci, saves a normalised copy and returns
     * a reference to the saved location
     *
     * @return absolute path to the normalised images
     */
    protected String prepareInputs() {
        if (inputs[1] == null) {
            return null;
        }
        return normaliseStacks(inputs[1].getImageStack(), null);
    }

    /**
     *
     * @param distances
     * @param centroids
     * @param assignments
     * @throws IOException
     */
    void outputFociDistanceData(double[][] distances) throws IOException {
        int N = distances.length;
        TextWindow tw = new TextWindow(title + " Results", resultsHeadings, new String(), 640, 480);
        DecimalFormat df1 = new DecimalFormat("00");
        DecimalFormat df2 = new DecimalFormat("0.000");
        for (int i = 0; i < N; i++) {
            Cell c = cells[i];
            double[] centroid = c.getNucleus().getCentroid();
            ArrayList<Particle> particles = c.getParticles();
            String result = c.getID()+"\t"+df1.format(Math.round(centroid[0])) + "\t" + df1.format(Math.round(centroid[1]));
            if (particles != null) {
                int L = particles.size();
                DescriptiveStatistics intensitites = new DescriptiveStatistics();
                DescriptiveStatistics dist = new DescriptiveStatistics();
                for (int j = 0; j < L; j++) {
                    Particle p = particles.get(j);
                    intensitites.addValue(p.getC1Intensity());
                    dist.addValue(distances[i][j]);
                }
                result = result.concat("\t" + L + "\t" + df2.format(intensitites.getMean()) + "\t" + df2.format(dist.getMean()));
            } else {
                result = result.concat("\t 0\t " + df2.format(0) + "\t " + df2.format(0));
            }
            tw.append(result);
        }
        saveTextWindow(tw, new File(outputDirName + "/foci_distance_data.csv"), resultsHeadings);
    }

    /**
     *
     * @param binaryImage
     */
    ByteProcessor checkBinaryImage(ImageProcessor image) {
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
        double percentile = ds.getPercentile(threshold);
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
            Nucleus nucleus = (Nucleus) cell.getRegion(new Nucleus());
            Cytoplasm cyto = (Cytoplasm) cell.getRegion(new Cytoplasm());
            Roi nucRoi = nucleus.getRoi();
            ByteProcessor nucMask = (ByteProcessor) nucRoi.getMask();
            image.setRoi(nucRoi);
            image.setMask(nucMask);
            ImageStatistics nucstats = ImageStatistics.getStatistics(image, measurements, null);
            Roi cytoRoi = cyto.getRoi();
            ByteProcessor cytoMask = (ByteProcessor) cytoRoi.getMask();
            int xc = nucRoi.getBounds().x - cytoRoi.getBounds().x;
            int yc = nucRoi.getBounds().y - cytoRoi.getBounds().y;
            (new ByteBlitter(cytoMask)).copyBits(nucMask, xc, yc, Blitter.SUBTRACT);
            image.setRoi(cytoRoi);
            image.setMask(cytoMask);
            ImageStatistics cytostats = ImageStatistics.getStatistics(image, measurements, null);
            vals[i] = new double[]{cell.getID(), nucstats.mean, nucstats.stdDev,
                cytostats.mean, cytostats.stdDev,
                nucstats.mean / cytostats.mean,
                nucstats.stdDev / cytostats.stdDev};
        }
        return vals;
    }
}
