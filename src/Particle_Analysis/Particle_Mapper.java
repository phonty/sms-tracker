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
import Math.Histogram;
import Particle.Particle;
import Particle.ParticleArray;
import ParticleTracking.UserVariables;
import Revision.Revision;
import UtilClasses.GenUtils;
import UtilClasses.GenVariables;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
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
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.ShortProcessor;
import ij.text.TextWindow;
import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import ui.DetectionGUI;

public class Particle_Mapper extends Particle_Tracker {

    String outputDirName;
    private ArrayList<Cell> cells;
    /**
     * Title of the application
     */
    protected String title = "Particle Mapper";
    private final String resultsHeadings = "X\t Y\t Number of Foci\t Mean Intensity of Foci\t Mean Foci Distance To Nuclear Boundary (" + IJ.micronSymbol + "m)";

    /**
     * Default constructor.
     */
    public Particle_Mapper() {

    }

    public void run(String arg) {
        Prefs.blackBackground = false;
        File inputDir = null;
        title = title + "_v" + Revision.VERSION + "." + intFormat.format(Revision.revisionNumber);
        stacks = new ImageStack[2];
//        if (IJ.getInstance() != null) {
//            if (!getActiveImages()) {
//                return;
//            }
//        } else {
        inputDir = buildStacks();
        if (inputDir == null) {
            return;
        }
//        }
        ImageProcessor nuclei = IJ.openImage((new OpenDialog("Specify Nuclei Image", null)).getPath()).getProcessor();
        if (!showDialog()) {
            return;
        }
        File outputDir = Utilities.getFolder(inputDir, "Specify directory for output files...", true);
        outputDirName = GenUtils.openResultsDirectory(outputDir + delimiter + title);
        try {
            checkBinaryImage(nuclei);
            IJ.saveAs(new ImagePlus("", nuclei), "PNG", outputDirName + "/Nuclei Mask");
            findCells(nuclei.duplicate());
            ParticleArray pa = findParticles();
            ImageProcessor cellMap = buildTerritories(nuclei.duplicate()).getProcessor();
            IJ.saveAs(new ImagePlus("", cellMap), "PNG", outputDirName + "/Region IDs");
            Collections.sort(cells);
            filterCells(IJ.openImage((new OpenDialog("Specify Region Index Image", null)).getPath()).getProcessor(),
                    new Cytoplasm(), 50.0, Measurements.MEAN);
            assignParticlesToCells(pa, cellMap);
            drawDetections(pa, stacks[0].getWidth(), stacks[0].getHeight());
            double[][] distances = calcDistances(buildDistanceMap(nuclei));
            buildHistograms(distances, 40, 20, -5);
            outputData(distances);
            cleanUp();
        } catch (IOException e) {
            IJ.error(e.getMessage());
        }
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
        cells = new ArrayList();
        ImagePlus imp = new ImagePlus("", image);
        ResultsTable rt = Analyzer.getResultsTable();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE + ParticleAnalyzer.CLEAR_WORKSHEET + ParticleAnalyzer.ADD_TO_MANAGER, Measurements.CENTROID, rt, 0, Double.MAX_VALUE);
        pa.setHideOutputImage(true);
        pa.analyze(imp);
        int n = rt.size();
        RoiManager roimanager = RoiManager.getInstance();
        roimanager.setVisible(false);
        Roi[] rois = roimanager.getRoisAsArray();
        int x = rt.getColumnIndex("X");
        int y = rt.getColumnIndex("Y");
        for (int i = 0; i < n; i++) {
            double[] centroid = new double[]{rt.getValueAsDouble(x, i), rt.getValueAsDouble(y, i)};
            cells.add(new Cell(new Nucleus(rois[i], centroid)));
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
        IJ.saveAs(new ImagePlus("", image), "PNG", outputDirName + "/Cell Boundaries");
        rt.reset();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_ROI_MASKS + ParticleAnalyzer.ADD_TO_MANAGER, 0, rt, 0, Double.MAX_VALUE);
        pa.setHideOutputImage(true);
        pa.analyze(imp);
        int n = rt.size();
        RoiManager roimanager = RoiManager.getInstance();
        roimanager.setVisible(false);
        Roi[] rois = roimanager.getRoisAsArray();
        ImageProcessor cellMap = pa.getOutputImage().getProcessor();
        for (int i = 0; i < n; i++) {
            Cell c = cells.get(i);
            double[] centroid = c.getNucleus().getCentroid();
            int xc = (int) Math.round(centroid[0]);
            int yc = (int) Math.round(centroid[1]);
            int id = cellMap.getPixel(xc, yc) - 1;
            c.setID(id);
            c.addCellRegion(new Cytoplasm(rois[i]));
        }
        return pa.getOutputImage();
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
        for (int i = 0; i < cells.size(); i++) {
            idToIndexMap.put(cells.get(i).getID(), i);
        }
        for (int i = 0; i < N; i++) {
            Particle p = detections.get(i);
            int xp = (int) Math.round(p.getX() / UserVariables.getSpatialRes());
            int yp = (int) Math.round(p.getY() / UserVariables.getSpatialRes());
            Integer pIndex = idToIndexMap.get(cellMap.getPixel(xp, yp) - 1);
            if (pIndex != null) {
                Cell c = cells.get(pIndex);
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
        int N = cells.size();
        double[][] distances = new double[N][];
        double res = UserVariables.getSpatialRes();
        for (int i = 0; i < N; i++) {
            Cell c = cells.get(i);
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
        GenericDialog gd = new GenericDialog("Specify parameters for histogram");
        gd.addNumericField("Minimum Value:", min, 1);
        gd.addNumericField("Maximum Value:", max, 1);
        gd.addNumericField("Number of Bins:", nBins, 0);
        gd.setMinimumSize(new Dimension(400, 100));
        gd.showDialog();
        if (gd.wasOKed()) {
            min = gd.getNextNumber();
            max = gd.getNextNumber();
            nBins = (int) Math.round(nBins);
        }
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
        rt.save(outputDirName + "/distance_histogram_data.csv");
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
        ImagePlus cytoImp = IJ.openImage((new OpenDialog("Specify Foci Image", null)).getPath());
        if (cytoImp == null) {
            return null;
        }

        return normaliseStacks(cytoImp.getImageStack(), null);
    }

    /**
     *
     * @param distances
     * @param centroids
     * @param assignments
     * @throws IOException
     */
    void outputData(double[][] distances) throws IOException {
        int N = distances.length;
        TextWindow tw = new TextWindow(title + " Results", resultsHeadings, new String(), 640, 480);
        DecimalFormat df1 = new DecimalFormat("00");
        DecimalFormat df2 = new DecimalFormat("0.000");
        for (int i = 0; i < N; i++) {
            Cell c = cells.get(i);
            double[] centroid = c.getNucleus().getCentroid();
            ArrayList<Particle> particles = c.getParticles();
            String result = df1.format(Math.round(centroid[0])) + "\t" + df1.format(Math.round(centroid[1]));
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
        saveTextWindow(tw);
    }

    /**
     *
     * @param tw
     * @throws IOException
     */
    void saveTextWindow(TextWindow tw) throws IOException {
        File dataFile = new File(outputDirName + "/cell_data.csv");
        CSVPrinter printer = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(dataFile), GenVariables.ISO), CSVFormat.EXCEL);
        int L = tw.getTextPanel().getLineCount();
        printer.printRecord(resultsHeadings.replace("\t", ","));
        for (int l = 0; l < L; l++) {
            printer.printRecord(tw.getTextPanel().getLine(l).replace("\t", ","));
        }
        printer.close();
    }

    /**
     *
     * @param binaryImage
     */
    void checkBinaryImage(ImageProcessor binaryImage) {
        ImageStatistics stats = binaryImage.getStatistics();
        if (stats.histogram[0] + stats.histogram[255] < binaryImage.getWidth() * binaryImage.getHeight()) {
            binaryImage.autoThreshold();
        }
        if (binaryImage.isInvertedLut()) {
//            binaryImage.invertLut();
        }
        stats = binaryImage.getStatistics();
        if (stats.histogram[0] > stats.histogram[255]) {
            binaryImage.invert();
        }
    }

    void filterCells(ImageProcessor image, CellRegion region, double threshold, int measurement) {
        DescriptiveStatistics ds = new DescriptiveStatistics();
        for (Cell cell : cells) {
            CellRegion cr = cell.getRegionOfType(region);
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
        double percentile = ds.getPercentile(threshold);
        double[] measures = ds.getValues();
        ArrayList<Cell> cells2 = new ArrayList();
        for (int i = 0; i < cells.size(); i++) {
            if (measures[i] > percentile) {
                cells2.add(cells.get(i));
            }
        }
        cells = cells2;
    }

    void analyseCellFluorescenceDistribution(ImageProcessor image) {

    }
}
