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
package ParticleMapping;

import IAClasses.Utils;
import Math.Histogram;
import ParticleTracking.Analyse_;
import ParticleTracking.Particle;
import ParticleTracking.ParticleArray;
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
import ij.io.OpenDialog;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.EDM;
import ij.plugin.filter.ParticleAnalyzer;
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
import java.util.LinkedHashMap;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import ui.DetectionGUI;

/**
 *
 * @author Dave Barry <david.barry at crick.ac.uk>
 */
public class Particle_Mapper extends Analyse_ {

    LinkedHashMap<Integer, double[]> regionCentroidMap = new LinkedHashMap();
    String outputDirName;
    protected String title = "Particle Mapper";
    private final String resultsHeadings = "X\t Y\t Number of Foci\t Mean Intensity of Foci\t Mean Foci Distance To Nuclear Boundary (" + IJ.micronSymbol + "m)";

    public Particle_Mapper() {

    }

    /**
     *
     * @param arg
     */
    public void run(String arg) {
        Prefs.blackBackground = false;
        File inputDir = null;
        title = title + "_v" + Revision.VERSION + "." + intFormat.format(Revision.revisionNumber);
        stacks = new ImageStack[2];
        if (IJ.getInstance() != null) {
            if (!getActiveImages()) {
                return;
            }
        } else {
            inputDir = buildStacks();
            if (inputDir == null) {
                return;
            }
        }
        if (!showDialog()) {
            return;
        }
        File outputDir = Utilities.getFolder(inputDir, "Specify directory for output files...", true);
        outputDirName = GenUtils.openResultsDirectory(outputDir + delimiter + title, delimiter);
        try {
            ParticleArray pa = findParticles();
            drawDetections(pa, stacks[0].getWidth(), stacks[0].getHeight());
            ImageProcessor nuclei = IJ.openImage((new OpenDialog("Specify Nuclei Image", null)).getPath()).getProcessor();
            checkBinaryImage(nuclei);
            IJ.saveAs(new ImagePlus("", nuclei), "PNG", outputDirName + "/Nuclei Mask");
            double[][] centroids = getCentroids(nuclei.duplicate());
            ArrayList<Particle>[] assignments = assignParticlesToCells(
                    pa, buildTerritories(nuclei.duplicate()).getProcessor(), centroids
            );
            double[][] distances = calcDistances(assignments, buildDistanceMap(nuclei));
            buildHistograms(distances, 40, 20, -5);
            outputData(distances, centroids, assignments);
            cleanUp();
        } catch (IOException e) {
            IJ.error("IOException");
        }
    }

    public double[][] getCentroids(ImageProcessor image) {
        ImagePlus imp = new ImagePlus("", image);
        ResultsTable rt = Analyzer.getResultsTable();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE + ParticleAnalyzer.CLEAR_WORKSHEET, Measurements.CENTROID, rt, 0, Double.MAX_VALUE);
        pa.setHideOutputImage(true);
        pa.analyze(imp);
        double[][] centroids = {rt.getColumnAsDoubles(rt.getColumnIndex("X")), rt.getColumnAsDoubles(rt.getColumnIndex("Y"))};
        return centroids;
    }

    public ImagePlus buildTerritories(ImageProcessor image) {
        ImagePlus imp = new ImagePlus("", image);
        ResultsTable rt = Analyzer.getResultsTable();
        EDM edm = new EDM();
        edm.setup("voronoi", imp);
        edm.run(image);
        image.threshold(1);
        IJ.saveAs(new ImagePlus("", image), "PNG", outputDirName + "/Cell Boundaries");
        rt.reset();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_ROI_MASKS, 0, rt, 0, Double.MAX_VALUE);
        pa.setHideOutputImage(true);
        pa.analyze(imp);
        return pa.getOutputImage();
    }

    ArrayList<Particle>[] assignParticlesToCells(ParticleArray pa, ImageProcessor cellMap, double[][] cellCentroids) {
        ByteProcessor map = new ByteProcessor(cellMap.getWidth(), cellMap.getHeight());
        map.setValue(0);
        map.fill();
        map.setValue(255);
        ArrayList<Particle> detections = pa.getLevel(0);
        int M = cellCentroids[0].length;
        for (int i = 0; i < M; i++) {
            int xc = (int) Math.round(cellCentroids[0][i]);
            int yc = (int) Math.round(cellCentroids[1][i]);
            int pIndex = cellMap.getPixel(xc, yc) - 1;
            regionCentroidMap.put(pIndex, new double[]{cellCentroids[0][i], cellCentroids[1][i]});
        }
        int N = detections.size();
        ArrayList<Particle>[] assignments = new ArrayList[M];
        for (int i = 0; i < N; i++) {
            Particle p = detections.get(i);
            int xp = (int) Math.round(p.getX() / UserVariables.getSpatialRes());
            int yp = (int) Math.round(p.getY() / UserVariables.getSpatialRes());
            int pIndex = cellMap.getPixel(xp, yp) - 1;
            if (pIndex >= 0) {
                if (assignments[pIndex] == null) {
                    assignments[pIndex] = new ArrayList();
                }
                assignments[pIndex].add(p);
                double[] centroid = regionCentroidMap.get(pIndex);
                map.drawLine(xp, yp,
                        (int) Math.round(centroid[0]), (int) Math.round(centroid[1]));
            }
        }
        IJ.saveAs(new ImagePlus("", map), "PNG", outputDirName + "/Foci-Nuclei Associations");
        return assignments;
    }

    double[][] calcDistances(ArrayList<Particle>[] assignments, ImageProcessor distanceMap) {
        int N = assignments.length;
        double[][] distances = new double[N][];
        double res = UserVariables.getSpatialRes();
        for (int i = 0; i < N; i++) {
            if (assignments[i] != null) {
                int M = assignments[i].size();
                distances[i] = new double[M];
                for (int j = 0; j < M; j++) {
                    Particle p = assignments[i].get(j);
                    distances[i][j] = res * distanceMap.getInterpolatedValue(p.getX() / res, p.getY() / res);
                }
            }
        }
        return distances;
    }

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

    public boolean showDialog() {
        DetectionGUI ui = new DetectionGUI(null, true, title, this);
        ui.setVisible(true);
        return ui.isWasOKed();
    }

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

    protected String prepareInputs() {
        ImagePlus cytoImp = IJ.openImage((new OpenDialog("Specify Foci Image", null)).getPath());
        if (cytoImp == null) {
            return null;
        }

        return normaliseStacks(cytoImp.getImageStack(), null);
    }

    void outputData(double[][] distances, double[][] centroids, ArrayList<Particle>[] assignments) throws IOException {
        int N = distances.length;
        TextWindow tw = new TextWindow(title + " Results", resultsHeadings, new String(), 640, 480);
        DecimalFormat df1 = new DecimalFormat("00");
        DecimalFormat df2 = new DecimalFormat("0.000");
        for (int i = 0; i < N; i++) {
            String result = df1.format(Math.round(centroids[0][i])) + "\t" + df1.format(Math.round(centroids[1][i]));
            if (assignments[i] != null) {
                int L = assignments[i].size();
                DescriptiveStatistics intensitites = new DescriptiveStatistics();
                DescriptiveStatistics dist = new DescriptiveStatistics();
                for (int j = 0; j < L; j++) {
                    Particle p = assignments[i].get(j);
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

    void checkBinaryImage(ImageProcessor binaryImage) {
        ImageStatistics stats = binaryImage.getStatistics();
        if (stats.histogram[0] + stats.histogram[255] < binaryImage.getWidth() * binaryImage.getHeight()) {
            binaryImage.autoThreshold();
        }
        if (binaryImage.isInvertedLut()) {
            binaryImage.invertLut();
        }
        stats = binaryImage.getStatistics();
        if (stats.histogram[0] > stats.histogram[255]) {
            binaryImage.invert();
        }
    }
}
