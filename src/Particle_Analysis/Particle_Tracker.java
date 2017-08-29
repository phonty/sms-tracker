package Particle_Analysis;

import Adapt.Analyse_Movie;
import Cell.CellData;
import Detection.Blob_Detector;
import ui.ResultsPreviewInterface;
import ui.UserInterface;
import Particle.IsoGaussian;
import IAClasses.ProgressDialog;
import IAClasses.Region;
import IAClasses.Utils;
import MacroWriter.MacroWriter;
import Optimisation.IsoGaussianFitter;
import Particle.Particle;
import Particle.ParticleArray;
import ParticleTracking.ParticleTrajectory;
import ParticleTracking.TrajectoryBuilder;
import ParticleTracking.UserVariables;
import Revision.Revision;
import Segmentation.RegionGrower;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import MetaData.ParamsReader;
import UtilClasses.GenVariables;
import goshtasby.Multi_Goshtasby;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.Straightener;
import ij.plugin.TextReader;
import ij.plugin.filter.GaussianBlur;
import ij.process.AutoThresholder;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.FloatStatistics;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.StackStatistics;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Color;
import java.awt.Rectangle;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import javax.swing.JFileChooser;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

/**
 * Timelapse_Analysis seeks to identify individual particles in a stack of
 * images and plot their trajectories. Particles are first identified in
 * individual colour bands at different time-points by searching for local
 * maxima above a certain threshold, fitting IsoGaussian curves about these
 * maxima, associating detected particles between images and constructing
 * trajectories.
 *
 * @author David J Barry
 * @version 2.0, FEB 2011
 */
public class Particle_Tracker implements PlugIn {

//    protected double[] SIGMAS;
    public final int GOSHTASBY_M = 2, GOSHTASBY_N = 4;
//    public static final int VERSION = 5;
    private final double NORM_VAL = 0.99999;
    protected ArrayList<ParticleTrajectory> trajectories; //Trajectories of the detected particles
    protected long startTime;
    protected DecimalFormat numFormat = new DecimalFormat("0.000");
    protected DecimalFormat intFormat = new DecimalFormat("000");
    protected String title = "Particle Tracker", ext;
    protected final String C_0 = "C0", C_1 = "C1";
    private final double TRACK_WIDTH = 4.0;
    public final float TRACK_EXT = 1.0f;
    public final float TRACK_OFFSET = 0.75f;
    protected static File c0Dir, c1Dir,
            calDir = new File("C:\\Users\\barry05\\Desktop\\2016.03.10_PP1\\Cal");
    protected final String delimiter = GenUtils.getDelimiter();
    protected ImagePlus[] inputs;
    protected final String labels[] = {"Channel 1", "Channel 2"};
    protected static  boolean gpuEnabled = false;
    private final double PART_RAD = 2.0;
    private static File calFile;
    protected double normFactor;

    public Particle_Tracker() {
    }

    public Particle_Tracker(ImageStack[] stacks) {
        this.inputs = new ImagePlus[]{new ImagePlus("C1", stacks[0]), new ImagePlus("C1", stacks[1])};
    }

    public Particle_Tracker(ImagePlus imp, String ext) {
//        this.imp = imp;
//        this.stacks = imp.getImageStack();
        this.ext = ext;
    }

    public void run(String arg) {
        MacroWriter.write();
        File inputDir = null;
        title = title + "_v" + Revision.VERSION + "." + intFormat.format(Revision.revisionNumber);
        inputs = new ImagePlus[2];
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
//        readParamsFromImage();
        try {
            if (showDialog()) {
                analyse(inputDir);
            }
            if (IJ.getInstance() == null) {
                cleanUp();
            }
        } catch (Exception e) {
            GenUtils.error(e.getMessage());
        }

    }

    protected File buildStacks() {
        String dirName = prepareInputs();
        File inputDir;
        if (dirName != null) {
            inputDir = new File(dirName);
        } else {
            return null;
        }
        c0Dir = new File(inputDir.getAbsolutePath() + delimiter + C_0);
        inputs[0] = Utils.buildStack(c0Dir);
        this.ext = inputs[0].getTitle();
        c1Dir = new File(inputDir.getAbsolutePath() + delimiter + C_1);
        if (c1Dir.exists()) {
            inputs[1] = Utils.buildStack(c1Dir);
        }
        return inputDir;
    }

    protected void cleanUp() {
        try {
            FileUtils.deleteDirectory(c0Dir);
            FileUtils.deleteDirectory(c1Dir);
        } catch (Exception e) {
            IJ.error(e.toString());
        }
    }

    protected String prepareInputs() {
        ImagePlus cytoImp = IJ.openImage();
        if (cytoImp == null) {
            return null;
        }
        ImagePlus sigImp = IJ.openImage();
        ImageStack cytoStack = cytoImp.getImageStack();
        int cytoSize = cytoStack.getSize();
        ImageStack sigStack = null;
        if (sigImp != null) {
            sigStack = sigImp.getImageStack();
            if (sigStack.getSize() != cytoSize) {
                IJ.error("Stacks must contain same number of slices.");
                return null;
            }
        }

        return normaliseStacks(cytoStack, sigStack);
    }

    protected String normaliseStacks(ImageStack cytoStack, ImageStack sigStack) {
        ImagePlus cytoImp = new ImagePlus("", cytoStack);
        int cytoSize = cytoStack.getSize();
        StackStatistics cytoStats = new StackStatistics(cytoImp);
        double stackMin = cytoStats.min;
        for (int i = 1; i <= cytoSize; i++) {
            cytoStack.getProcessor(i).subtract(stackMin);
        }
        StackStatistics cytoStats2 = new StackStatistics(cytoImp);
        int histogram[] = cytoStats2.histogram16;
        if (histogram == null) {
            histogram = cytoStats2.histogram;
        }
        double sum = 0.0;
        int nPixels = cytoStats2.pixelCount;
        int thresh = 0;
        while (sum < nPixels * NORM_VAL && thresh < histogram.length) {
            sum += histogram[thresh++];
        }
        normFactor = 100.0 / thresh;
        GenUtils.convertStack(cytoImp, 32);
        cytoStack = cytoImp.getImageStack();
        for (int i = 1; i <= cytoSize; i++) {
            cytoStack.getProcessor(i).multiply(normFactor);
        }
        String parent = IJ.getDirectory("current");
        File seriesFolder = GenUtils.createDirectory(parent + FilenameUtils.getBaseName(cytoImp.getTitle()), false);
        File c0Folder = GenUtils.createDirectory(seriesFolder.getAbsolutePath() + delimiter + C_0, true);
        Utils.saveStackAsSeries(cytoStack, c0Folder.getAbsolutePath() + delimiter,
                "TIF", intFormat);
        if (sigStack != null) {
            File c1Folder = GenUtils.createDirectory(seriesFolder.getAbsolutePath() + delimiter + C_1, true);
            Utils.saveStackAsSeries(sigStack, c1Folder.getAbsolutePath() + delimiter,
                    "TIF", intFormat);
        }
        return seriesFolder.getAbsolutePath();
    }

    public boolean showDialog() {
        UserInterface ui = new UserInterface(null, true, title, this);
        ui.setVisible(true);
        return ui.isWasOKed();
    }

    /**
     * Analyses the {@link ImageStack} specified by <code>stack</code>.
     */
    public void analyse(File inputDir) throws Exception {
        ImageStack[] stacks = getStacks();
        File outputDir = Utilities.getFolder(inputDir, "Specify directory for output files...", true);
        String parentDir = GenUtils.openResultsDirectory(outputDir + delimiter + title);
        String sigc0Dir = GenUtils.openResultsDirectory(parentDir + delimiter + "C0");
        String sigc1Dir = GenUtils.openResultsDirectory(parentDir + delimiter + "C1");
        ParticleTrajectory.resetMSDPlot();
        if (!(stacks[1] == null) && UserVariables.isUseCals()) {
            JFileChooser fileChooser = new JFileChooser(calDir);
            fileChooser.setDialogTitle("Specify file containing bead calibration data");
            fileChooser.showOpenDialog(null);
            calFile = fileChooser.getSelectedFile();
            calDir = calFile.getParentFile();
            if (!(new Multi_Goshtasby()).run(calFile, Bead_Calibration.HEADER_SIZE)) {
                UserVariables.setUseCals(false);
                System.out.println("Calculation of calibration coefficients failed - proceeding without calibration.");
            }
        }
        if (stacks != null) {
            IJ.register(this.getClass());
            startTime = System.currentTimeMillis();
            trajectories = new ArrayList();
            findParticles();
            TextWindow results = new TextWindow(title + " Results", "X\tY\tFrame\tChannel 1\tChannel 2\tChannel 1 " + '\u03C3'
                    + "\tChannel 2 " + '\u03C3',
                    new String(), 1000, 500);
            TextWindow resultSummary = new TextWindow(title + " Results Summary",
                    "Particle\tDuration (s)\tDisplacement (" + IJ.micronSymbol
                    + "m)\tVelocity (" + IJ.micronSymbol + "m/s)\tDirectionality\tDiffusion Coefficient ("
                    + IJ.micronSymbol + "m^2/s)",
                    new String(), 1200, 500);
            int n = trajectories.size();
            for (int i = 0; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (!(traj.getDisplacement(traj.getEnd(), traj.getSize()) > UserVariables.getMinTrajDist()
                        && traj.getDuration() > UserVariables.getMinTrajLength()
                        && checkTrajColocalisation(traj, UserVariables.getColocalThresh(), UserVariables.isColocal()))) {
                    trajectories.remove(i);
                    i--;
                    n--;
                }
            }
            if (UserVariables.isPrevRes()) {
                ArrayList<Integer> includeList = previewResults();
                if (includeList != null) {
                    ArrayList<ParticleTrajectory> temps = new ArrayList();
                    for (Integer e : includeList) {
                        temps.add(trajectories.get(e));
                    }
                    trajectories = new ArrayList();
                    trajectories.addAll(temps);
                }
            }
            n = trajectories.size();
            ProgressDialog trajProg = new ProgressDialog(null, "Analysing trajectories...", false, title, false);
            trajProg.setVisible(true);
            for (int i = 0; i < n; i++) {
                trajProg.updateProgress(i, n);
                boolean remove = false;
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (traj != null) {
                    traj.calcMSD(-1, i + 1);
                    printData(i, resultSummary, i + 1);
                    traj.printTrajectory(i + 1, results, numFormat, title);
                    if (stacks[1] != null && UserVariables.isExtractsigs()) {
                        ImageStack signals[] = extractTrajSignalValues(traj,
                                (int) Math.round(UserVariables.getTrackLength() / UserVariables.getSpatialRes()),
                                (int) Math.round(TRACK_WIDTH / UserVariables.getSpatialRes()),
                                TRACK_EXT / ((float) UserVariables.getSpatialRes()), stacks[0].getWidth(), stacks[0].getHeight(), i + 1);
                        if (signals[0].getSize() > 0) {
                            for (int j = 1; j <= signals[0].getSize(); j++) {
                                IJ.saveAs((new ImagePlus("", signals[0].getProcessor(j))),
                                        "TIF", sigc0Dir + delimiter + "C0-" + String.valueOf(i + 1)
                                        + "-" + signals[0].getSliceLabel(j));
                                IJ.saveAs((new ImagePlus("", signals[1].getProcessor(j))),
                                        "TIF", sigc1Dir + delimiter + "C1-" + String.valueOf(i + 1)
                                        + "-" + signals[1].getSliceLabel(j));
                            }
                        } else {
                            remove = true;
                        }
                    }
                } else {
                    remove = true;
                }
                if (remove) {
                    trajectories.remove(i);
                    i--;
                    n--;
                }
            }
            ParticleTrajectory.drawGlobalMSDPlot();
            trajProg.dispose();
            if (trajectories.size() > 0) {
                ImageStack maps = mapTrajectories((new RGBStackMerge()).mergeStacks(stacks[0].getWidth(), stacks[0].getHeight(), stacks[0].getSize(), stacks[0], stacks[1], null, true),
                        trajectories, UserVariables.getSpatialRes(), UserVariables.getMinTrajLength(),
                        UserVariables.getTimeRes(), true, 0, trajectories.size() - 1, 1, false, calcParticleRadius(UserVariables.getSpatialRes(), UserVariables.getSigEstRed()));
                resultSummary.append("\nAnalysis Time (s): " + numFormat.format((System.currentTimeMillis() - startTime) / 1000.0));
                results.setVisible(true);
                resultSummary.setVisible(true);
                IJ.saveString(results.getTextPanel().getText(), parentDir + "/results.txt");
                IJ.saveString(resultSummary.getTextPanel().getText(), parentDir + "/resultsSummary.txt");
                Plot msdPlot = ParticleTrajectory.getMsdPlot();
                try {
                    printTrajectories(trajectories, new File(String.format("%s%s%s", parentDir, File.separator, "AllParticleData.csv")), stacks[0].size());
                    if (msdPlot != null) {
                        msdPlot.getResultsTable().saveAs(parentDir + "/MSD_Plot.csv");
                        IJ.saveAs(msdPlot.makeHighResolution("", 10.0f, true, false), "PNG", parentDir + "/MSD_Plot");
                    }
                } catch (IOException e) {
                }
                if (maps != null) {
                    (new ImagePlus("Trajectory Maps", maps)).show();
                    IJ.saveAs((new ImagePlus("", maps)), "TIF", parentDir + "/trajectories.tif");
                }
            } else {
                IJ.log("No Particle Trajectories Constructed.");
            }
        }
        printParams(parentDir);
    }

    protected ParticleArray findParticles() {
        ImageStack[] stacks = getStacks();
        return findParticles(true, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1], UserVariables.isBlobs());
    }

    /**
     * Median filter and IsoGaussian filter the image specified by
     * <code>processor</code>.
     *
     * @param processor the image to be pre-processed.
     */
    protected ImageProcessor preProcess(ImageProcessor processor, double sigma) {
        if (processor == null) {
            return null;
        }
        ImageProcessor fp = (new TypeConverter(processor, false)).convertToFloat(null);
        fp.noise(0.0001);
        if (UserVariables.isPreProcess()) {
            (new GaussianBlur()).blurGaussian(fp, sigma, sigma, 0.1);
        }
        return fp;
    }

    public ParticleArray findParticles(boolean update, int startSlice, int endSlice,
            double c1FitTol, ImageStack channel1, ImageStack channel2, boolean blobs) {
        return findParticles(update, startSlice, endSlice, c1FitTol,
                channel1, channel2, true, false, UserVariables.isFitC2(), blobs);
    }

    public ParticleArray findParticles(boolean update, int startSlice, int endSlice,
            double c1FitTol, ImageStack channel1, ImageStack channel2, boolean showProgress,
            boolean floatingSigma, boolean fitC2, boolean blobs) {
        if (channel1 == null) {
            return null;
        }
        int i, noOfImages = channel1.getSize(), width = channel1.getWidth(), height = channel1.getHeight(),
                arraySize = endSlice - startSlice + 1;
        int searchRad = calcParticleRadius(UserVariables.getSpatialRes(), UserVariables.getSigEstRed());
        int fitRad = searchRad;
        int pSize = 2 * fitRad + 1;
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
        ProgressDialog progress = new ProgressDialog(null, "Finding Particles...", false, title, false);
        if (showProgress) {
            progress.setVisible(true);
        }
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
            IJ.freeMemory();
            progress.updateProgress(i - startSlice, arraySize);
            if (UserVariables.isBlobs()) {
                detectBlobs(startSlice, channel1, channel2, i, pSize,
                        width, height, searchRad, fitRad, particles, fitC2);
            } else {
                detectParticles(startSlice, channel1, channel2, i,
                        width, height, searchRad, fitRad, particles,
                        floatingSigma, spatialRes, fitC2, c1FitTol);
            }
        }
        progress.dispose();
        if (update) {
            updateTrajs(particles, spatialRes, update);
        }
        return particles;
    }

    public void detectBlobs(int startSlice, ImageStack channel1, ImageStack channel2,
            int i, int pSize, int width, int height, int searchRad, int fitRad,
            ParticleArray particles, boolean fitC2) {
        ImageProcessor chan1Proc = channel1.getProcessor(i + 1).convertToFloat();
        ImageProcessor chan2Proc = (channel2 != null) ? channel2.getProcessor(i + 1).convertToFloat() : null;
        ByteProcessor C2Max = null, C1Max = getMaxima(pSize, chan1Proc.duplicate(), searchRad, fitRad, UserVariables.getChan1MaxThresh());
        if (chan2Proc != null && fitC2) {
            C2Max = getMaxima(pSize, chan2Proc.duplicate(), searchRad, fitRad, UserVariables.getChan2MaxThresh());
        }
        for (int c1X = 0; c1X < width; c1X++) {
            for (int c1Y = 0; c1Y < height; c1Y++) {
                if (C1Max.getPixel(c1X, c1Y) == UserVariables.FOREGROUND) {
                    double px = c1X * UserVariables.getSpatialRes();
                    double py = c1Y * UserVariables.getSpatialRes();
                    Particle p1 = new Particle(i - startSlice, px, py, chan1Proc.getPixelValue(c1X, c1Y));
                    Particle p2 = chan2Proc != null ? new Particle(i - startSlice, px, py, chan2Proc.getPixelValue(c1X, c1Y)) : null;
                    if (chan2Proc != null && fitC2) {
                        int[][] c2Points = Utils.searchNeighbourhood(c1X, c1Y, searchRad, UserVariables.FOREGROUND, C2Max);
                        if (c2Points != null) {
                            px = c2Points[0][0] * UserVariables.getSpatialRes();
                            py = c2Points[0][1] * UserVariables.getSpatialRes();
                            p2 = new Particle(i - startSlice, px, py, chan2Proc.getPixelValue(c2Points[0][0], c2Points[0][1]));
                        } else {
                            p2 = null;
                        }
                    }
                    p1.setColocalisedParticle(p2);
                    particles.addDetection(i - startSlice, p1);
                }
            }
        }
        if (UserVariables.isTrackRegions()) {
            findRegions(particles.getLevel(i - startSlice), i - startSlice, chan1Proc);
        }
    }

    void findRegions(ArrayList<Particle> particles, int frame, ImageProcessor inputProc) {
        Analyse_Movie am = new Analyse_Movie();
        PointRoi points = null;
        int N = particles.size();
        for (Particle p : particles) {
            if (points == null) {
                points = new PointRoi(p.getX(), p.getY());
            } else {
                points.addPoint(p.getX(), p.getY());
            }
        }
        ArrayList<CellData> cd = new ArrayList();
        RegionGrower.initialiseROIs(null, -1, frame + 1, null, points, inputProc.getWidth(), inputProc.getHeight(), frame + 1, cd, null, false);
        int t = RegionGrower.getThreshold(inputProc, true, -1, AutoThresholder.Method.Default.toString());
        ArrayList<Region> regions = RegionGrower.findCellRegions(inputProc, t, am.getCellData());
        for (int i = 0; i < cd.size(); i++) {
            Region[] allRegions = new Region[inputs[0].getImageStackSize()];
            allRegions[frame] = regions.get(i);
            cd.get(i).setCellRegions(allRegions);
            cd.get(i).setGreyThresholds(new int[]{t});
        }
        try {
            am.getMorphologyData(cd, false, Measurements.MEAN + Measurements.SHAPE_DESCRIPTORS, inputProc.duplicate(), UserVariables.getSigEstRed());
        } catch (Exception e) {
        }
        for (int i = 0; i < N; i++) {
            particles.get(i).setRegion(regions.get(i));
        }
    }

    ByteProcessor getMaxima(int pSize, ImageProcessor ip, int searchRad, int fitRad, double pThreshold) {
        Blob_Detector bd1 = new Blob_Detector(UserVariables.getSigEstRed(), pSize);
        ImageProcessor log1 = bd1.laplacianOfGaussian(ip);
        log1.invert();
        double t = Utils.getPercentileThresh(log1, pThreshold);
        return Utils.findLocalMaxima(searchRad, searchRad, UserVariables.FOREGROUND, log1, t, false, fitRad - searchRad);
    }

    public void detectParticles(int startSlice, ImageStack channel1, ImageStack channel2, int i,
            int width, int height, int searchRad, int fitRad, ParticleArray particles,
            boolean floatingSigma, double spatialRes, boolean fitC2Gaussian, double c1FitTol) {
        int pSize = 2 * fitRad + 1;
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        FloatProcessor chan1Proc = (FloatProcessor) preProcess(channel1.getProcessor(i + 1).duplicate(), UserVariables.getSigEstRed());
        FloatProcessor chan2Proc = (channel2 != null) ? (FloatProcessor) preProcess(channel2.getProcessor(i + 1).duplicate(), UserVariables.getSigEstGreen()) : null;
        double c1Threshold = Utils.getPercentileThresh(chan1Proc, UserVariables.getChan1MaxThresh());
        double c2Threshold = chan2Proc == null ? 0.0 : Utils.getPercentileThresh(chan2Proc, UserVariables.getChan2MaxThresh());
        ByteProcessor thisC1Max = Utils.findLocalMaxima(searchRad, searchRad, UserVariables.FOREGROUND, chan1Proc, c1Threshold, false, fitRad - searchRad);
        for (int c1X = 0; c1X < width; c1X++) {
            for (int c1Y = 0; c1Y < height; c1Y++) {
                if (thisC1Max.getPixel(c1X, c1Y) == UserVariables.FOREGROUND) {
                    Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
                    ArrayList<IsoGaussian> c1Fits = doFitting(xCoords, yCoords, pixValues,
                            floatingSigma, c1X, c1Y, fitRad, spatialRes, c1Threshold, i - startSlice,
                            UserVariables.getSigEstRed() / UserVariables.getSpatialRes());
                    Particle p2 = null;
                    if (chan2Proc != null && chan2Proc.getPixelValue(c1X, c1Y) > c2Threshold) {
                        p2 = new Particle(i - startSlice, c1X * UserVariables.getSpatialRes(), c1Y * UserVariables.getSpatialRes(), chan2Proc.getPixelValue(c1X, c1Y));
                        if (fitC2Gaussian) {
                            Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan2Proc);
                            ArrayList<IsoGaussian> c2Fits = doFitting(xCoords, yCoords, pixValues,
                                    floatingSigma, c1X, c1Y, fitRad, spatialRes,
                                    Utils.getPercentileThresh(chan2Proc, UserVariables.getChan2MaxThresh()),
                                    i - startSlice, UserVariables.getSigEstGreen() / UserVariables.getSpatialRes());
                            p2 = c2Fits.get(0);
                            if (((IsoGaussian) p2).getFit() < c1FitTol) {
                                p2 = null;
                            }
                        }
                    }
                    if (c1Fits != null) {
                        for (IsoGaussian c1Fit : c1Fits) {
                            if (c1Fit.getFit() > c1FitTol) {
                                c1Fit.setColocalisedParticle(p2);
                                particles.addDetection(i - startSlice, c1Fit);
                            }
                        }
                    }
                }
            }
        }
    }

    protected ArrayList<IsoGaussian> doFitting(double[] xCoords, double[] yCoords,
            double[][] pixValues, boolean floatingSigma, int x0, int y0, int fitRad,
            double spatialRes, double c1Threshold, int t, double sigma) {
        IsoGaussianFitter fitter = new IsoGaussianFitter(xCoords, yCoords, pixValues, floatingSigma, sigma);
        fitter.doFit();
        ArrayList<IsoGaussian> fits = new ArrayList();
        double x = (x0 - fitRad + fitter.getX0()) * spatialRes;
        double y = (y0 - fitRad + fitter.getY0()) * spatialRes;
        fits.add(new IsoGaussian(t, x, y, fitter.getMag(), fitter.getXsig(),
                fitter.getYsig(), fitter.getRSquared(), null, -1, null));
        return fits;
    }

    protected void updateTrajs(ParticleArray particles, double spatialRes, boolean update) {
        if (update) {
            TrajectoryBuilder.updateTrajectories(particles, UserVariables.getTimeRes(), UserVariables.getTrajMaxStep(), spatialRes, true, Utils.getStackMinMax(inputs[0].getImageStack())[1], trajectories, UserVariables.isTrackRegions());
        }
    }

    /**
     * Produces a {@link Plot} of the trajectory specified by
     * <code>particleNumber</code>.
     *
     * @param width the width of the images from which the trajectory was
     * extracted
     * @param height the height of the images from which the trajectory was
     * extracted
     * @param particleNumber the trajectory index
     */
    public boolean plotTrajectory(int width, int height, int particleNumber, int label) {
        if (width * height == 0 || trajectories.size() < 1) {
            return false;
        }
        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
        if (traj == null) {
            return false;
        }
        Particle current = traj.getEnd();
        int size = traj.getSize();
        double xValues[] = new double[size];
        double yValues[] = new double[size];
        width *= UserVariables.getSpatialRes();
        height *= UserVariables.getSpatialRes();

        for (int i = size - 1; i >= 0; i--) {
            xValues[i] = (double) current.getX() / width;
            /*
             * Y-coordinates are inverted so trajectory is displayed as per the
             * image
             */
            yValues[i] = -(double) current.getY() / height;
            current = current.getLink();
        }

        Plot plot = new Plot("Particle " + label + " Trajectory",
                "X", "Y", xValues, yValues,
                (Plot.X_TICKS + Plot.Y_TICKS + Plot.X_NUMBERS + Plot.Y_NUMBERS));
        plot.setLimits(0, 1.0, -1.0, 0);
        plot.setLineWidth(2);
        plot.setColor(Color.BLUE);
        plot.draw();
        plot.show();

        return true;
    }

    /**
     * Outputs velocity and directionality data on the particle specified by
     * <code>particleNumber</code>. Directionality ( <code>D</code>) is
     * calculated according to: <br> <br>
     * <code>D = 1 / (1 + &lambda;<sub>1</sub>
     * &lambda;<sub>2</sub><sup>-1</sup>)</code>
     * <br> <br> where <code>&lambda;<sub>1</sub></code> and
     * <code>&lambda;<sub>2</sub></code> are the eigenvalues of the trajectory
     * data and      <code>&lambda;<sub>1</sub>&lambda;<sub>2</sub></code>.
     *
     */
    public boolean printData(int particleNumber, TextWindow output, int label) {
        if (trajectories.size() < 1) {
            return false;
        }
        DecimalFormat decFormat = new DecimalFormat("0.000");
        DecimalFormat msdFormat = new DecimalFormat("0.000000");
        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
        if (traj == null) {
            return false;
        }
        double points[][] = traj.getPoints();
        traj.calcDirectionality(points[0], points[1]);
        double displacement = traj.getDisplacement(traj.getEnd(), traj.getSize());
        double duration = traj.getDuration();
        output.append(label + "\t"
                + decFormat.format(duration) + "\t"
                + decFormat.format(displacement)
                + "\t" + decFormat.format(displacement / duration) + "\t"
                + decFormat.format(traj.getDirectionality()) + "\t"
                + msdFormat.format(traj.getDiffCoeff()));
        return true;
    }

    /**
     * Constructed trajectories are drawn onto the original image sequence and
     * displayed as a stack sequence.
     */
    public ImageStack mapTrajectories(ImageStack stack, ArrayList<ParticleTrajectory> trajectories, double spatialRes, double minTrajLength, double timeRes, boolean tracks, int startT, int endT, int index, boolean preview, int radius) {
        if (stack == null) {
            return null;
        }
        int i, j, width = stack.getWidth(), height = stack.getHeight(),
                frames = stack.getSize();
        double lastX, lastY;
        ImageStack outputStack = new ImageStack(width, height);
        Particle current;
        ParticleTrajectory traj;
        int n = trajectories.size();
        ImageProcessor processor;
        int lastTP;
        if (n < 1) {
            return null;
        }
        for (i = 0; i < frames; i++) {
            processor = (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToRGB());
            processor.setInterpolationMethod(ImageProcessor.BICUBIC);
            processor.setInterpolate(true);
            outputStack.addSlice("" + i, processor.resize(width, height));
        }
        Random r = new Random();
        int tLength = (int) Math.round(UserVariables.getTrackLength() / UserVariables.getSpatialRes());
        ProgressDialog progress = new ProgressDialog(null, "Mapping Output...", false, title, false);
        progress.setVisible(true);
        for (i = startT; i <= endT && i < n; i++) {
            progress.updateProgress(i, n);
            traj = (ParticleTrajectory) (trajectories.get(i));
            if (traj != null) {
                Color thiscolor = new Color(r.nextInt(256), r.nextInt(256), r.nextInt(256));
                current = traj.getEnd();
                lastX = current.getX();
                lastY = current.getY();
                lastTP = current.getTimePoint();
                current = current.getLink();
                while (current != null) {
                    for (j = frames - 1; j >= lastTP; j--) {
                        processor = outputStack.getProcessor(j + 1);
                        processor.setColor(thiscolor);
                        if (j - 1 < lastTP) {
                            markParticle(processor, (int) Math.round(lastX / spatialRes) - radius,
                                    (int) Math.round(lastY / spatialRes) - radius, radius, true, "" + index);
                        }
                        if (tracks && j <= lastTP + tLength) {
                            int x = (int) (Math.round(current.getX() / spatialRes));
                            int y = (int) (Math.round(current.getY() / spatialRes));
                            processor.drawLine(x, y, (int) Math.round(lastX / spatialRes),
                                    (int) Math.round(lastY / spatialRes));
                        }
                    }
                    lastX = current.getX();
                    lastY = current.getY();
                    lastTP = current.getTimePoint();
                    current = current.getLink();
                }
                processor = outputStack.getProcessor(lastTP + 1);
                processor.setColor(thiscolor);
                markParticle(processor, (int) Math.round(lastX / spatialRes) - radius,
                        (int) Math.round(lastY / spatialRes) - radius, radius, true, "" + index);
                index++;
            }
        }
        progress.dispose();
        return outputStack;
    }

    void markParticle(ImageProcessor processor, int x, int y, int radius, boolean string, String label) {
        processor.drawRect(x, y, 2 * radius + 1, 2 * radius + 1);
        if (string) {
            processor.drawString(label, x, y);
        }
    }

    public int calcParticleRadius(double spatialRes) {
        return calcParticleRadius(spatialRes, UserVariables.getSigEstRed());
    }

    public int calcParticleRadius(double spatialRes, double sigEst) {
        return (int) Math.ceil(PART_RAD * sigEst / spatialRes);
    }

    public ArrayList<Integer> previewResults() {
        if (trajectories.size() < 1) {
            return null;
        }
        ResultsPreviewInterface previewUI = new ResultsPreviewInterface(IJ.getInstance(), true, title, this);
        previewUI.setVisible(true);
        if (previewUI.isWasOKed()) {
            return previewUI.getRemoveList();
        } else {
            return null;
        }
    }

    public ArrayList<ParticleTrajectory> getTrajectories() {
        return trajectories;
    }

    public ImageStack[] getStacks() {
        ImageStack[] stacks = new ImageStack[inputs.length];
        for (int i = 0; i < inputs.length; i++) {
            stacks[i] = inputs[i] != null ? inputs[i].getImageStack() : null;
        }
        return stacks;
    }

    public boolean isMonoChrome() {
        return inputs[1] == null;
    }

    /**
     * Procedure for obtaining Goshtasby coefficients - see Goshtasby (1988),
     * IEEE Transactions on Geoscience and Remote Sensing, 26:60-64
     *
     * 1. "Track" beads in two channels to provide list of coordinates 2. Input
     * list of coordinates for each channel to MATLAB 3. Run goshtasby.m to
     * obtain translation coefficients for both x and y, mapping 'green'
     * coordinates onto 'red'. 4. Export x-coeffs, y-coeffs and original green
     * channel coords (assuming virus is in red).
     *
     * @param ptraj
     * @param signalLength
     * @param signalWidth
     * @return
     */
    ImageStack[] extractTrajSignalValues(ParticleTrajectory ptraj, int signalLength, int signalWidth, float offset, int width, int height, int count) {
        ImageStack[] stacks = getStacks();
        String calParent = calDir + delimiter + "goshtasby" + delimiter + GOSHTASBY_M + "_" + GOSHTASBY_N;
        TextReader reader = new TextReader();
        double xdiv = width * UserVariables.getSpatialRes() / GOSHTASBY_M;
        double ydiv = height * UserVariables.getSpatialRes() / GOSHTASBY_N;
        Particle sigStartP = ptraj.getEnd();
        if (signalWidth % 2 == 0) {
            signalWidth++;
        }
        int size = ptraj.getSize();
        float xSigArray[], ySigArray[], xVirArray[], yVirArray[];
        ArrayList<Float> xSigPoints, ySigPoints, xVirPoints, yVirPoints;
        ImageProcessor[] sigTemps = new ImageProcessor[size];
        ImageProcessor[] virTemps = new ImageProcessor[size];
        Particle last = null, next;
//        ProgressDialog progress = new ProgressDialog(null, "Extracting Signal Areas " + count + "...", false, title, false);
//        progress.setVisible(true);
        double radius = calcParticleRadius(UserVariables.getSpatialRes(),
                UserVariables.getSigEstRed()) * UserVariables.getSpatialRes();
        for (int i = 0; i < size; i++) {
//            progress.updateProgress(i, size);
            Particle current = sigStartP;
            next = sigStartP.getLink();
            xSigPoints = new ArrayList();
            ySigPoints = new ArrayList();
            xVirPoints = new ArrayList();
            yVirPoints = new ArrayList();
            double x1, y1, x2, y2, t1, t2;
            if (last != null) {
                x1 = last.getX();
                y1 = last.getY();
                t1 = last.getTimePoint() / UserVariables.getTimeRes();
            } else {
                x1 = sigStartP.getX();
                y1 = sigStartP.getY();
                t1 = sigStartP.getTimePoint() / UserVariables.getTimeRes();
            }
            if (next != null) {
                x2 = next.getX();
                y2 = next.getY();
                t2 = next.getTimePoint() / UserVariables.getTimeRes();
            } else {
                x2 = sigStartP.getX();
                y2 = sigStartP.getY();
                t2 = sigStartP.getTimePoint() / UserVariables.getTimeRes();
            }
            double vel;
            if (t1 - t2 > 0.0) {
                vel = Utils.calcDistance(x1, y1, x2, y2) / (t1 - t2);
            } else {
                vel = 0.0;
            }
            for (int index = 1; index <= size && current != null; index++) {
                double xg = current.getX();
                double yg = current.getY();
                if (UserVariables.isUseCals()) {
                    double x = current.getX(), y = current.getY();
                    int xi = 1 + (int) Math.floor(x / xdiv);
                    int yi = 1 + (int) Math.floor(y / ydiv);
                    if (new File(calParent + delimiter + "xcoeffs" + xi + "_" + yi + ".txt").exists()) {
                        ImageProcessor xcoeffs = reader.open(calParent + delimiter + "xcoeffs" + xi + "_" + yi + ".txt");
                        ImageProcessor ycoeffs = reader.open(calParent + delimiter + "ycoeffs" + xi + "_" + yi + ".txt");
                        ImageProcessor coords = reader.open(calParent + delimiter + "coords" + xi + "_" + yi + ".txt");
                        if (xcoeffs == null || ycoeffs == null || coords == null) {
                            GenUtils.error("Incomplete calibration parameters - proceeding without.");
                            UserVariables.setUseCals(false);
                        } else {
                            xg = goshtasbyEval(xcoeffs, coords, x, y);
                            yg = goshtasbyEval(ycoeffs, coords, x, y);
                        }
                    }
                }
                xSigPoints.add((float) (xg / UserVariables.getSpatialRes()));
                ySigPoints.add((float) (yg / UserVariables.getSpatialRes()));
                xVirPoints.add((float) (current.getX() / UserVariables.getSpatialRes()));
                yVirPoints.add((float) (current.getY() / UserVariables.getSpatialRes()));
                current = current.getLink();
            }
            xSigArray = new float[xSigPoints.size() + 2];
            ySigArray = new float[ySigPoints.size() + 2];
            xVirArray = new float[xVirPoints.size() + 2];
            yVirArray = new float[yVirPoints.size() + 2];
            for (int j = 1; j < xSigPoints.size(); j++) {
                xSigArray[j] = xSigPoints.get(j - 1);
                ySigArray[j] = ySigPoints.get(j - 1);
                xVirArray[j] = xVirPoints.get(j - 1);
                yVirArray[j] = yVirPoints.get(j - 1);
            }
            extendSignalArea(xSigArray, ySigArray, offset, 1);
            extendSignalArea(xVirArray, yVirArray, offset, 1);
            PolygonRoi sigProi = new PolygonRoi(xSigArray, ySigArray, xSigArray.length, Roi.POLYLINE);
            PolygonRoi virProi = new PolygonRoi(xVirArray, yVirArray, xVirArray.length, Roi.POLYLINE);
            Straightener straightener = new Straightener();
            ImagePlus sigImp = new ImagePlus("", stacks[1].getProcessor(sigStartP.getTimePoint() + 1));
            ImagePlus virImp = new ImagePlus("", stacks[0].getProcessor(sigStartP.getTimePoint() + 1));
            sigImp.setRoi(sigProi);
            virImp.setRoi(virProi);
            sigTemps[size - 1 - i] = straightener.straighten(sigImp, sigProi, signalWidth);
            virTemps[size - 1 - i] = straightener.straighten(virImp, virProi, signalWidth);
            if (virTemps[size - 1 - i] != null) {
                virTemps[size - 1 - i].putPixelValue(0, 0, sigStartP.getTimePoint());
                virTemps[size - 1 - i].putPixelValue(1, 0, vel);
            }
            last = sigStartP;
            sigStartP = sigStartP.getLink();
        }
//        progress.dispose();
        int xc = (int) Math.ceil(TRACK_OFFSET * offset);
        int yc = (signalWidth - 1) / 2;
        int outputWidth = (int) Math.round(signalLength + offset);
        ImageStack output[] = new ImageStack[2];
        output[0] = new ImageStack(outputWidth, signalWidth);
        output[1] = new ImageStack(outputWidth, signalWidth);

//        progress = new ProgressDialog(null, "Aligning Signal Areas " + count + "...", false, title, false);
//        progress.setVisible(true);
        for (int j = 0; j < size; j++) {
//            progress.updateProgress(j, size);
            if (virTemps[j] != null && sigTemps[j] != null && sigTemps[j].getWidth() >= outputWidth) {
                Particle alignmentParticle = null;
                if (UserVariables.isUseCals()) {
                    ImageStack virStack = new ImageStack(virTemps[j].getWidth(), virTemps[j].getHeight());
                    virStack.addSlice(virTemps[j]);
                    ParticleArray particles = findParticles(false, 0, 0, 0.0, virStack, null, false, false, false, UserVariables.isBlobs());
                    ArrayList<Particle> detections = particles.getLevel(0);
                    double mindist = Double.MAX_VALUE;
                    int minindex = -1;
                    for (int k = 0; k < detections.size(); k++) {
                        Particle p = detections.get(k);
                        double dist = Utils.calcDistance(p.getX(), p.getY(), xc * UserVariables.getSpatialRes(), yc * UserVariables.getSpatialRes());
                        if (dist < mindist && dist < radius) {
                            mindist = dist;
                            minindex = k;
                        }
                    }
                    if (minindex > -1) {
                        alignmentParticle = detections.get(minindex);
                    }
                }
                if (!UserVariables.isUseCals() || alignmentParticle != null) {
                    String label = Float.toString(virTemps[j].getPixelValue(0, 0))
                            + "-" + numFormat.format(virTemps[j].getPixelValue(1, 0));
                    virTemps[j].setInterpolate(true);
                    virTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    sigTemps[j].setInterpolate(true);
                    sigTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    double xinc = 0.0;
                    double yinc = 0.0;
                    if (UserVariables.isUseCals()) {
                        xinc = alignmentParticle.getX() / UserVariables.getSpatialRes() - xc;
                        yinc = alignmentParticle.getY() / UserVariables.getSpatialRes() - yc;
                    }
                    virTemps[j].translate(-xinc, -yinc);
                    sigTemps[j].translate(-xinc, -yinc);
                    FloatProcessor sigSlice = new FloatProcessor(outputWidth, signalWidth);
                    FloatBlitter sigBlitter = new FloatBlitter(sigSlice);
                    sigBlitter.copyBits(sigTemps[j], 0, 0, Blitter.COPY);
                    output[1].addSlice(label, sigSlice);
                    FloatProcessor virSlice = new FloatProcessor(outputWidth, signalWidth);
                    FloatBlitter virBlitter = new FloatBlitter(virSlice);
                    virBlitter.copyBits(virTemps[j], 0, 0, Blitter.COPY);
                    output[0].addSlice(label, virSlice);
                }
            }
        }
//        progress.dispose();
        return output;
    }

    void extendSignalArea(float[] xpoints, float[] ypoints, float dist, int window) {
        xpoints[0] = xpoints[1];
        ypoints[0] = ypoints[1];
        int l = xpoints.length;
        xpoints[l - 1] = xpoints[l - 2];
        ypoints[l - 1] = ypoints[l - 2];
        if (xpoints.length - 1 <= window || ypoints.length - 1 <= window) {
            return;
        }
        getExtension(xpoints, ypoints, dist, window, 0, 1);
        getExtension(xpoints, ypoints, dist, window, l - 1, l - 2);
    }

    void getExtension(float[] xpoints, float[] ypoints, float dist, int window, int index1, int index2) {
        float xdiff = 0.0f, ydiff = 0.0f;
        int a, b, inc;
        if (index1 < index2) {
            a = 1;
            b = a + window;
            inc = 1;
        } else {
            b = xpoints.length;
            a = xpoints.length - window;
            inc = -1;
        }
        for (int i = a; i < b; i++) {
            xdiff += xpoints[i] - xpoints[i + inc];
            ydiff += ypoints[i] - ypoints[i + inc];
        }
        float newX, newY;
        if (xdiff != 0.0) {
            float ratio = Math.abs(ydiff / xdiff);
            newX = dist / (float) Math.sqrt(1.0f + (float) Math.pow(ratio, 2.0f));
            newY = newX * ratio;
        } else {
            newX = 0.0f;
            newY = ydiff;
        }
        if (xdiff < 0.0) {
            xpoints[index1] = xpoints[index2] - newX;
        } else if (xdiff > 0.0) {
            xpoints[index1] = xpoints[index2] + newX;
        }
        if (ydiff < 0.0) {
            ypoints[index1] = ypoints[index2] - newY;
        } else if (ydiff > 0.0) {
            ypoints[index1] = ypoints[index2] + newY;
        }
    }

    void goshtasbyShiftEval(ImageProcessor xcoeffs, ImageProcessor ycoeffs, ImageProcessor coords) {
        for (int y = 0; y < 512; y++) {
            for (int x = 0; x < 256; x++) {
                double xg = goshtasbyEval(xcoeffs, coords, x * UserVariables.getSpatialRes(), y * UserVariables.getSpatialRes());
                System.out.print(" " + (xg - x * UserVariables.getSpatialRes()) + " ");
            }
            System.out.println();
        }
        System.out.println();
        for (int y = 0; y < 512; y++) {
            for (int x = 0; x < 256; x++) {
                double yg = goshtasbyEval(ycoeffs, coords, x * UserVariables.getSpatialRes(), y * UserVariables.getSpatialRes());
                System.out.print(" " + (yg - y * UserVariables.getSpatialRes()) + " ");
            }
            System.out.println();
        }
    }

    void goshtasbyErrorEval() {
        TextReader reader = new TextReader();
        ImageProcessor tcoords = reader.open(calDir + delimiter + "testcoords.txt");
        ImageProcessor xcoeffs = reader.open(calDir + delimiter + "xcoeffs.txt");
        ImageProcessor ycoeffs = reader.open(calDir + delimiter + "ycoeffs.txt");
        ImageProcessor coords = reader.open(calDir + delimiter + "coords.txt");
        int size = tcoords.getHeight();
        System.out.println("x,y,xg,yg");
        for (int i = 1; i <= size; i++) {
            double x = tcoords.getPixelValue(0, i - 1);
            double y = tcoords.getPixelValue(1, i - 1);
            double xg = goshtasbyEval(xcoeffs, coords, x, y);
            double yg = goshtasbyEval(ycoeffs, coords, x, y);
            System.out.println(x + "," + y + "," + xg + "," + yg + "");
        }
    }

    void multiGoshtasbyErrorEval(int m, int n, int width, int height) {
        double xdiv = width * UserVariables.getSpatialRes() / m;
        double ydiv = height * UserVariables.getSpatialRes() / n;
        TextReader reader = new TextReader();
        ImageProcessor tcoords = reader.open(calDir + delimiter + "testcoords.txt");
        int size = tcoords.getHeight();
        System.out.println(m + "_" + n);
        System.out.println("x,y,xg,yg");
        for (int i = 1; i <= size; i++) {
            double x = tcoords.getPixelValue(0, i - 1);
            double y = tcoords.getPixelValue(1, i - 1);
            int xi = 1 + (int) Math.floor(x / xdiv);
            int yi = 1 + (int) Math.floor(y / ydiv);
            ImageProcessor xcoeffs = reader.open(calDir + delimiter + "goshtasby"
                    + delimiter + m + "_" + n + delimiter + "xcoeffs" + xi + "_" + yi + ".txt");
            ImageProcessor ycoeffs = reader.open(calDir + delimiter + "goshtasby"
                    + delimiter + m + "_" + n + delimiter + "ycoeffs" + xi + "_" + yi + ".txt");
            ImageProcessor coords = reader.open(calDir + delimiter + "goshtasby"
                    + delimiter + m + "_" + n + delimiter + "coords" + xi + "_" + yi + ".txt");
            double xg = goshtasbyEval(xcoeffs, coords, x, y);
            double yg = goshtasbyEval(ycoeffs, coords, x, y);
            System.out.println(x + "," + y + "," + xg + "," + yg + "");
        }
    }

    double goshtasbyEval(ImageProcessor coeffs, ImageProcessor coords, double x, double y) {
        int l = coeffs.getHeight();
        double sum = 0.0;
        for (int i = 3; i < l; i++) {
//            if (Utils.calcDistance(x, y, coords.getPixelValue(0, i - 3), coords.getPixelValue(1, i - 3)) < 20.0) {
            double r = Math.pow((x - coords.getPixelValue(0, i - 3)), 2.0) + Math.pow((y - coords.getPixelValue(1, i - 3)), 2.0);
            if (r > 0.0) {
                double R = r * Math.log(r);
                sum = sum + coeffs.getPixelValue(0, i) * R;
//                System.out.println("x," + coords.getPixelValue(0, i - 3) + ",y," + coords.getPixelValue(1, i - 3) + ",c," + (coeffs.getPixelValue(0, i) * R));
            }
//            }
        }
        return coeffs.getPixelValue(0, 0) + coeffs.getPixelValue(0, 1) * x + coeffs.getPixelValue(0, 2) * y + sum;
    }

    protected void printParams(String dir) {
        File paramFile;
        PrintWriter paramStream;
        try {
            paramFile = new File(dir + delimiter + "params.csv");
            paramStream = new PrintWriter(new FileOutputStream(paramFile));
        } catch (FileNotFoundException e) {
            System.out.println("Error: Failed to create parameter file.\n");
            System.out.println(e.toString());
            return;
        }
        paramStream.println(title);
        paramStream.println(Utilities.getDate("dd/MM/yyyy HH:mm:ss"));
        paramStream.println();
        paramStream.println(UserInterface.getRedSigEstText() + "," + UserVariables.getSigEstRed());
        paramStream.println(UserInterface.getGreenSigEstText() + "," + UserVariables.getSigEstGreen());
        paramStream.println(UserInterface.getSpatResLabelText() + "," + UserVariables.getSpatialRes());
        paramStream.println(UserInterface.getChan1MaxThreshLabelText() + "," + UserVariables.getChan1MaxThresh());
        paramStream.println(UserInterface.getChan2MaxThreshLabelText() + "," + UserVariables.getChan2MaxThresh());
        paramStream.println(UserInterface.getCurveFitTolLabelText() + "," + UserVariables.getCurveFitTol());
        paramStream.println(UserInterface.getPreprocessToggleText() + "," + UserVariables.isPreProcess());
        paramStream.println(UserInterface.getGpuToggleText() + "," + UserVariables.isGpu());
        paramStream.println(UserInterface.getFpsLabelText() + "," + UserVariables.getTimeRes());
        paramStream.println(UserInterface.getMinTrajLengthLabelText() + "," + UserVariables.getMinTrajLength());
        paramStream.println(UserInterface.getMinTrajDistLabelText() + "," + UserVariables.getMinTrajDist());
        paramStream.println(UserInterface.getMinTrajMSDLabelText() + "," + UserVariables.getMsdThresh());
        paramStream.println(UserInterface.getMaxLinkDistLabelText() + "," + UserVariables.getTrajMaxStep());
        paramStream.println(UserInterface.getTrackLengthText() + "," + UserVariables.getTrackLength());
        paramStream.println(UserInterface.getColocalToggleText() + "," + UserVariables.isColocal());
        paramStream.println(UserInterface.getPrevResToggleText() + "," + UserVariables.isPrevRes());
        paramStream.println(UserInterface.getExtractSigsToggleText() + "," + UserVariables.isExtractsigs());
        paramStream.println(UserInterface.getUseCalToggleText() + "," + UserVariables.isUseCals());
        paramStream.close();
    }

    public boolean isGpuEnabled() {
        return gpuEnabled;
    }

    protected IsoGaussian getC2Gaussian(double x0, double y0, ImageProcessor ip2) {
        if (ip2 == null) {
            return null;
        }
        int rx = (int) Math.round(x0 / UserVariables.getSpatialRes());
        int ry = (int) Math.round(y0 / UserVariables.getSpatialRes());
        int width = (int) Math.round(1.0 / UserVariables.getSpatialRes());
        int width2 = width * 2;
        Rectangle r = new Rectangle(rx - width, ry - width, width2, width2);
        IsoGaussian c2Gaussian = null;
        ip2.setRoi(r);
        ImageProcessor ip3 = ip2.crop();
        FloatStatistics stats1 = new FloatStatistics(ip3, ImageStatistics.MEDIAN, null);
        ip3.multiply(1.0 / stats1.median);
        FloatStatistics stats2 = new FloatStatistics(ip3, ImageStatistics.MIN_MAX, null);
        if (stats2.max > UserVariables.getChan2MaxThresh()) {
            c2Gaussian = new IsoGaussian(x0, y0, stats2.max * stats1.median, UserVariables.getSigEstGreen(), UserVariables.getSigEstGreen(), 0.0);
        }
        return c2Gaussian;
    }

    boolean checkTrajColocalisation(ParticleTrajectory traj, double thresh, boolean colocal) {
        if (!colocal) {
            return true;
        } else if (traj.getType(thresh) == ParticleTrajectory.COLOCAL) {
            return true;
        } else {
            return false;
        }
    }

    protected boolean getActiveImages() {
        ImageStack[] stacks = new ImageStack[2];
        if (IJ.getInstance() != null) {
            inputs = GenUtils.specifyInputs(labels);
            if (inputs == null) {
                return false;
            }
            stacks[0] = inputs[0].getImageStack();
            if (inputs[1] == null) {
            } else {
                stacks[1] = inputs[1].getImageStack();
            }
        }
        if (stacks[0].getProcessor(1).getNChannels() > 1 || (!(stacks[1] == null) && stacks[1].getProcessor(1).getNChannels() > 1)) {
            GenUtils.error("Monochrome images required.");
            return false;
        }
        if (!(stacks[1] == null) && stacks[0].getSize() != stacks[1].getSize()) {
            GenUtils.error("Stacks must have same number of slices.");
            return false;
        }
        return true;
    }

    void readParamsFromImage() {
        ParamsReader reader = new ParamsReader(inputs[0]);
        UserVariables.setSpatialRes(reader.getSpatialRes());
        UserVariables.setTimeRes(reader.getFrameRate());
    }

    public void printTrajectories(ArrayList<ParticleTrajectory> trajectories, File output, int length) throws IOException, FileNotFoundException {
        String headings[] = new String[]{"X", "Y", "Time (s)", "Channel 1", "Channel 2"};
        CSVPrinter printer = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(output), GenVariables.ISO), CSVFormat.EXCEL);
        int n = trajectories.size();
        Particle[][] particles = new Particle[n][length];
        for (int i = 0; i < n; i++) {
            ParticleTrajectory pt = trajectories.get(i);
            printer.print(String.format("Particle %d", i));
            for (int j = 1; j < headings.length; j++) {
                printer.print("");
            }
            Arrays.fill(particles[i], null);
            Particle current = pt.getEnd();
            while (current != null) {
                particles[i][current.getTimePoint()] = current;
                current = current.getLink();
            }
        }
        printer.println();
        for (int i = 0; i < n; i++) {
            for (String heading : headings) {
                printer.print(heading);
            }
        }
        printer.println();
        for (int t = 0; t < length; t++) {
            for (int i = 0; i < n; i++) {
                if (particles[i][t] != null) {
                    printer.print(particles[i][t].getX());
                    printer.print(particles[i][t].getY());
                    printer.print(t / UserVariables.getTimeRes());
                    printer.print(particles[i][t].getMagnitude());
                    Particle p2 = particles[i][t].getColocalisedParticle();
                    if (p2 != null) {
                        printer.print(p2.getMagnitude());
                    } else {
                        printer.print("");
                    }
                } else {
                    for (String heading : headings) {
                        printer.print("");
                    }
                }
            }
            printer.println();
        }
        printer.close();
    }
}
