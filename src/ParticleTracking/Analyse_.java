package ParticleTracking;

import ui.ResultsPreviewInterface;
import ui.UserInterface;
import IAClasses.IsoGaussian;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.Straightener;
import ij.plugin.TextReader;
import ij.plugin.filter.GaussianBlur;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

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
public class Analyse_ implements PlugIn {

//    protected static double hystDiff = 1.25;
    protected final double SIG_EST_RED = 0.158, SIG_EST_GREEN = 0.168;
    protected final double sigmas[] = {SIG_EST_RED, SIG_EST_GREEN};
//    protected int xyPartRad; //Radius over which to draw particles in visualisation
    public final int GOSHTASBY_M = 2, GOSHTASBY_N = 4;
    public static final int VERSION = 5;
//    protected final double LAMBDA = 650.0, //Wavelength of light
//            NUM_AP = 1.4; //Numerical aperture of system
    protected static double colocalThresh = 0.1;
    protected ArrayList<ParticleTrajectory> trajectories = new ArrayList(); //Trajectories of the detected particles
    protected ImagePlus imp; //The active image stack
    protected ImageStack stacks[];
    private long startTime;
    protected DecimalFormat numFormat = new DecimalFormat("0.000");
    protected DecimalFormat intFormat = new DecimalFormat("000");
    protected DecimalFormat floatFormat = new DecimalFormat("0.00");
    String title = "Particle Tracker", ext;
    protected static boolean msdPlot = false, intensPlot = false, trajPlot = false;
    protected boolean monoChrome, useCals = true;
//    private final double IMAGE_MAX = 255.0;
    protected final double SEARCH_SCALE = 0.5;
//    private final double TRACK_LENGTH = 5.0;
    private final double TRACK_WIDTH = 4.0;
    public static final float TRACK_EXT = 1.0f;
    public static final float TRACK_OFFSET = 0.75f;
    protected static File inputDir = new File("C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails"),
            c0Dir, c1Dir, outputDir,
            calDir = new File("C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails");
    private final String delimiter = GenUtils.getDelimiter();
    String parentDir;
    protected ImagePlus[] inputs;
    protected final String labels[] = {"Channel 1", "Channel 2"};
    protected boolean gpuEnabled = false;

//    public static void main(String args[]) {
////        if (imp != null) {
//        Analyse_ instance = new Analyse_();
//        instance.run(null);
////        }
//        System.exit(0);
//    }
    public Analyse_(double spatialRes, double timeRes, double trajMaxStep, double chan1MaxThresh, boolean monoChrome, ImagePlus imp, double scale, double minTrajLength) {
        UserVariables.setSpatialRes(spatialRes);
        UserVariables.setTimeRes(timeRes);
        UserVariables.setTrajMaxStep(trajMaxStep);
        UserVariables.setChan1MaxThresh(chan1MaxThresh);
//        Timelapse_Analysis.hystDiff = hystDiff;
        UserVariables.setMinTrajLength(minTrajLength);
        this.monoChrome = monoChrome;
//        this.imp = imp;
//        this.stacks = imp.getStack();
    }

    public Analyse_() {
    }

    public Analyse_(ImageStack[] stacks) {
        this.stacks = stacks;
    }

    public Analyse_(ImagePlus imp, String ext) {
//        this.imp = imp;
//        this.stacks = imp.getImageStack();
        this.ext = ext;
    }

    /**
     * Implements run method from {@link PlugIn}.
     */
    public void run(String arg) {
        Utilities.setLookAndFeel(UserInterface.class);
        title = title + "_v" + VERSION + "." + intFormat.format(Revision.Revision.revisionNumber);
        stacks = new ImageStack[2];
        if (IJ.getInstance() != null) {
            if (!getActiveImages()) {
                return;
            }
        } else {
            inputDir = Utilities.getFolder(inputDir, null, true);
            if (inputDir == null) {
                return;
            }
            c0Dir = new File(inputDir.getAbsolutePath() + delimiter + "C0");
            ImagePlus imp = Utils.buildStack(c0Dir);
            stacks[0] = imp.getImageStack();
//        Utils.normaliseStack(stacks[0], IMAGE_MAX);
            this.ext = imp.getTitle();
            c1Dir = new File(inputDir.getAbsolutePath() + delimiter + "C1");
            if (c1Dir.exists()) {
                stacks[1] = (Utils.buildStack(c1Dir)).getImageStack();
//            Utils.normaliseStack(stacks[1], IMAGE_MAX);
            }
        }
//        if (IJ.getInstance() != null) {
//            imp = WindowManager.getCurrentImage();
//        }
        if (showDialog()) {
            analyse();
        }
    }

    public boolean showDialog() {
//        if (imp == null) {
//            Toolkit.getDefaultToolkit().beep();
//            IJ.error("No image stack open.");
//            return false;
//        }
//        stack = imp.getImageStack();
        monoChrome = (stacks[1] == null);
        UserInterface ui = new UserInterface(null, true, title, this);
        ui.setVisible(true);
        return ui.isWasOKed();
    }

    /**
     * Analyses the {@link ImageStack} specified by <code>stack</code>.
     */
    public void analyse() {
        outputDir = Utilities.getFolder(inputDir, "Specify directory for output files...", true);
        parentDir = GenUtils.openResultsDirectory(outputDir + delimiter + title, delimiter);
        String sigc0Dir = GenUtils.openResultsDirectory(parentDir + delimiter + "C0", delimiter);
        String sigc1Dir = GenUtils.openResultsDirectory(parentDir + delimiter + "C1", delimiter);
        if (!monoChrome && useCals) {
            calDir = Utilities.getFolder(outputDir, "Specify directory containing calibrations...", true);
        }
        if (stacks != null) {
            IJ.register(this.getClass());
            startTime = System.currentTimeMillis();
            //gaussianRadius = 0.139d / spatialRes; // IsoGaussian filter radius set to 139 nm
            calcParticleRadius(UserVariables.getSpatialRes(), sigmas[UserVariables.getC1Index()]);
            int i, count;
//            int width = stacks[0].getWidth(), height = stacks[0].getHeight();
            findParticles();
            TextWindow results = new TextWindow(title + " Results", "X\tY\tFrame\tChannel 1 ("
                    + UserVariables.channels[UserVariables.getC1Index()]
                    + ")\tChannel 2 (" + UserVariables.channels[UserVariables.getC2Index()]
                    + ")\tChannel 2 " + '\u03C3' + "x\tChannel 2 " + '\u03C3' + "y\t" + '\u03B8',
                    new String(), 1000, 500);
//            results.append(imp.getTitle() + "\n\n");
            TextWindow resultSummary = new TextWindow(title + " Results Summary",
                    "Particle\tType\t% Colocalisation\tDuration (s)\tDisplacement (" + IJ.micronSymbol
                    + "m)\tVelocity (" + IJ.micronSymbol + "m/s)\tDirectionality\tDiffusion Coefficient ("
                    + IJ.micronSymbol + "m^2/s)" + "\tFractal Dimension"
                    + "\tFluorescence Ratio ("
                    + UserVariables.channels[UserVariables.getC2Index()] + "/"
                    + UserVariables.channels[UserVariables.getC1Index()]
                    + ")\tAngle Spread\tStep Spread\tDC\tCurvature\tC2 Fluor Area\tC2 Fluor Skew",
                    new String(), 1200, 500);
//            resultSummary.append(imp.getTitle() + "\n\n");

            int n = trajectories.size();
            for (i = 0; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (!(traj.getSize() / UserVariables.getTimeRes() > UserVariables.getMinTrajLength()
                        && traj.getDisplacement(traj.getEnd(), traj.getSize()) > UserVariables.getMinTrajDist()
                        && ((traj.getType(colocalThresh) == ParticleTrajectory.COLOCAL)
                        || ((traj.getType(colocalThresh) == ParticleTrajectory.NON_COLOCAL) && !UserVariables.isColocal())))) {
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
            ImageStack maps = mapTrajectories((new RGBStackMerge()).mergeStacks(stacks[0].getWidth(), stacks[0].getHeight(), stacks[0].getSize(), stacks[0], stacks[1], null, true),
                    trajectories, UserVariables.getSpatialRes(), UserVariables.getMinTrajLength(),
                    UserVariables.getTimeRes(), true, 0, trajectories.size() - 1, 1, false);
            for (i = 0, count = 1; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (traj != null) {
                    int type = traj.getType(colocalThresh);
                    if (traj.getSize() > UserVariables.getMinTrajLength() && ((type == ParticleTrajectory.COLOCAL)
                            || ((type == ParticleTrajectory.NON_COLOCAL) && !UserVariables.isColocal()))) {
//                        if (intensPlot) {
//                            plotIntensity(i, count);
//                        }
//                        if (trajPlot) {
//                            plotTrajectory(width, height, i, count);
//                        }
                        printData(i, resultSummary, count);
                        traj.printTrajectory(count, results, numFormat, title);
                        if (type == ParticleTrajectory.COLOCAL) {
                            ImageStack signals[] = extractTrajSignalValues(traj,
                                    (int) Math.round(UserVariables.getTrackLength() / UserVariables.getSpatialRes()),
                                    (int) Math.round(TRACK_WIDTH / UserVariables.getSpatialRes()),
                                    TRACK_EXT / ((float) UserVariables.getSpatialRes()), stacks[0].getWidth(), stacks[0].getHeight(), count);
                            if (signals[0].getSize() > 0) {
                                for (int j = 1; j <= signals[0].getSize(); j++) {
                                    IJ.saveAs((new ImagePlus("", signals[0].getProcessor(j))),
                                            "TIF", sigc0Dir + delimiter + "C0-" + String.valueOf(count)
                                            + "-" + signals[0].getSliceLabel(j));
                                    IJ.saveAs((new ImagePlus("", signals[1].getProcessor(j))),
                                            "TIF", sigc1Dir + delimiter + "C1-" + String.valueOf(count)
                                            + "-" + signals[1].getSliceLabel(j));
                                }
                            }
                        }
                        count++;
                    }
                }
            }
            resultSummary.append("\nAnalysis Time (s): " + numFormat.format((System.currentTimeMillis() - startTime) / 1000.0));
            results.append(toString());
            results.setVisible(true);
            resultSummary.setVisible(true);
            IJ.saveString(results.getTextPanel().getText(), parentDir + "/results.txt");
            IJ.saveString(resultSummary.getTextPanel().getText(), parentDir + "/resultsSummary.txt");
            if (maps != null) {
                (new ImagePlus("Trajectory Maps", maps)).show();
                IJ.saveAs((new ImagePlus("", maps)), "TIF", parentDir + "/trajectories.tif");
            }
        }
        printParams(parentDir);
    }

    protected ParticleArray findParticles() {
        return findParticles(SEARCH_SCALE, true, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1], false, sigmas[UserVariables.getC1Index()], sigmas[1 - UserVariables.getC1Index()], UserVariables.isColocal(), true, true, UserVariables.getC2CurveFitTol());
    }

    /**
     * Median filter and IsoGaussian filter the image specified by
     * <code>processor</code>.
     *
     * @param processor the image to be pre-processed.
     */
    public FloatProcessor preProcess(ImageProcessor processor) {
        if (processor == null) {
            return null;
        }
        FloatProcessor fp;
        if (UserVariables.isPreProcess()) {
            TypeConverter tc = new TypeConverter(processor, false);
            fp = (FloatProcessor) tc.convertToFloat(null);
            (new GaussianBlur()).blur(fp, sigmas[UserVariables.getC1Index()]);
        } else {
            TypeConverter tc = new TypeConverter(processor, false);
            fp = (FloatProcessor) tc.convertToFloat(null);
        }
        return fp;
    }

    public ParticleArray findParticles(boolean update, int startSlice, int endSlice, double c1FitTol, ImageStack channel1, ImageStack channel2, boolean fitC2Gaussian, boolean colocal, boolean floatingSigma, double c2FitTol) {
        return findParticles(SEARCH_SCALE, update, startSlice, endSlice, c1FitTol,
                channel1, channel2, fitC2Gaussian, sigmas[UserVariables.getC1Index()], sigmas[1 - UserVariables.getC1Index()], colocal, true, floatingSigma, c2FitTol);
    }

    public ParticleArray findParticles(double searchScale, boolean update, int startSlice, int endSlice, double c1FitTol, ImageStack channel1, ImageStack channel2, boolean fitC2Gaussian, double sigEst1, double sigEst2, boolean colocal, boolean showProgress, boolean floatingSigma, double c2FitTol) {
        if (channel1 == null) {
            return null;
        }
        int i, noOfImages = channel1.getSize(), width = channel1.getWidth(), height = channel1.getHeight(),
                arraySize = endSlice - startSlice + 1;
        int xyPartRad = calcParticleRadius(UserVariables.getSpatialRes(), sigEst1);
        int fitRad = (int) Math.ceil(xyPartRad);
        int c1X, c1Y, pSize = 2 * fitRad + 1;
        int c2Points[][];
        int radius = (int) Math.round(fitRad * searchScale);
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
//        ImageStack detect_output = new ImageStack(stack.getWidth(), stack.getHeight());
//        ImageStack maxima = new ImageStack(stack.getWidth(), stack.getHeight());
//        ImageStack input_output = new ImageStack(stack.getWidth(), stack.getHeight());
        ProgressDialog progress = new ProgressDialog(null, "Finding Particles...", false, title, false);
        if (showProgress) {
            progress.setVisible(true);
        }
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
//            ByteProcessor oslice = new ByteProcessor(detect_output.getWidth(), detect_output.getHeight());
            IJ.freeMemory();
            progress.updateProgress(i - startSlice, arraySize);
            FloatProcessor chan1Proc = preProcess(channel1.getProcessor(i + 1).duplicate());
            Utils.normalise(chan1Proc, 1.0);
            FloatProcessor chan2Proc = (channel2 != null) ? preProcess(channel2.getProcessor(i + 1).duplicate()) : null;
            if (chan2Proc != null) {
                Utils.normalise(chan2Proc, 1.0);
            }
            double c1Threshold = Utils.getPercentileThresh(chan1Proc, UserVariables.getChan1MaxThresh());
            ByteProcessor thisC1Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, UserVariables.FOREGROUND, chan1Proc, c1Threshold, true);
//            maxima.addSlice(thisC1Max);
            double c2Threshold = Utils.getPercentileThresh(chan2Proc, UserVariables.getChan2MaxThresh());
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (thisC1Max.getPixel(c1X, c1Y) == UserVariables.FOREGROUND) {
                        IsoGaussian c2Gaussian = null;
                        c2Points = Utils.searchNeighbourhood(c1X, c1Y, radius,
                                c2Threshold, (ImageProcessor) chan2Proc);
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
//                        MultiGaussFitter c1Fitter = new MultiGaussFitter(UserVariables.getnMax(), fitRad, pSize);
//                        c1Fitter.fit(pixValues, sigEst1 / UserVariables.getSpatialRes());
//                        ArrayList<IsoGaussian> c1Fits = c1Fitter.getFits(spatialRes, c1X - fitRad, c1Y - fitRad, c1Threshold, fitTol);
                        IsoGaussianFitter c1Fitter = new IsoGaussianFitter(xCoords, yCoords, pixValues, false);
                        c1Fitter.doFit(sigEst1 / UserVariables.getSpatialRes());
                        ArrayList<IsoGaussian> c1Fits = new ArrayList();
                        c1Fits.add(new IsoGaussian((c1X - fitRad + c1Fitter.getX0()) * spatialRes, (c1Y - fitRad + c1Fitter.getY0()) * spatialRes,
                                c1Fitter.getMag(), sigEst1 / UserVariables.getSpatialRes(), sigEst1 / UserVariables.getSpatialRes(), c1Fitter.getRSquared()));
                        if (c2Points != null) {
                            c2Gaussian = findC2Particle(c1X, c1Y, radius, pSize, c2Threshold, chan2Proc, fitC2Gaussian, sigEst2, fitRad, spatialRes, floatingSigma);
                        }

                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        for (IsoGaussian c1Fit : c1Fits) {
                            if (c1Fit.getFit() > c1FitTol
                                    && (c2Gaussian == null && !colocal)
                                    || (c2Gaussian != null && colocal
                                    && c2Gaussian.getFit() > c2FitTol)) {
                                particles.addDetection(i - startSlice, new Particle(i - startSlice, c1Fit, c2Gaussian, null, -1));
                            }
//                                Utils.draw2DGaussian(oslice, c1Fit, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), false, false);
//                                Utils.draw2DGaussian(chan1Proc, c1Fit, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), false, true);
                        }
                    }
                }
            }

//            detect_output.addSlice(oslice);
//            input_output.addSlice(chan1Proc.duplicate());
        }
        progress.dispose();
//        IJ.saveAs(new ImagePlus("", detect_output), "TIF", parentDir + "/all_detections.tif");
//        IJ.saveAs(new ImagePlus("", maxima), "TIF", parentDir + "/all_maxima.tif");
//        IJ.saveAs(new ImagePlus("", input_output), "TIF", parentDir + "/input_output.tif");
        if (update) {
            TrajectoryBuilder.updateTrajectories(particles, UserVariables.getTimeRes(), UserVariables.getTrajMaxStep(), spatialRes, true, 1.0, trajectories);
        }
        return particles;
    }

    IsoGaussian findC2Particle(int c1X, int c1Y, int radius, int pSize, double c2Threshold, FloatProcessor chan2Proc, boolean fitGaussian, double sigEst2, int fitRad, double spatialRes, boolean floatingSigma) {
        IsoGaussian c2Gaussian = null;
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        int c2Points[][] = Utils.searchNeighbourhood(c1X, c1Y, radius, c2Threshold, (ImageProcessor) chan2Proc);
        if (c2Points != null) {
            if (fitGaussian) {
                Utils.extractValues(xCoords, yCoords, pixValues, c2Points[0][0], c2Points[0][1], chan2Proc);
                IsoGaussianFitter c2Fitter = new IsoGaussianFitter(xCoords, yCoords, pixValues, floatingSigma);
                c2Fitter.doFit(sigEst2 / UserVariables.getSpatialRes());
                c2Gaussian = new IsoGaussian((c2Points[0][0] - fitRad + c2Fitter.getX0()) * spatialRes, (c2Points[0][1] - fitRad + c2Fitter.getY0()) * spatialRes,
                        c2Fitter.getMag(), c2Fitter.getXsig(), c2Fitter.getYsig(), c2Fitter.getRSquared());
            } else {
                c2Gaussian = new IsoGaussian(c2Points[0][0] * spatialRes, c2Points[0][1] * spatialRes,
                        chan2Proc.getPixelValue(c2Points[0][0], c2Points[0][1]), sigEst2, sigEst2, 0.0);
            }
        }
        return c2Gaussian;
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
     * <code>D = 1 / (1 + &lambda<sub>1</sub>
     * &lambda<sub>2</sub><sup>-1</sup>)</code>
     * <br> <br> where <code>&lambda<sub>1</sub></code> and
     * <code>&lambda<sub>2</sub></code> are the eigenvalues of the trajectory
     * data and      <code>&lambda<sub>1</sub> < &lambda<sub>2</sub></code>.
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
//        traj.smooth();
        double points[][] = traj.getPoints();
        traj.calcMSD(label, -1, msdPlot, points[0], points[1]);
//        traj.calcAngleSpread();
//        traj.calcStepSpread();
        traj.calcDirectionality(points[0], points[1]);
        double displacement = traj.getDisplacement(traj.getEnd(), traj.getSize());
        double duration = traj.getDuration();
        int type = traj.getType(colocalThresh);
        String trajType = null;
        switch (type) {
            case ParticleTrajectory.COLOCAL:
                trajType = "Colocalised";
                break;
            case ParticleTrajectory.NON_COLOCAL:
                trajType = "Non-Colocalised";
                break;
            case ParticleTrajectory.UNKNOWN:
                trajType = "Unknown";
        }
        output.append(label + "\t" + trajType + "\t"
                + decFormat.format(traj.getDualScore() * 100.0 / traj.getSize()) + "\t"
                + decFormat.format(duration) + "\t"
                + decFormat.format(displacement)
                + "\t" + decFormat.format(displacement / duration) + "\t"
                + decFormat.format(traj.getDirectionality()) + "\t"
                + msdFormat.format(traj.getDiffCoeff()) + "\t"
                + decFormat.format(traj.getBoxCountFD()) + "\t"
                + decFormat.format(traj.getFluorRatio()) + "\t"
                + decFormat.format(traj.getAngleSpread()) + "\t"
                + decFormat.format(traj.getStepSpread()) + "\t"
                + decFormat.format(traj.getLogDC()) + "\t"
                + decFormat.format(traj.getMeanKappa()) + "\t");
        return true;
    }

//    /**
//     * Produces a {@link Plot} of normalised intensity in the red and green
//     * channels, together with a ratio of red:green intensity, for the
//     * trajectory denoted by <code>particleNumber</code>.
//     */
//    public boolean plotIntensity(int particleNumber, int label) {
//        if (trajectories.size() < 1) {
//            return false;
//        }
//        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
//        Particle current = traj.getEnd();
//        int size = traj.getSize();
//        double xvalues[] = new double[size];
//        double redValues[] = new double[size];
//        double greenValues[] = new double[size];
//        double ratios[] = new double[size];
//        double temp, maxVal = -Double.MAX_VALUE, minVal = Double.MAX_VALUE;
//
//        if (traj.getType(colocalThresh) != ParticleTrajectory.COLOCAL) {
//            return false;
//        }
//
//        for (int i = size - 1; i >= 0; i--) {
//            xvalues[i] = current.getTimePoint() * UserVariables.getTimeRes();
//            redValues[i] = current.getC1Intensity() / 255.0d;
//            greenValues[i] = current.getC2Intensity() / 255.0d;
//            ratios[i] = greenValues[i] / redValues[i];
//            if (ratios[i] > maxVal) {
//                maxVal = ratios[i];
//            }
//            temp = Math.min(redValues[i], greenValues[i]);
//            if (temp < minVal) {
//                minVal = temp;
//            }
//            current = current.getLink();
//        }
//
//        Plot plot = new Plot("Particle " + label + " Intensities",
//                "Time (s)", "Normalised Intensity", xvalues, ratios,
//                (Plot.X_TICKS + Plot.Y_TICKS + Plot.X_NUMBERS + Plot.Y_NUMBERS));
//        plot.changeFont(new Font("Serif", Font.BOLD, 14));
//        plot.setLimits(0.0, (stacks.getSize() - 1.0) * UserVariables.getTimeRes(), minVal, maxVal);
//        plot.setLineWidth(2);
//        plot.setColor(Color.BLUE);
//        plot.draw();
//        plot.setColor(Color.RED);
//        plot.addPoints(xvalues, redValues, Plot.LINE);
//        plot.setColor(Color.GREEN);
//        plot.addPoints(xvalues, greenValues, Plot.LINE);
//        plot.show();
//
//        return true;
//    }
    /**
     * Constructed trajectories are drawn onto the original image sequence and
     * displayed as a stack sequence.
     */
    public ImageStack mapTrajectories(ImageStack stack, ArrayList<ParticleTrajectory> trajectories, double spatialRes, double minTrajLength, double timeRes, boolean tracks, int startT, int endT, int index, boolean preview) {
        if (stack == null) {
            return null;
        }
        int radius = calcParticleRadius(spatialRes, sigmas[UserVariables.getC1Index()]);
        int i, j, width = stack.getWidth(), height = stack.getHeight(),
                type, frames = stack.getSize();
        double lastX, lastY;
        ImageStack outputStack = new ImageStack(width, height);
        Particle current;
        ParticleTrajectory traj;
        int length, n = trajectories.size();
        ImageProcessor processor;
//        Rectangle bounds;
//        ByteProcessor trajImage;
        int lastTP;
//        double ptScale = ParticleTrajectory.scale;

        if (n < 1) {
            return null;
        }

        for (i = 0; i < frames; i++) {
            if (monoChrome) {
                processor = (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte());
            } else {
                processor = (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToRGB());
            }
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
                length = traj.getSize();
                type = traj.getType(colocalThresh);
//                bounds = traj.getBounds();
//                trajImage = new ByteProcessor(bounds.width, bounds.height);
//                trajImage.setColor(Color.white);
//                trajImage.fill();
//                trajImage.setColor(Color.black);
                if (length > minTrajLength && ((type == ParticleTrajectory.COLOCAL)
                        || ((type == ParticleTrajectory.NON_COLOCAL) && !UserVariables.isColocal()))) {
                    current = traj.getEnd();
                    lastX = current.getX();
                    lastY = current.getY();
                    lastTP = current.getTimePoint();
                    current = current.getLink();
                    while (current != null) {
                        for (j = frames - 1; j >= lastTP; j--) {
                            processor = outputStack.getProcessor(j + 1);
                            if (!monoChrome && !preview) {
                                processor.setColor(thiscolor);
                            } else {
                                processor.setColor(Color.white);
                            }
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
//                        if (tracks) {
//                            trajImage.drawLine((int) Math.round(ptScale * current.getX() - traj.getBounds().x),
//                                    (int) Math.round(ptScale * current.getY() - traj.getBounds().y),
//                                    (int) Math.round(ptScale * lastX - traj.getBounds().x),
//                                    (int) Math.round(ptScale * lastY - traj.getBounds().y));
//                        }
                        lastX = current.getX();
                        lastY = current.getY();
                        lastTP = current.getTimePoint();
                        current = current.getLink();
                    }
                    processor = outputStack.getProcessor(lastTP + 1);
                    if (!monoChrome && !preview) {
                        processor.setColor(thiscolor);
                    } else {
                        processor.setColor(Color.white);
                    }
                    markParticle(processor, (int) Math.round(lastX / spatialRes) - radius,
                            (int) Math.round(lastY / spatialRes) - radius, radius, true, "" + index);
                    index++;
                }
            }
        }
        progress.dispose();
        return outputStack;
    }

    void markParticle(ImageProcessor processor, int x, int y, int radius, boolean string, String label) {
        processor.drawRect(x, y, 2 * radius + 1, 2 * radius + 1);
//        processor.drawOval(x, y, 2 * radius, 2 * radius);
        if (string) {
            processor.drawString(label, x, y);
        }
    }

    public int calcParticleRadius(double spatialRes) {
        return calcParticleRadius(spatialRes, sigmas[UserVariables.getC1Index()]);
    }

    public int calcParticleRadius(double spatialRes, double sigEst) {
//        double airyRad = 1.22 * LAMBDA / (2.0 * NUM_AP); //Airy radius
//        xyPartRad = (int) Math.ceil(2.0*airyRad / (spatialRes * 1000.0));
//        sigEst = airyRad / (2.0 * spatialRes * 1000.0);
        return (int) Math.ceil(3.0 * sigEst / spatialRes);
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

    public Color getDrawColor(int key) {
        switch (key) {
            case UserVariables.RED:
                return Color.red;
            case UserVariables.GREEN:
                return Color.green;
//            case UserVariables.BLUE:
//                return Color.blue;
            default:
                return Color.white;
        }
    }

    public ArrayList<ParticleTrajectory> getTrajectories() {
        return trajectories;
    }

    public ImageStack[] getStacks() {
        return stacks;
    }

    public boolean isMonoChrome() {
        return monoChrome;
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
        ProgressDialog progress = new ProgressDialog(null, "Extracting Signal Areas " + count + "...", false, title, false);
        progress.setVisible(true);
        for (int i = 0; i < size; i++) {
            progress.updateProgress(i, size);
            Particle current = sigStartP;
            next = sigStartP.getLink();
            xSigPoints = new ArrayList();
            ySigPoints = new ArrayList();
            xVirPoints = new ArrayList();
            yVirPoints = new ArrayList();
            double x1, y1, x2, y2, t = 2.0 / UserVariables.getTimeRes();
            if (last != null) {
                x1 = last.getX();
                y1 = last.getY();
            } else {
                x1 = sigStartP.getX();
                y1 = sigStartP.getY();
                t /= 2.0;
            }
            if (next != null) {
                x2 = next.getX();
                y2 = next.getY();
            } else {
                x2 = sigStartP.getX();
                y2 = sigStartP.getY();
                t /= 2.0;
            }
            double vel = Utils.calcDistance(x1, y1, x2, y2) / t;
            for (int index = 1; index <= size && current != null; index++) {
                double xg = current.getC1Gaussian().getX();
                double yg = current.getC1Gaussian().getY();
                if (useCals) {
                    double x = current.getC1Gaussian().getX(), y = current.getC1Gaussian().getY();
                    int xi = 1 + (int) Math.floor(x / xdiv);
                    int yi = 1 + (int) Math.floor(y / ydiv);
                    ImageProcessor xcoeffs = reader.open(calDir + delimiter + "goshtasby"
                            + delimiter + GOSHTASBY_M + "_" + GOSHTASBY_N + delimiter + "xcoeffs" + xi + "_" + yi + ".txt");
                    ImageProcessor ycoeffs = reader.open(calDir + delimiter + "goshtasby"
                            + delimiter + GOSHTASBY_M + "_" + GOSHTASBY_N + delimiter + "ycoeffs" + xi + "_" + yi + ".txt");
                    ImageProcessor coords = reader.open(calDir + delimiter + "goshtasby"
                            + delimiter + GOSHTASBY_M + "_" + GOSHTASBY_N + delimiter + "coords" + xi + "_" + yi + ".txt");
                    xg = goshtasbyEval(xcoeffs, coords, x, y);
                    yg = goshtasbyEval(ycoeffs, coords, x, y);
                }
                xSigPoints.add((float) (xg / UserVariables.getSpatialRes()));
                ySigPoints.add((float) (yg / UserVariables.getSpatialRes()));
                xVirPoints.add((float) (current.getC1Gaussian().getX() / UserVariables.getSpatialRes()));
                yVirPoints.add((float) (current.getC1Gaussian().getY() / UserVariables.getSpatialRes()));
                current = current.getLink();
            }
            xSigArray = new float[xSigPoints.size() + 1];
            ySigArray = new float[ySigPoints.size() + 1];
            xVirArray = new float[xVirPoints.size() + 1];
            yVirArray = new float[yVirPoints.size() + 1];
            for (int j = 1; j <= xSigPoints.size(); j++) {
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
        progress.dispose();
        int xc = (int) Math.ceil(TRACK_OFFSET * offset);
        int yc = (signalWidth - 1) / 2;
        int outputWidth = (int) Math.round(signalLength + offset);
        ImageStack output[] = new ImageStack[2];
        output[0] = new ImageStack(outputWidth, signalWidth);
        output[1] = new ImageStack(outputWidth, signalWidth);

        progress = new ProgressDialog(null, "Aligning Signal Areas " + count + "...", false, title, false);
        progress.setVisible(true);
        for (int j = 0; j < size; j++) {
            progress.updateProgress(j, size);
            if (virTemps[j] != null && sigTemps[j] != null) {
                Particle alignmentParticle = null;
                if (useCals) {
                    ImageStack virStack = new ImageStack(virTemps[j].getWidth(), virTemps[j].getHeight());
                    virStack.addSlice(virTemps[j]);
                    ParticleArray particles = findParticles(0.0, false, 0, 0, 0.0, virStack, null, true, sigmas[UserVariables.getC1Index()], sigmas[1 - UserVariables.getC1Index()], false, false, false, 0.0);
                    ArrayList<Particle> detections = particles.getLevel(0);
                    double mindist = Double.MAX_VALUE;
                    int minindex = -1;
                    for (int k = 0; k < detections.size(); k++) {
                        Particle p = detections.get(k);
                        double dist = Utils.calcDistance(p.getX(), p.getY(), xc * UserVariables.getSpatialRes(), yc * UserVariables.getSpatialRes());
                        if (dist < mindist) {
                            mindist = dist;
                            minindex = k;
                        }
                    }
                    if (minindex > -1) {
                        alignmentParticle = detections.get(minindex);
                    }
                }
                if (!useCals || alignmentParticle != null) {
                    String label = Float.toString(virTemps[j].getPixelValue(0, 0))
                            + "-" + floatFormat.format(virTemps[j].getPixelValue(1, 0));
                    virTemps[j].setInterpolate(true);
                    virTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    sigTemps[j].setInterpolate(true);
                    sigTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    double xinc = 0.0;
                    double yinc = 0.0;
                    if (useCals) {
                        xinc = alignmentParticle.getC1Gaussian().getX() / UserVariables.getSpatialRes() - xc;
                        yinc = alignmentParticle.getC1Gaussian().getY() / UserVariables.getSpatialRes() - yc;
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
        progress.dispose();
        return output;
    }

    ImageStack[] extractStaticSignalValues(ParticleTrajectory ptraj, int signalWidth) {
        TextReader reader = new TextReader();
        ImageProcessor xcoeffs = null, ycoeffs = null, coords = null;
        if (useCals) {
            xcoeffs = reader.open(calDir + delimiter + "xcoeffs.txt");
            ycoeffs = reader.open(calDir + delimiter + "ycoeffs.txt");
            coords = reader.open(calDir + delimiter + "coords.txt");
        }
        Particle sigStartP = ptraj.getEnd();
        if (signalWidth % 2 == 0) {
            signalWidth++;
        }
        int size = ptraj.getSize();
        ImageProcessor[] sigTemps = new ImageProcessor[size];
        ImageProcessor[] virTemps = new ImageProcessor[size];
        Particle current = sigStartP;
        for (int f = size - 1; f >= 0; f--) {
            double xr = current.getX();
            double yr = current.getY();
            double xg = xr;
            double yg = yr;
            if (useCals) {
                xg = goshtasbyEval(xcoeffs, coords, xr, yr);
                yg = goshtasbyEval(ycoeffs, coords, xr, yr);
            }
            ImageProcessor sigIP = stacks[1].getProcessor(current.getTimePoint() + 1);
            ImageProcessor virIP = stacks[0].getProcessor(current.getTimePoint() + 1);
            FloatProcessor sigRegion = new FloatProcessor(2 * signalWidth + 1, 2 * signalWidth + 1);
            FloatProcessor virRegion = new FloatProcessor(2 * signalWidth + 1, 2 * signalWidth + 1);
            sigIP.setInterpolate(true);
            sigIP.setInterpolationMethod(ImageProcessor.BICUBIC);
            virIP.setInterpolate(true);
            virIP.setInterpolationMethod(ImageProcessor.BICUBIC);
            for (int j = 0; j < sigRegion.getHeight(); j++) {
                for (int i = 0; i < sigRegion.getWidth(); i++) {
                    sigRegion.putPixelValue(i, j, sigIP.getInterpolatedValue(xg / UserVariables.getSpatialRes() - signalWidth + i,
                            yg / UserVariables.getSpatialRes() - signalWidth + j));
                    virRegion.putPixelValue(i, j, virIP.getInterpolatedValue(xr / UserVariables.getSpatialRes() - signalWidth + i,
                            yr / UserVariables.getSpatialRes() - signalWidth + j));
                }
            }
            sigTemps[f] = sigRegion;
            virTemps[f] = virRegion;
            IJ.saveAs(new ImagePlus("", virRegion), "TIF", "C:\\Users\\barry05\\Desktop\\virRegions\\" + current.getTimePoint() + ".tif");
            if (virTemps[f] != null) {
                virTemps[f].putPixelValue(0, 0, current.getTimePoint());
            }
            current = current.getLink();
        }
        int xc = signalWidth;
        int yc = signalWidth;
        ImageStack output[] = new ImageStack[2];
        output[0] = new ImageStack(2 * signalWidth + 1, 2 * signalWidth + 1);
        output[1] = new ImageStack(2 * signalWidth + 1, 2 * signalWidth + 1);
        for (int j = 0; j < size; j++) {
            if (virTemps[j] != null && sigTemps[j] != null) {
                ParticleArray particles = null;
                if (useCals) {
                    ImageStack virStack = new ImageStack(virTemps[j].getWidth(), virTemps[j].getHeight());
                    virStack.addSlice(virTemps[j]);
                    particles = findParticles(0.0, false, 0, 0, 0.0, virStack, null, true, sigmas[UserVariables.getC1Index()], sigmas[1 - UserVariables.getC1Index()], false, true, false, 0.0);
                }
                if (!useCals || !particles.getLevel(0).isEmpty()) {
                    String timepoint = Float.toString(virTemps[j].getPixelValue(0, 0));
                    virTemps[j].setInterpolate(true);
                    virTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    sigTemps[j].setInterpolate(true);
                    sigTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    double xinc = 0.0;
                    double yinc = 0.0;
                    if (useCals) {
                        Particle p = particles.getLevel(0).get(0);
                        xinc = p.getC1Gaussian().getX() / UserVariables.getSpatialRes() - xc;
                        yinc = p.getC1Gaussian().getY() / UserVariables.getSpatialRes() - yc;
                    }
                    virTemps[j].translate(-xinc, -yinc);
                    sigTemps[j].translate(-xinc, -yinc);
                    FloatProcessor sigSlice = new FloatProcessor(2 * signalWidth + 1, 2 * signalWidth + 1);
                    FloatBlitter sigBlitter = new FloatBlitter(sigSlice);
                    sigBlitter.copyBits(sigTemps[j], 0, 0, Blitter.COPY);
                    output[1].addSlice(timepoint, sigSlice);
                    FloatProcessor virSlice = new FloatProcessor(2 * signalWidth + 1, 2 * signalWidth + 1);
                    FloatBlitter virBlitter = new FloatBlitter(virSlice);
                    virBlitter.copyBits(virTemps[j], 0, 0, Blitter.COPY);
                    output[0].addSlice(timepoint, virSlice);
                }
            }
        }
        return output;
    }

    void extendSignalArea(float[] xpoints, float[] ypoints, float dist, int window) {
        float xdiff = 0.0f, ydiff = 0.0f;
        xpoints[0] = xpoints[1];
        ypoints[0] = ypoints[1];
        if (xpoints.length - 1 <= window || ypoints.length - 1 <= window) {
            return;
        }
        for (int i = 1; i <= window; i++) {
            xdiff += xpoints[i] - xpoints[i + 1];
            ydiff += ypoints[i] - ypoints[i + 1];
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
            xpoints[0] = xpoints[1] - newX;
        } else if (xdiff > 0.0) {
            xpoints[0] = xpoints[1] + newX;
        }
        if (ydiff < 0.0) {
            ypoints[0] = ypoints[1] - newY;
        } else if (ydiff > 0.0) {
            ypoints[0] = ypoints[1] + newY;
        }
    }

//    public int getXyPartRad() {
//        return xyPartRad;
//    }
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

    void printParams(String dir) {
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
        paramStream.println(UserInterface.getChannel1LabelText() + "," + UserVariables.getC1Index());
        paramStream.println(UserInterface.getChannel2LabelText() + "," + UserVariables.getC2Index());
        paramStream.println(UserInterface.getSpatResLabelText() + "," + UserVariables.getSpatialRes());
        paramStream.println(UserInterface.getFpsLabelText() + "," + UserVariables.getTimeRes());
        paramStream.println(UserInterface.getMinTrajLengthLabelText() + "," + UserVariables.getMinTrajLength());
        paramStream.println(UserInterface.getMinTrajDistLabelText() + "," + UserVariables.getMinTrajDist());
        paramStream.println(UserInterface.getMaxLinkDistLabelText() + "," + UserVariables.getTrajMaxStep());
        paramStream.println(UserInterface.getTrackLengthText() + "," + UserVariables.getTrackLength());
        paramStream.println(UserInterface.getChan1MaxThreshLabelText() + "," + UserVariables.getChan1MaxThresh());
        paramStream.println(UserInterface.getChan2MaxThreshLabelText() + "," + UserVariables.getChan2MaxThresh());
        paramStream.println(UserInterface.getCurveFitTolLabelText() + "," + UserVariables.getCurveFitTol());
        paramStream.println(UserInterface.getC2CurveFitTolLabelText() + "," + UserVariables.getC2CurveFitTol());
        paramStream.println(UserInterface.getColocalToggleText() + "," + UserVariables.isColocal());
        paramStream.println(UserInterface.getPreprocessToggleText() + "," + UserVariables.isPreProcess());
        paramStream.println(UserInterface.getGpuToggleText() + "," + UserVariables.isGpu());
        paramStream.println(UserInterface.getPrevResToggleText() + "," + UserVariables.isPrevRes());
        paramStream.close();
    }

    public static File getDirectory() {
        return outputDir;
    }

    public boolean isGpuEnabled() {
        return gpuEnabled;
    }

//    public double getSigEst() {
//        return sigEst;
//    }
    protected boolean getActiveImages() {
        if (IJ.getInstance() != null) {
            inputs = GenUtils.specifyInputs(labels);
            if (inputs == null) {
                return false;
            }
            stacks[0] = inputs[0].getImageStack();
            if (inputs[1] == null) {
                monoChrome = true;
            } else {
                stacks[1] = inputs[1].getImageStack();
            }
        }
        if (stacks[0].getProcessor(1).getNChannels() > 1 || (!monoChrome && stacks[1].getProcessor(1).getNChannels() > 1)) {
            GenUtils.error("Monochrome images required.");
            return false;
        }
        if (!monoChrome && stacks[0].getSize() != stacks[1].getSize()) {
            GenUtils.error("Stacks must have same number of slices.");
            return false;
        }
        return true;
    }
}
