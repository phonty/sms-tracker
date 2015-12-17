/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ChannelSplitter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Rectangle;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author barry05
 */
public class VolumeAnalysis extends Analyse_ {

    private int outputsize = 51, midpoint = (outputsize - 1) / 2;
    double spatialRes = UserVariables.getSpatialRes();
    double timeRes = UserVariables.getTimeRes();
    double chan1MaxThresh = UserVariables.getChan1MaxThresh();
    double minTrajLength = UserVariables.getMinTrajLength();
    double c1CurveFitTol = UserVariables.getC1CurveFitTol();
    boolean colocal = UserVariables.isColocal();
    private int xyPartRad;
    ImagePlus imp;

//    public static void main(String args[]) {
//        File image = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\Tracking Test Sequences"), null);
//        ImageStack stack = Utils.buildStack(image);
//        ImagePlus imp = new ImagePlus("Stack", stack);
//        Volume_Analysis instance = new Volume_Analysis(imp);
//        if (instance.showDialog()) {
//            instance.analyse();
//        }
//        return;
//    }
    public VolumeAnalysis() {
        super();
    }

    public VolumeAnalysis(ImagePlus imp) {
        super(imp, null);
        this.imp = imp;
        ImagePlus tempImps[] = new ImagePlus[3];
        tempImps = ChannelSplitter.split(imp);
        stacks[0] = tempImps[0].getImageStack();
        stacks[1] = tempImps[1].getImageStack();
        this.stacks[0] = imp.getImageStack();
    }

    public void analyse(File inputDir) {
        if (stacks[0] != null) {
            IJ.register(this.getClass());
            int i, count;
            int width = stacks[0].getWidth(), height = stacks[0].getHeight();

            findParticles(1.0, true, 0, stacks[0].getSize() - 1, UserVariables.getC1CurveFitTol(), stacks[0], stacks[1] == null);

            TextWindow results = new TextWindow(title + " Results", "X\tY\tFrame\tChannel 1 ("
                    + UserVariables.channels[UserVariables.getC1Index()]
                    + ")\tChannel 2 (" + UserVariables.channels[UserVariables.getC2Index()]
                    + ")\tChannel 2 " + '\u03C3' + "x\tChannel 2 " + '\u03C3' + "y\t" + '\u03B8',
                    new String(), 1000, 500);
            results.append(imp.getTitle() + "\n\n");
            TextWindow resultSummary = new TextWindow(title + " Results Summary",
                    "Particle\tType\t% Colocalisation\tDuration (s)\tDisplacement (" + IJ.micronSymbol
                    + "m)\tVelocity (" + IJ.micronSymbol + "m/s)\tDirectionality\tDiffusion Coefficient ("
                    + IJ.micronSymbol + "m^2/s)" + "\tFractal Dimension"
                    + "\tFluorescence Ratio ("
                    + UserVariables.channels[UserVariables.getC2Index()] + "/"
                    + UserVariables.channels[UserVariables.getC1Index()]
                    + ")\tAngle Spread\tStep Spread\tDC\tCurvature\tC2 Fluor Area\tC2 Fluor Skew",
                    new String(), 1200, 500);
            resultSummary.append(imp.getTitle() + "\n\n");

            int n = trajectories.size();
            for (i = 0; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (!(traj.getSize() > minTrajLength && ((traj.getType(0.1) == ParticleTrajectory.COLOCAL)
                        || ((traj.getType(0.1) == ParticleTrajectory.NON_COLOCAL) && !colocal)))) {
                    trajectories.remove(i);
                    i--;
                    n--;
                }
            }
            n = trajectories.size();
            mapTrajectories(stacks[0], trajectories, spatialRes, minTrajLength, timeRes, true, 0, trajectories.size() - 1, 1, false, calcParticleRadius(spatialRes, SIGMAS[UserVariables.getC1Index()]));
            ArrayList distributions = new ArrayList();
            xyPartRad = calcParticleRadius(spatialRes, SIG_EST_RED);
            int cropRad = 4 * xyPartRad + 1;
            for (i = 0, count = 1; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                int s = traj.getSize();
                int t = traj.getType(0.1);
                if (s > minTrajLength && ((t == ParticleTrajectory.COLOCAL)
                        || ((t == ParticleTrajectory.NON_COLOCAL) && !colocal))) {
                    Particle current = traj.getEnd();
                    double xsum = 0.0, ysum = 0.0;
                    int peakTime = (int) Math.round(traj.getPeakTime() / timeRes);
                    while (current != null) {
                        xsum += current.getX() / spatialRes;
                        ysum += current.getY() / spatialRes;
                        current = current.getLink();
                    }
                    int x = (int) Math.round(xsum / s);
                    int y = (int) Math.round(ysum / s);
                    if (!((x < cropRad) || (y < cropRad) || (stacks[0].getWidth() - x < cropRad)
                            || (stacks[0].getHeight() - y < cropRad))) {
                        ImageStack dist = new ImageStack(cropRad, cropRad);
                        Rectangle roi = new Rectangle(x - 2 * xyPartRad, y - 2 * xyPartRad, cropRad, cropRad);
                        for (int k = peakTime + 1; k <= stacks[0].getSize() && k - peakTime < midpoint + 1; k++) {
                            ImageProcessor currentIP = stacks[0].getProcessor(k);
                            currentIP.setRoi(roi);
                            dist.addSlice("" + count, currentIP.crop());
                        }
                        while (dist.getSize() < midpoint + 1) {
                            dist.addSlice("" + count, new ColorProcessor(cropRad, cropRad));
                        }
                        for (int k = peakTime; k > 0 && peakTime - k < midpoint - 1; k--) {
                            ImageProcessor currentIP = stacks[0].getProcessor(k);
                            currentIP.setRoi(roi);
                            dist.addSlice("" + count, currentIP.crop(), 0);
                        }
                        while (dist.getSize() < outputsize) {
                            dist.addSlice("" + count, new ColorProcessor(cropRad, cropRad), 0);
                        }
                        distributions.add(dist);
                        count++;
//                        if (intensPlot) {
//                            plotIntensity(i, count);
//                        }
//                        if (trajPlot) {
//                            plotTrajectory(width, height, i, count);
//                        }
                        printData(i, resultSummary, count);
                        traj.printTrajectory(count, results, numFormat, title);
                    }
                }
            }
            double distVals[][] = new double[outputsize][cropRad * cropRad];
            for (int d = 0; d < outputsize; d++) {
                Arrays.fill(distVals[d], 0.0);
            }
            double maxval = -Double.MAX_VALUE;
            for (int m = 0; m < distributions.size(); m++) {
                ImageStack thisStack = (ImageStack) distributions.get(m);
                for (int c = 0; c < thisStack.getSize(); c++) {
                    ImageProcessor slice = (ColorProcessor) thisStack.getProcessor(c + 1);
                    for (int x = 0; x < cropRad; x++) {
                        for (int y = 0; y < cropRad; y++) {
                            distVals[c][x + y * cropRad] += (slice.getPixel(x, y, null))[1];
                            if (distVals[c][x + y * cropRad] > maxval) {
                                maxval = distVals[c][x + y * cropRad];
                            }
                        }
                    }
                }
            }
            ImageStack output = new ImageStack(cropRad, cropRad);
            for (int b = 0; b < distVals.length; b++) {
                output.addSlice("" + b, new FloatProcessor(cropRad, cropRad, distVals[b]));
            }
            ImagePlus outimp = new ImagePlus("Sum of Distributions", output);
            outimp.setDisplayRange(0.0, maxval);
            outimp.show();
            results.append(toString());
            results.setVisible(true);
            resultSummary.setVisible(true);
        }
        return;
    }

    public ParticleArray findParticles(double searchScale, boolean update, int startSlice, int endSlice, double fitTol, ImageStack stack, boolean monoChrome) {
        if (stack == null) {
            return null;
        }
//        SIG_EST_RED = (0.21 * 650.0 / 1.4) / (spatialRes * 1000.0);
        xyPartRad = (int) Math.round(2.0 * SIG_EST_RED / 0.95);
        int i, noOfImages = stack.getSize(), width = stack.getWidth(), height = stack.getHeight(),
                size = width * height, arraySize = endSlice - startSlice + 1;
        byte c1Pix[] = new byte[size], c2Pix[] = new byte[size],
                c3Pix[] = new byte[size];
        int c1X, c1Y, pSize = 2 * xyPartRad + 1;
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        ParticleArray particles = new ParticleArray(arraySize);
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
            IJ.freeMemory();
            IJ.showStatus("Analysing Frame " + i);
            IJ.showProgress(i, noOfImages);
            if (!monoChrome) {
                byte[][] tempPix = new byte[3][size];
                ((ColorProcessor) stack.getProcessor(i + 1)).getRGB(tempPix[0], tempPix[1], tempPix[2]);
                c1Pix = tempPix[UserVariables.getC1Index()];
                c2Pix = tempPix[UserVariables.getC2Index()];
            } else {
                c1Pix = (byte[]) (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte().getPixels());
                c2Pix = null;
            }
            FloatProcessor chan1Proc = (FloatProcessor) preProcess(new ByteProcessor(width, height, c1Pix, null), SIG_EST_RED);
            ByteProcessor thisC1Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, UserVariables.FOREGROUND, chan1Proc, chan1MaxThresh, true);
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (thisC1Max.getPixel(c1X, c1Y) == UserVariables.FOREGROUND) {
                        IsoGaussian c1Gaussian = null;
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
                        /*
                         * Remove adjacent Gaussians
                         */
                        IsoGaussianFitter c1GF = new IsoGaussianFitter(xCoords, yCoords, pixValues, false);
                        c1GF.doFit(SIG_EST_RED);
                        //if (c1GF.getXsig() < (c1SigmaTol * xySigEst)) {
                        if (c1GF.getRSquared() > c1CurveFitTol) {
                            c1Gaussian = new IsoGaussian((c1GF.getX0() + c1X - xyPartRad) * spatialRes,
                                    (c1GF.getY0() + c1Y - xyPartRad) * spatialRes, c1GF.getMag(),
                                    c1GF.getXsig(), c1GF.getYsig(), c1GF.getRSquared() - c1CurveFitTol);
                        } else {
                            c1Gaussian = new IsoGaussian(c1X * spatialRes, c1Y * spatialRes, chan1Proc.getPixelValue(c1X, c1Y),
                                    SIG_EST_RED, SIG_EST_RED, c1GF.getRSquared() - c1CurveFitTol);
                        }
                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        if (c1Gaussian != null) {
                            particles.addDetection(i - startSlice, new Particle(i - startSlice, c1Gaussian, null, null, -1));
                        }
                        //}
                    }
                }
            }
        }
        if (update) {
            TrajectoryBuilder.updateTrajectories(particles, timeRes, UserVariables.getTrajMaxStep(), spatialRes, true, 1.0, trajectories);
        }
        return particles;
    }

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
        traj.smooth();
        traj.calcMSD(-1);
        traj.calcAngleSpread();
        traj.calcStepSpread();
        traj.calcDirectionality(traj.getPoints()[0], traj.getPoints()[1]);
        traj.calcFluorSpread();
        double fluorMaj = Math.max(traj.getxFluorSpread(), traj.getyFluorSpread());
        double fluorMin = Math.min(traj.getxFluorSpread(), traj.getyFluorSpread());
        double fluorArea = Math.PI * 4.0 * fluorMin * fluorMaj * spatialRes * spatialRes;
        double displacement = traj.getDisplacement(traj.getEnd(), traj.getSize());
        double duration = traj.getSize() * timeRes;
        int type = traj.getType(0.1);
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
                + decFormat.format(traj.getMeanKappa()) + "\t"
                + decFormat.format(fluorArea) + "\t"
                + decFormat.format(fluorMin / fluorMaj));
        return true;
    }
}
