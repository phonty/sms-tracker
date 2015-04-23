/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import ParticleTracking.FloatingMultiGaussFitter;
import ParticleTracking.Particle;
import ParticleTracking.ParticleArray;
import ParticleTracking.UserVariables;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Toolkit;
import java.io.File;
import java.util.ArrayList;

/**
 *
 * @author barry05
 */
public class PSFEstimator extends Analyse_ {

    private TextWindow results;
    private String psfTitle = "PSF Estimator v1.0";

//    public static void main(String args[]) {
//        File image = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\Test Data Sets\\PSF Estimator\\Test 1"), null, true);
////        ImageStack stack = Utils.buildStack(image);
//        ImagePlus imp = Utils.buildStack(image);
//        PSFEstimator instance = new PSFEstimator(imp);
//        if (instance.showDialog()) {
//            instance.analyse();
//        }
//    }

    public PSFEstimator() {
        super();
    }

    public PSFEstimator(ImagePlus imp) {
        super(imp, null);
        this.imp = imp;
    }

    public void analyse() {
        stacks = new ImageStack[2];
        stacks[0] = imp.getImageStack();
        if (stacks[0] != null) {
            calcParticleRadius(UserVariables.getSpatialRes(), SIG_EST_RED);
            IJ.register(this.getClass());
            results = new TextWindow(psfTitle + " Results", "frame\tx (" + IJ.micronSymbol + "m)\ty (" + IJ.micronSymbol + "m)\tA\tsigma (nm)\tR^2",
                    new String(), 1000, 500);
            results.append(imp.getTitle() + "\n\n");
            UserVariables.setnMax(1);
            ParticleArray particles = findParticles(1.0, true, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], monoChrome);
            for (int i = 0; i < particles.getDepth(); i++) {
                ArrayList detections = particles.getLevel(i);
                for (int j = 0; j < detections.size(); j++) {
                    IsoGaussian c1 = ((Particle) detections.get(j)).getC1Gaussian();
                    results.append(String.valueOf(i) + "\t" + String.valueOf(c1.getX() * UserVariables.getSpatialRes())
                            + "\t" + String.valueOf(c1.getY() * UserVariables.getSpatialRes()) + "\t" + String.valueOf(c1.getMagnitude())
                            + "\t" + String.valueOf(c1.getXSigma() * UserVariables.getSpatialRes() * 1000.0) + "\t" + String.valueOf(c1.getFit()));
                }
            }
            results.setVisible(true);
        }
    }

    public boolean showDialog() {
        if (imp == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image stack open.");
            return false;
        }
        GenericDialog gd = new GenericDialog(psfTitle);
        gd.addNumericField("Spatial Resolution", UserVariables.getSpatialRes() * 1000.0, 5, 5, "nm");
        gd.addNumericField("Peak Threshold", UserVariables.getChan1MaxThresh(), 5);
        gd.addNumericField("Fit Tolerance", UserVariables.getCurveFitTol(), 5);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        UserVariables.setSpatialRes(gd.getNextNumber() / 1000.0);
        UserVariables.setChan1MaxThresh(gd.getNextNumber());
        UserVariables.setCurveFitTol(gd.getNextNumber());
        return true;
    }

    public ParticleArray findParticles(double searchScale, boolean update, int startSlice, int endSlice, double fitTol, ImageStack stack, boolean monoChrome) {
        if (stack == null) {
            return null;
        }
        int i, noOfImages = stack.getSize(), width = stack.getWidth(), height = stack.getHeight(),
                arraySize = endSlice - startSlice + 1;
        byte c1Pix[];
        int xyPartRad = calcParticleRadius(UserVariables.getSpatialRes(), SIG_EST_RED);
        int fitRad = (int) Math.ceil(xyPartRad * 4.0 / 3.0);
        int c1X, c1Y, pSize = 2 * fitRad + 1;
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
        ProgressDialog progress = new ProgressDialog(null, "Finding Particles...", false, title, false);
        progress.setVisible(true);
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
            IJ.freeMemory();
            progress.updateProgress(i - startSlice, arraySize);
            c1Pix = (byte[]) (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte().getPixels());
            FloatProcessor chan1Proc = preProcess(new ByteProcessor(width, height, c1Pix, null));
            double c1Threshold = Utils.getPercentileThresh(chan1Proc, UserVariables.getChan1MaxThresh());
            ByteProcessor thisC1Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, UserVariables.FOREGROUND, chan1Proc, c1Threshold, true);
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (thisC1Max.getPixel(c1X, c1Y) == UserVariables.FOREGROUND) {
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
                        FloatingMultiGaussFitter c1Fitter = new FloatingMultiGaussFitter(UserVariables.getnMax(), fitRad, pSize);
                        c1Fitter.fit(pixValues, SIG_EST_RED/UserVariables.getSpatialRes());
                        ArrayList<IsoGaussian> c1Fits = c1Fitter.getFits(spatialRes, c1X - fitRad, c1Y - fitRad, c1Threshold, fitTol);
                        if (c1Fits != null) {
                            for (IsoGaussian c1Fit : c1Fits) {
                                particles.addDetection(i - startSlice, new Particle(i - startSlice, c1Fit, null, null, -1));
                            }
                        }
                    }
                }
            }
        }
        progress.dispose();
        return particles;
    }
}
