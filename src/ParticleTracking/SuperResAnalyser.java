/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.ProgressDialog;
import ParticleTracking.ParticleArray;
import ParticleTracking.UserVariables;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Toolkit;
import java.text.DecimalFormat;
import java.util.ArrayList;

/**
 *
 * @author barry05
 */
public class SuperResAnalyser extends Colocalisation_Analysis {

    private int scaleFactor = 10;
    private final String TITLE = this.getClass().getName();

//    public static void main(String args[]) {
//        (new SuperResAnalyser(new ImagePlus("C:\\Users\\barry05\\Desktop\\SuperResTest.tif"))).run(null);
//    }
//    public SuperResAnalyser(ImagePlus imp) {
//        super(imp);
//        imp.show();
//        this.imp = imp;
//        if (imp != null) {
//            this.stacks = imp.getStack();
//        } else {
//            stacks = null;
//        }
//    }

    public void run(String arg) {
        if (IJ.getInstance() != null) {
            imp = IJ.getImage();
//            stacks = imp.getImageStack();
        }
        if (imp == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image stack open.");
            return;
        }
        if (showDialog()) {
            resultsHeadings = "Image\tChannel 1 (" + UserVariables.channels[UserVariables.getC1Index()]
                    + ") Detections\tColocalised Channel 2 (" + UserVariables.channels[UserVariables.getC2Index()]
                    + ") Detections\t% Colocalisation";
            UserVariables.setPreProcess(true);
            Analyse_ analyser = new Analyse_(stacks);
            analyser.calcParticleRadius(UserVariables.getSpatialRes(), SIG_EST_RED);
            //Timelapse_Analysis.setGaussianRadius(0.139 / Timelapse_Analysis.getSpatialRes());
            //IJ.saveAs(buildOutput(analyser), "TIF", "C:\\Users\\barry05\\Desktop\\SuperResTestOutputII.tif");
            buildOutput(analyser);
        }
    }

    public boolean draw2DGaussian(ImageProcessor image, IsoGaussian g, double tol, double par3) {
        if (image == null || g == null) {
            return false;
        }
        double res = UserVariables.getSpatialRes();
        int x, y, drawRad;
        double x0 = scaleFactor * g.getX() / res;
        double y0 = scaleFactor * g.getY() / res;
        double xSigma = g.getXSigma();
        double value;
        drawRad = (int) Math.round(xSigma * 3.0);
        double fit = g.getFit();
        if (fit < 0.0) {
            fit = 0.0;
        }
        if (g.getMagnitude() > 0.0 && g.getMagnitude() < 255.0) {
            for (x = (int) Math.floor(x0 - drawRad); x <= x0 + drawRad; x++) {
                for (y = (int) Math.floor(y0 - drawRad); y <= y0 + drawRad; y++) {
                    /*
                     * The current pixel value is added so as not to "overwrite"
                     * other Gaussians in close proximity:
                     */
                    value = g.getMagnitude() * Math.exp(-(((x - x0) * (x - x0))
                            + ((y - y0) * (y - y0))) / (2 * xSigma * xSigma));
                    value += image.getPixelValue(x, y);
                    image.putPixelValue(x, y, value);
                }
            }
        } else {
            return false;
        }
        return true;
    }

    ImagePlus buildOutput(Analyse_ analyser) {
        if (stacks == null) {
            return null;
        }
        if (analyser == null) {
            analyser = new Analyse_(stacks);
        }
        DecimalFormat format = new DecimalFormat("000");
        int width = imp.getWidth() * scaleFactor, height = imp.getHeight() * scaleFactor;
        //ImageStack outStack = new ImageStack(width, height);
        ProgressDialog dialog = new ProgressDialog(null, "Processing...", false, TITLE, false);
        dialog.setVisible(true);
        FloatProcessor ch1proc = new FloatProcessor(width, height);
        for (int i = 0; i < stacks[0].getSize(); i++) {
            dialog.updateProgress(i, stacks[0].getSize());
            ParticleArray curves = analyser.findParticles(coFactor, false, i, i, UserVariables.getCurveFitTol(), stacks[0], stacks[1], true, SIG_EST_RED, SIG_EST_GREEN, UserVariables.isColocal(), true, true, UserVariables.getC2CurveFitTol());
            //ImagePlus temp = new ImagePlus("", ch1proc);
            //temp.show();
            //temp.setDisplayRange(0.0, 255.0);
            ArrayList detections = curves.getLevel(0);
            for (int j = 0; j < detections.size(); j++) {
                IsoGaussian c1 = ((IsoGaussian[]) detections.get(j))[0];
                if (draw2DGaussian(ch1proc, c1, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes())) {
                }
                //temp.updateAndDraw();
            }
            IJ.saveAs(new ImagePlus("", ch1proc.duplicate()), "TIF", "C:\\Users\\barry05\\Desktop\\SuperResTest\\Output_" + format.format(i) + ".tif");
            //outStack.addSlice("" + i, ch1proc.duplicate());
        }
        dialog.dispose();
        if (results != null) {
            results.append("\n" + toString());
            results.setVisible(true);
        }
        //ImagePlus output = new ImagePlus("Detected Particles", outStack);
        //output.setDisplayRange(0.0, 255.0);
        return null;
    }
}
