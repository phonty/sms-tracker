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

import IAClasses.ProgressDialog;
import IAClasses.Utils;
import IO.DataWriter;
import Particle.IsoGaussian;
import Particle.ParticleArray;
import ParticleTracking.GPUAnalyse;
import ParticleTracking.UserVariables;
import static ParticleTracking.UserVariables.BLUE;
import static ParticleTracking.UserVariables.GREEN;
import static ParticleTracking.UserVariables.RED;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.io.File;
import java.util.ArrayList;

public class Particle_Colocaliser extends GPUAnalyse {

    private String resultsHeadings = String.format("Image\tChannel 1 Detections\tColocalised Channel 2 Detections\t%% Colocalisation\t%c (nm)", '\u0394'),
            coordHeadings = "C0_X\tC0_Y\tC1_X\tC1_Y\tC0_\u03c3\tC1_\u03c3\tC0 Fit\tC1 Fit";

    public Particle_Colocaliser() {
        super();
        gpuEnabled = false;
    }

    public void colocalise() {

    }

    public ParticleArray findParticles(boolean update, int startSlice, int endSlice, double c1FitTol, ImageStack channel1, ImageStack channel2) {
        return findParticles(update, startSlice, endSlice, c1FitTol,
                channel1, channel2, true, false, true);
    }

    protected ParticleArray findParticles() {
        ImageStack[] stacks = getStacks();
        return findParticles(false, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1]);
    }

    public void analyse(File inputDir) throws Exception {
        ImageStack[] stacks = getStacks();
        if (stacks != null) {
            startTime = System.currentTimeMillis();
            buildOutput(findParticles());
        }
    }

    void buildOutput(ParticleArray curves) {
        ImageStack[] stacks = getStacks();
        int colocalisation, count;
        int width = stacks[0].getWidth(), height = stacks[0].getHeight();
        ImageStack[] outStack = new ImageStack[2];
        outStack[0] = new ImageStack(width, height);
        outStack[1] = new ImageStack(width, height);
        double sepsum;
        TextWindow results = new TextWindow(title + " Results", resultsHeadings, new String(), 1000, 500);
        TextWindow particleCoords = new TextWindow(title + " Particle Coordinates", coordHeadings, new String(), 1000, 500);
        ProgressDialog progress = new ProgressDialog(null, "Analysing Stacks...", false, title, false);
        progress.setVisible(true);
        for (int i = 0; i < stacks[0].getSize(); i++) {
            progress.updateProgress(i, stacks[0].getSize());
            colocalisation = 0;
            count = 0;
            sepsum = 0.0;
            FloatProcessor ch1proc = new FloatProcessor(width, height);
            FloatProcessor ch2proc = new FloatProcessor(width, height);
            ArrayList detections = curves.getLevel(i);
            for (int j = 0; j < detections.size(); j++) {
                IsoGaussian c1 = (IsoGaussian) detections.get(j);
                String coordString;
                if (Utils.draw2DGaussian(ch1proc, c1, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), false)) {
                    count++;
                    IsoGaussian c2 = (IsoGaussian) c1.getColocalisedParticle();
                    if (Utils.draw2DGaussian(ch2proc, c2, UserVariables.getCurveFitTol(),
                            UserVariables.getSpatialRes(), false)) {
                        colocalisation++;
                        sepsum += Utils.calcDistance(c1.getX(), c1.getY(), c2.getX(), c2.getY());
                        coordString = String.valueOf(c1.getX()) + "\t" + String.valueOf(c1.getY())
                                + "\t" + String.valueOf(c2.getX()) + "\t" + String.valueOf(c2.getY())
                                + "\t" + String.valueOf(c1.getXSigma() * UserVariables.getSpatialRes())
                                + "\t" + String.valueOf(c2.getXSigma() * UserVariables.getSpatialRes())
                                + "\t" + String.valueOf(c1.getFit()) + "\t" + String.valueOf(c2.getFit());
                    } else {
                        coordString = String.valueOf(c1.getX()) + "\t" + String.valueOf(c1.getY())
                                + "\t \t \t" + String.valueOf(c1.getXSigma()) + "\t \t"
                                + String.valueOf(c1.getFit());
                    }
                    particleCoords.append(coordString);
                }
            }
            results.append("Slice " + i + "\t" + count + "\t" + colocalisation
                    + "\t" + numFormat.format(100.0 * colocalisation / count)
                    + "\t" + numFormat.format(1000.0 * sepsum / count));

            outStack[0].addSlice("" + i, ch1proc);
            outStack[1].addSlice("" + i, ch2proc);
        }
        progress.dispose();
        results.append("\n" + toString());
        results.setVisible(true);
        particleCoords.setVisible(true);
        try {
            String resultsDir = GenUtils.openResultsDirectory(String.format("%s%s%s", c0Dir.getParent(), File.separator, title));
//            String c1Title = String.format("%s_Detections.tif", inputs[0].getTitle());
            String c1Title = String.format("%s_Detections.tif", "C1");
            IJ.save(new ImagePlus(c1Title, outStack[0]), String.format("%s%s%s", resultsDir, File.separator, c1Title));
//            String c2Title = String.format("%s_Detections.tif", inputs[1].getTitle());
            String c2Title = String.format("%s_Detections.tif", "C2");
            IJ.save(new ImagePlus(c2Title, outStack[1]), String.format("%s%s%s", resultsDir, File.separator, c2Title));
            DataWriter.saveTextWindow(results, new File(String.format("%s%s%s", resultsDir, File.separator, "Results.csv")), resultsHeadings);
        } catch (Exception e) {
            GenUtils.error(String.format("Failed to generate output files: %s", e.toString()));
        }
    }

    byte[] outPix(ImageProcessor ch1, int w, int h) {
        if (ch1 == null) {
            return (byte[]) (new ByteProcessor(w, h)).getPixels();
        }
        return (byte[]) ((new TypeConverter(ch1, true)).convertToByte()).getPixels();
    }

}
