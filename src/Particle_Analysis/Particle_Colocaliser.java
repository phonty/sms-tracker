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
import IO.PropertyWriter;
import Particle.Particle;
import Particle.ParticleArray;
import ParticleTracking.GPUAnalyse;
import ParticleTracking.UserVariables;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.io.File;
import java.util.ArrayList;
import java.util.Properties;
import ui.DetectionGUI;

public class Particle_Colocaliser extends GPUAnalyse {

    private String resultsHeadings = String.format("Image\tChannel 1 Detections\tColocalised Channel 2 Detections\t%% Colocalisation\t%c (nm)", '\u0394'),
            coordHeadings = "C0_X\tC0_Y\tC1_X\tC1_Y";
    private Properties props;

    public Particle_Colocaliser() {
        super();
        title = "Particle Colocaliser";
        gpuEnabled = false;
    }

    public ParticleArray findParticles(boolean update, int startSlice, int endSlice, double c1FitTol, ImageStack channel1, ImageStack channel2) {
        return findParticles(update, startSlice, endSlice, c1FitTol,
                channel1, channel2, true, false, UserVariables.isFitC2());
    }

    protected ParticleArray findParticles() {
        ImageStack[] stacks = getStacks();
        return findParticles(false, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1]);
    }

    public void analyse(File inputDir) throws Exception {
        ImageStack[] stacks = getStacks();
        File outputDir = Utilities.getFolder(inputDir, "Specify directory for output files...", true);
        if (stacks != null) {
            startTime = System.currentTimeMillis();
            buildOutput(findParticles(), outputDir);
        }
    }

    void buildOutput(ParticleArray curves, File outputDir) {
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
            ArrayList<Particle> detections = curves.getLevel(i);
            for (int j = 0; j < detections.size(); j++) {
                Particle p1 = detections.get(j);
                String coordString = String.format("%3.3f\t%3.3f", p1.getX(), p1.getY());
                Utils.drawParticle(ch1proc, p1, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), false);
                count++;
                Particle p2 = p1.getColocalisedParticle();
                if (p2 != null) {
                    Utils.drawParticle(ch2proc, p2, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), false);
                    colocalisation++;
                    sepsum += Utils.calcDistance(p1.getX(), p1.getY(), p2.getX(), p2.getY());
                    coordString = String.format("%s\t%3.3f\t%3.3f", coordString, p2.getX(), p2.getY());
                }
                particleCoords.append(coordString);
            }
            results.append(String.format("Slice %d\t%d\t%d\t%3.3f\t%3.3f", i, count, colocalisation, (100.0 * colocalisation / count), (1000.0 * sepsum / count)));
            outStack[0].addSlice("" + i, ch1proc);
            outStack[1].addSlice("" + i, ch2proc);
        }
        progress.dispose();
        results.setVisible(true);
        particleCoords.setVisible(true);
        try {
            String resultsDir = GenUtils.openResultsDirectory(String.format("%s%s%s", outputDir, File.separator, title));
//            String c1Title = String.format("%s_Detections.tif", inputs[0].getTitle());
            String c1Title = String.format("%s_Detections.tif", "C1");
            IJ.save(new ImagePlus(c1Title, outStack[0]), String.format("%s%s%s", resultsDir, File.separator, c1Title));
//            String c2Title = String.format("%s_Detections.tif", inputs[1].getTitle());
            String c2Title = String.format("%s_Detections.tif", "C2");
            IJ.save(new ImagePlus(c2Title, outStack[1]), String.format("%s%s%s", resultsDir, File.separator, c2Title));
            DataWriter.saveTextWindow(results, new File(String.format("%s%s%s", resultsDir, File.separator, "Results.csv")), resultsHeadings);
            PropertyWriter.printProperties(props, resultsDir, title, true);
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

    public boolean showDialog() {
        DetectionGUI ui = new DetectionGUI(null, true, title, this, false);
        ui.setVisible(true);
        props = ui.getProperties();
        return ui.isWasOKed();
    }

}
