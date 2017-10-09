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
package ParticleTracking;

import IAClasses.Utils;
import IO.DataWriter;
import MacroWriter.MacroWriter;
import Particle.IsoGaussian;
import Particle.Particle;
import Particle.ParticleArray;
import Revision.Revision;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import ui.DetectionGUI;

/**
 *
 * @author Dave Barry <david.barry at crick.ac.uk>
 */
public class FLAP_ extends GPUAnalyse {

    private final String RESULTS_HEADINGS = String.format("X\tY\tFrame\tChannel 1\tChannel 2\tChannel 1 %c\tChannel 2 %c", '\u03C3', '\u03C3');

    public FLAP_() {
        super();
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
            inputDir = buildStacks(false);
            if (inputDir == null) {
                return;
            }
        }
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

    public void analyse(File inputDir) throws Exception {
        ImageStack[] stacks = getStacks();
        File outputDir = Utilities.getFolder(inputDir, "Specify directory for output files...", true);
        String parentDir = GenUtils.openResultsDirectory(outputDir + delimiter + title);
        if (stacks == null) {
            return;
        }
        startTime = System.currentTimeMillis();
        trajectories = new ArrayList();
        int nC2 = stacks[1].getSize();
        ParticleArray c1Particles = findC1Particles();
        ParticleArray particles = duplicateDetections(c1Particles, nC2);
        detectParticles(0, 0, particles, false, true, UserVariables.getCurveFitTol(), stacks[1]);
        TrajectoryBuilder.updateTrajectories(particles, UserVariables.getTimeRes(), UserVariables.getTrajMaxStep(),
                UserVariables.getSpatialRes(), false, Utils.getStackMinMax(inputs[0].getImageStack())[1], trajectories, false);
        TextWindow results = new TextWindow(title + " Results", RESULTS_HEADINGS,
                new String(), 1000, 500);
        int n = trajectories.size();
        for (int i = 0; i < n; i++) {
            ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
            if (traj != null) {
                traj.printTrajectory(i + 1, results, numFormat, title);
            }
        }
        if (trajectories.size() > 0) {
            ImageStack maps = mapTrajectories((new RGBStackMerge()).mergeStacks(stacks[1].getWidth(), stacks[1].getHeight(), stacks[1].getSize(), null, stacks[1], null, true),
                    trajectories, UserVariables.getSpatialRes(), UserVariables.getMinTrajLength(),
                    UserVariables.getTimeRes(), true, 0, trajectories.size() - 1, 1, false, calcParticleRadius(UserVariables.getSpatialRes(), UserVariables.getSigEstRed()));
            results.append("\nAnalysis Time (s): " + numFormat.format((System.currentTimeMillis() - startTime) / 1000.0));
            results.setVisible(true);
            DataWriter.saveTextWindow(results, new File(String.format("%s%s%s", parentDir, File.separator, "results.csv")), RESULTS_HEADINGS);
            String[] colHeadings = new String[trajectories.size()];
            for (int t = 0; t < colHeadings.length; t++) {
                colHeadings[t] = String.format("Trajectory %d", t);
            }
            DataWriter.saveValues(extractFluorVals(trajectories, stacks[1].size()), new File(String.format("%s%s%s", parentDir, File.separator, "fluorVals.csv")), colHeadings, null);
            try {
                printTrajectories(trajectories, new File(String.format("%s%s%s", parentDir, File.separator, "AllParticleData.csv")), stacks[1].size());
            } catch (IOException e) {
            }
            if (maps != null) {
                (new ImagePlus("Trajectory Maps", maps)).show();
                IJ.saveAs((new ImagePlus("", maps)), "TIF", parentDir + "/trajectories.tif");
            }
        } else {
            IJ.log("No Particle Trajectories Constructed.");
        }
        printParams(parentDir);
    }

    protected ParticleArray findC1Particles() {
        ImageStack[] stacks = getStacks();
        if (UserVariables.isGpu()) {
            return cudaFindParticles(false, 0, 0, null);
        } else {
            return findParticles(false, 0, 0, UserVariables.getCurveFitTol(), stacks[0], null, false, false, UserVariables.isFitC2());
        }
    }

    ParticleArray duplicateDetections(ParticleArray particles, int nDuplicates) {
        ParticleArray output = new ParticleArray(nDuplicates);
        ArrayList<Particle> slice = particles.getLevel(0);
        for (int i = 0; i < nDuplicates; i++) {
            for (Particle p : slice) {
                if (p instanceof IsoGaussian) {
                    output.addDetection(i, new IsoGaussian(i, (IsoGaussian) p));
                }
            }
        }
        return output;
    }

    public void detectParticles(int startSlice, int i, ParticleArray particles, boolean floatingSigma, boolean fitC2Gaussian, double fitTol, ImageStack stack) {
        int fitRad = calcParticleRadius(UserVariables.getSpatialRes(), UserVariables.getSigEstRed());
        int pSize = 2 * fitRad + 1;
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
//        double c2Threshold = stack == null ? 0.0 : Utils.getPercentileThresh(stack, UserVariables.getChan2MaxThresh());
        int depth = stack.getSize();
        for (int d = 0; d < depth; d++) {
            ArrayList<Particle> slice = particles.getLevel(d);
            for (Particle p : slice) {
                int x0 = (int) Math.round(p.getX() / UserVariables.getSpatialRes());
                int y0 = (int) Math.round(p.getY() / UserVariables.getSpatialRes());
                ImageProcessor ip = stack.getProcessor(d + 1);
//                if (stack != null && ip.getPixelValue(x0, y0) > c2Threshold) {
                Particle p2 = new Particle(i - startSlice, x0 * UserVariables.getSpatialRes(), y0 * UserVariables.getSpatialRes(), ip.getPixelValue(x0, y0));
                if (fitC2Gaussian) {
                    Utils.extractValues(xCoords, yCoords, pixValues, x0, y0, ip);
                    ArrayList<IsoGaussian> c2Fits = doFitting(xCoords, yCoords, pixValues,
                            floatingSigma, x0, y0, fitRad, UserVariables.getSpatialRes(),
                            i - startSlice, UserVariables.getSigEstGreen() / UserVariables.getSpatialRes());
                    IsoGaussian c2 = c2Fits.get(0);
                    if (c2.getFit() >= fitTol) {
                        p2 = c2;
                    }
                }
                p.setColocalisedParticle(p2);
//                }
            }
        }
    }

    public boolean showDialog() {
        DetectionGUI ui = new DetectionGUI(null, true, title, this, true);
        ui.setVisible(true);
        return ui.isWasOKed();
    }

    double[][] extractFluorVals(ArrayList<ParticleTrajectory> trajectories, int nFrames) {
        int nTraj = trajectories.size();
        double[][] output = new double[nFrames][nTraj];
        for (int t = 0; t < nTraj; t++) {
            ParticleTrajectory traj = trajectories.get(t);
            Particle p = traj.getEnd();
            int f = nFrames - 1;
            while (p != null) {
                Particle p2 = p.getColocalisedParticle();
                double mag = p2 != null ? p2.getMagnitude() : 0.0;
                output[f--][t] = mag;
                p = p.getLink();
            }
        }
        return output;
    }
}
