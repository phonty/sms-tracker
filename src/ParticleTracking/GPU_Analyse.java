/*
 * Copyright (C) 2015 David Barry <david.barry at cancer.org.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import UtilClasses.GenUtils;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import java.io.File;
import java.util.ArrayList;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class GPU_Analyse extends Analyse_ {

    private final String delimiter = GenUtils.getDelimiter();
    private final int CUDA_FILE_COLS = 5;

    static {
        System.loadLibrary("cuda_gauss_tracker"); // Load native library at runtime
    }

    private native boolean cudaGaussFitter(String folder, String ext, float spatialRes, float sigmaEst, float maxthresh, float fitTol, int startSlice, int endSlice);

//    public static void main(String args[]) {
//        GPU_Analyse instance = new GPU_Analyse();
//        instance.run(null);
//        System.exit(0);
//    }

    public GPU_Analyse() {
        super();
        gpuEnabled = true;
    }

    protected ParticleArray findParticles() {
        if (UserVariables.isGpu()) {
            return cudaFindParticles(SEARCH_SCALE, true, 0, stacks[0].getSize() - 1, stacks[1]);
        } else {
            return findParticles(SEARCH_SCALE, true, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1], false, sigmas[UserVariables.getC1Index()], sigmas[1 - UserVariables.getC1Index()], UserVariables.isColocal());
        }
    }

    public ParticleArray cudaFindParticles(boolean update, int startSlice, int endSlice, ImageStack channel2) {
        return cudaFindParticles(SEARCH_SCALE, update, startSlice, endSlice, channel2);
    }

    public ParticleArray cudaFindParticles(double searchScale, boolean update, int startSlice, int endSlice, ImageStack channel2) {
        if (!cudaGaussFitter(c0Dir.getAbsolutePath(), ext, (float) UserVariables.getSpatialRes() * 1000.0f,
                (float) (sigmas[UserVariables.getC1Index()] / UserVariables.getSpatialRes()), (float) UserVariables.getChan1MaxThresh(),
                (float) UserVariables.getCurveFitTol(), startSlice, endSlice)) {
            IJ.log("CUDA Error");
            return null;
        }
        File cudaFile = new File(c0Dir + delimiter + "cudadata.txt");
        File fileList[] = {cudaFile};
        ArrayList<double[]>[] cudaData = GenUtils.readData(CUDA_FILE_COLS, fileList, delimiter);
        int arraySize = endSlice - startSlice + 1;
        int xyPartRad = calcParticleRadius(UserVariables.getSpatialRes(), sigmas[UserVariables.getC1Index()]);
        int fitRad = (int) Math.ceil(xyPartRad);
        int c2Points[][];
        int pSize = 2 * fitRad + 1;
        int radius = (int) Math.round(fitRad * searchScale);
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
        double c2Thresholds[] = new double[channel2.getSize()];
        ImageStack procChannel2 = new ImageStack(channel2.getWidth(), channel2.getHeight());
        for (int i = 0; i < channel2.getSize(); i++) {
            procChannel2.addSlice(Utils.normalise(preProcess(channel2.getProcessor(i + 1).duplicate()), 1.0));
            c2Thresholds[i] = Utils.getPercentileThresh(procChannel2.getProcessor(i + 1),
                    UserVariables.getChan2MaxThresh());
        }
        for (int f = 0; f < fileList.length; f++) {
            ProgressDialog progress = new ProgressDialog(null,
                    "Reading data for file " + f + " of " + fileList.length + "...",
                    false, title, false);
            progress.setVisible(true);
            int size = cudaData[f].size();
            for (int row = 0; row < size; row++) {
                progress.updateProgress(row, size);
                int t = (int) Math.round(cudaData[f].get(row)[0]);
                double x = cudaData[f].get(row)[1];
                double y = cudaData[f].get(row)[2];
                double mag = cudaData[f].get(row)[3];
                double fit = cudaData[f].get(row)[4];
                IsoGaussian c1Gaussian = new IsoGaussian(x, y, mag, sigmas[UserVariables.getC1Index()], sigmas[UserVariables.getC1Index()], fit);
                int c1X = (int) Math.round(x / UserVariables.getSpatialRes());
                int c1Y = (int) Math.round(y / UserVariables.getSpatialRes());
                IsoGaussian c2Gaussian = null;
                c2Points = Utils.searchNeighbourhood(c1X, c1Y,
                        (int) Math.round(fitRad * searchScale),
                        c2Thresholds[t],
                        procChannel2.getProcessor(t + 1));
                if (c2Points != null) {
                    c2Gaussian = findC2Particle(c1X, c1Y, radius, pSize, c2Thresholds[t],
                            (FloatProcessor) procChannel2.getProcessor(t + 1), true, sigmas[UserVariables.getC2Index()],
                            fitRad, spatialRes);
                }
                if ((c2Gaussian == null && !UserVariables.isColocal())
                        || (c2Gaussian != null && UserVariables.isColocal() && c2Gaussian.getFit() > UserVariables.getC2CurveFitTol())) {
                    particles.addDetection(t - startSlice, new Particle(t, c1Gaussian, c2Gaussian, null, -1));
                }
            }
            progress.dispose();
        }
        if (update) {
            TrajectoryBuilder.updateTrajectories(particles, UserVariables.getTimeRes(), UserVariables.getTrajMaxStep(), spatialRes, true, Utils.getStackMinMax(stacks[0])[1], trajectories);
        }
        return particles;
    }
}
