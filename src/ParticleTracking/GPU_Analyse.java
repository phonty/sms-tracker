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

    public static void main(String args[]) {
        GPU_Analyse instance = new GPU_Analyse();
        instance.run(null);
        System.exit(0);
    }

    public GPU_Analyse() {
        super();
        gpuEnabled = true;
    }

    protected ParticleArray findParticles() {
        if (UserVariables.isGpu()) {
            return GPU_Analyse.this.cudaFindParticles1C(SEARCH_SCALE, true, 0, stacks[0].getSize() - 1, stacks[1]);
        } else {
            return findParticles(SEARCH_SCALE, true, 0, stacks[0].getSize() - 1, UserVariables.getC1CurveFitTol(), stacks[0], stacks[1], false, sigmas[UserVariables.getC1Index()], sigmas[1 - UserVariables.getC1Index()], UserVariables.isColocal(), true, true, UserVariables.getC2CurveFitTol(), false);
        }
    }

    public ParticleArray cudaFindParticles2C(boolean update, int startSlice, int endSlice, ImageStack channel2) {
        return GPU_Analyse.this.cudaFindParticles2C(SEARCH_SCALE, update, startSlice, endSlice, channel2);
    }

    public ParticleArray cudaFindParticles2C(double searchScale, boolean update, int startSlice, int endSlice, ImageStack channel2) {
        if (!cudaGaussFitter(c0Dir.getAbsolutePath(), ext, (float) UserVariables.getSpatialRes() * 1000.0f,
                (float) (sigmas[UserVariables.getC1Index()] / UserVariables.getSpatialRes()), (float) UserVariables.getChan1MaxThresh(),
                (float) UserVariables.getC1CurveFitTol(), startSlice, endSlice)) {
            IJ.log("CUDA Error");
            return null;
        }
        if (!cudaGaussFitter(c1Dir.getAbsolutePath(), ext, (float) UserVariables.getSpatialRes() * 1000.0f,
                (float) (sigmas[UserVariables.getC2Index()] / UserVariables.getSpatialRes()), (float) UserVariables.getChan2MaxThresh(),
                (float) UserVariables.getC2CurveFitTol(), startSlice, endSlice)) {
            IJ.log("CUDA Error");
            return null;
        }
        File cudaC0File = new File(c0Dir + delimiter + "cudadata.txt");
        File c0FileList[] = {cudaC0File};
        ArrayList<double[]>[] c0CudaData = GenUtils.readData(CUDA_FILE_COLS, c0FileList, delimiter);
        File cudaC1File = new File(c1Dir + delimiter + "cudadata.txt");
        File c1FileList[] = {cudaC1File};
        ArrayList<double[]>[] c1CudaData = GenUtils.readData(CUDA_FILE_COLS, c1FileList, delimiter);
        int arraySize = endSlice - startSlice + 1;
        int xyPartRad = calcParticleRadius(UserVariables.getSpatialRes(), sigmas[UserVariables.getC1Index()]);
        int fitRad = (int) Math.ceil(xyPartRad);
        double radius = fitRad * searchScale * UserVariables.getSpatialRes();
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
        for (int f = 0; f < c0FileList.length; f++) {
            ProgressDialog progress = new ProgressDialog(null,
                    "Searching for colocalised particles...",
                    false, title, false);
            progress.setVisible(true);
            int c0Size = c0CudaData[f].size();
            for (int row = 0; row < c0Size; row++) {
                progress.updateProgress(row, c0Size);
                int c0t = (int) Math.round(c0CudaData[f].get(row)[0]);
                double x0 = c0CudaData[f].get(row)[1];
                double y0 = c0CudaData[f].get(row)[2];
                double mag = c0CudaData[f].get(row)[3];
                double fit = c0CudaData[f].get(row)[4];
                IsoGaussian c1Gaussian = new IsoGaussian(x0, y0, mag, sigmas[UserVariables.getC1Index()], sigmas[UserVariables.getC1Index()], fit);
                IsoGaussian c2Gaussian = null;
                double minDist = Double.MAX_VALUE;
                int minIndex = -1;
                int i = 0;
                int c1t = (int) Math.round(c1CudaData[f].get(i)[0]);
                int c1Size = c1CudaData[f].size();
                while (i < c1Size && c1t <= c0t) {
                    c1t = (int) Math.round(c1CudaData[f].get(i)[0]);
                    if (!(c1t < c0t)) {
                        double x1 = c1CudaData[f].get(i)[1];
                        double y1 = c1CudaData[f].get(i)[2];
                        double dist = Math.abs(x1 - x0) + Math.abs(y1 - y0);
                        if (dist < minDist && c0t == c1t) {
                            minDist = dist;
                            minIndex = i;
                        }
                    }
                    i++;
                }
                if (minIndex >= 0 && minDist < radius && c1CudaData[f].get(minIndex)[4] > UserVariables.getC2CurveFitTol()) {
                    c2Gaussian = new IsoGaussian(c1CudaData[f].get(minIndex)[1],
                            c1CudaData[f].get(minIndex)[2],
                            c1CudaData[f].get(minIndex)[3],
                            sigmas[UserVariables.getC2Index()],
                            sigmas[UserVariables.getC2Index()],
                            c1CudaData[f].get(minIndex)[4]);
                    c1CudaData[f].remove(minIndex);
                }
                if (c1Gaussian.getFit() > UserVariables.getC1CurveFitTol()) {
                    particles.addDetection(c0t - startSlice, new Particle(c0t, c1Gaussian, c2Gaussian, null, -1));
                }
            }
            progress.dispose();
        }
        if (update) {
            TrajectoryBuilder.updateTrajectories(particles, UserVariables.getTimeRes(), UserVariables.getTrajMaxStep(), spatialRes, true, Utils.getStackMinMax(stacks[0])[1], trajectories);
        }
        return particles;
    }

    public ParticleArray cudaFindParticles1C(boolean update, int startSlice, int endSlice, ImageStack channel2) {
        return GPU_Analyse.this.cudaFindParticles1C(SEARCH_SCALE, update, startSlice, endSlice, channel2);
    }

    public ParticleArray cudaFindParticles1C(double searchScale, boolean update, int startSlice, int endSlice, ImageStack channel2) {
        if (!cudaGaussFitter(c0Dir.getAbsolutePath(), ext, (float) UserVariables.getSpatialRes() * 1000.0f,
                (float) (sigmas[UserVariables.getC1Index()] / UserVariables.getSpatialRes()), (float) UserVariables.getChan1MaxThresh(),
                (float) UserVariables.getC1CurveFitTol(), startSlice, endSlice)) {
            IJ.log("CUDA Error");
            return null;
        }
        File cudaC0File = new File(c0Dir + delimiter + "cudadata.txt");
        File c0FileList[] = {cudaC0File};
        ArrayList<double[]>[] c0CudaData = GenUtils.readData(CUDA_FILE_COLS, c0FileList, delimiter);
        int arraySize = endSlice - startSlice + 1;
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
        for (int f = 0; f < c0FileList.length; f++) {
            ProgressDialog progress = new ProgressDialog(null,
                    "Searching for colocalised particles...",
                    false, title, false);
            progress.setVisible(true);
            int c0Size = c0CudaData[f].size();
            for (int row = 0; row < c0Size; row++) {
                progress.updateProgress(row, c0Size);
                int c0t = (int) Math.round(c0CudaData[f].get(row)[0]);
                double x0 = c0CudaData[f].get(row)[1];
                double y0 = c0CudaData[f].get(row)[2];
                double mag = c0CudaData[f].get(row)[3];
                double fit = c0CudaData[f].get(row)[4];
                IsoGaussian c1Gaussian = new IsoGaussian(x0, y0, mag, sigmas[UserVariables.getC1Index()], sigmas[UserVariables.getC1Index()], fit);
                if (c1Gaussian.getFit() > UserVariables.getC1CurveFitTol()) {
                    particles.addDetection(c0t - startSlice, new Particle(c0t, c1Gaussian, null, null, -1));
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
