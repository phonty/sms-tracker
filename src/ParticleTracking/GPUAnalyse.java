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

import Particle.ParticleArray;
import Particle_Analysis.Particle_Tracker;
import Particle.IsoGaussian;
import IAClasses.ProgressDialog;
import UtilClasses.GenUtils;
import ij.IJ;
import ij.ImageStack;
import ij.process.ImageProcessor;
import java.io.File;
import java.util.ArrayList;

public class GPUAnalyse extends Particle_Tracker {

    private final int CUDA_FILE_COLS = 5;

    static {
        try {
            System.loadLibrary("ParticleDetector"); // Load native library at runtime
        } catch (Exception | Error e) {
            IJ.log(String.format("No CUDA runtime library detected - GPU compute functionality disabled.\n"));
            gpuEnabled = false;
        }
        gpuEnabled = true;
    }

    private native boolean cudaGaussFitter(String folder, String ext, float spatialRes, float sigmaEst, float maxthresh, float fitTol, int startSlice, int endSlice);

//    public static void main(String args[]) {
//        GPUAnalyse instance = new GPUAnalyse();
//        instance.run(null);
//        System.exit(0);
//    }
    public GPUAnalyse() {
        super();
        gpuEnabled = true;
    }

    protected ParticleArray findParticles() {
        ImageStack[] stacks = getStacks();
        if (UserVariables.isGpu()) {
            return cudaFindParticles(true, 0, stacks[0].getSize() - 1, stacks[1]);
        } else {
            return findParticles(true, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1], true, false, UserVariables.isFitC2());
        }
    }

    public ParticleArray cudaFindParticles(boolean update, int startSlice, int endSlice, ImageStack channel2) {
        if (!cudaGaussFitter(c0Dir.getAbsolutePath(), ext, (float) UserVariables.getSpatialRes() * 1000.0f,
                (float) (UserVariables.getSigEstRed() / UserVariables.getSpatialRes()), (float) UserVariables.getChan1MaxThresh(),
                (float) UserVariables.getCurveFitTol(), startSlice, endSlice)) {
            IJ.log("CUDA Error");
            return null;
        }
        double c2sigma = (UserVariables.getSigEstGreen() / UserVariables.getSpatialRes());
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
            int currentT = -1;
            ImageProcessor ip2 = null;
            for (int row = 0; row < c0Size; row++) {
                progress.updateProgress(row, c0Size);
                int c0t = (int) Math.round(c0CudaData[f].get(row)[0]);
                if (c0t > currentT) {
                    currentT = c0t;
                    if (channel2 != null) {
                        ip2 = preProcess(channel2.getProcessor(currentT + 1), c2sigma);
                    }
                }
                double x0 = c0CudaData[f].get(row)[1];
                double y0 = c0CudaData[f].get(row)[2];
                double mag = c0CudaData[f].get(row)[3];
                double fit = c0CudaData[f].get(row)[4];
                if (fit > UserVariables.getCurveFitTol()) {
                    IsoGaussian g1 = new IsoGaussian(c0t - startSlice, x0, y0, mag, UserVariables.getSigEstRed(), UserVariables.getSigEstRed(), fit, null, -1, null);
                    g1.setColocalisedParticle(getC2Gaussian(x0, y0, ip2));
                    particles.addDetection(c0t - startSlice, g1);

                }
            }
            progress.dispose();
        }
        updateTrajs(particles, spatialRes, update);
        return particles;
    }
}
