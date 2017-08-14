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
package Detection;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.io.File;

/**
 *
 * @author Dave Barry <david.barry at crick.ac.uk>
 */
public class Filament_Detector {

    private static double sigma;
    private double sig2;
    private static int radius;

    public Filament_Detector() {

    }

    public Filament_Detector(double sigma, int radius) {
        this();
        this.sigma = sigma;
        this.radius = radius;
        sig2 = Math.pow(this.sigma, 2.0);
        createKernels();
    }

    public ImageProcessor convolve(ImageProcessor input, String directory) {
        ImageProcessor output = input.convertToFloat();
        int height = input.getHeight();
        int width = input.getWidth();
        double[][][] kernels = createKernels();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                output.putPixelValue(x, y, convolveNeighbourhood(input, x, y, kernels));
            }
        }
        IJ.saveAs(new ImagePlus("", output), "TIFF", String.format("%s%sfd_%1.3f_%d", directory, File.separator, sigma, radius));
        return output;
    }

    private double convolveNeighbourhood(ImageProcessor input, int x, int y, double[][][] kernels) {
        double output = 0.0;
        int x0 = x - radius;
        int y0 = y - radius;
        for (int j = y0; j <= y + radius; j++) {
            for (int i = x0; i <= x + radius; i++) {
                output += input.getPixelValue(i, j) * Math.sqrt(Math.pow(kernels[0][i - x0][j - y0], 2.0) + Math.pow(kernels[1][i - x0][j - y0], 2.0));
            }
        }
        return output;
    }

    private double[][][] createKernels() {
        int size = 2 * radius + 1;
        double[][][] kernels = new double[2][size][size];
        FloatProcessor[] kernelImages = new FloatProcessor[2];
        kernelImages[0] = new FloatProcessor(size, size);
        kernelImages[1] = new FloatProcessor(size, size);
        int x0 = radius, y0 = radius;
        double sumx = 0.0, sumy = 0.0;
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                kernels[0][x][y] = operator(x, x0);
                kernels[1][x][y] = operator(y, y0);
                sumx += kernels[0][x][y];
                sumy += kernels[1][x][y];
            }
        }
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                kernels[0][x][y] /= sumx;
                kernelImages[0].putPixelValue(x, y, kernels[0][x][y]);
                kernels[1][x][y] /= sumy;
                kernelImages[1].putPixelValue(x, y, kernels[1][x][y]);
            }
        }
        IJ.saveAs(new ImagePlus("", kernelImages[0]), "TIFF", String.format("C:/Users/barryd/filament_detector_debug/xkernel_%1.3f_%d", sigma, radius));
        IJ.saveAs(new ImagePlus("", kernelImages[1]), "TIFF", String.format("C:/Users/barryd/filament_detector_debug/ykernel_%1.3f_%d", sigma, radius));
        return kernels;
    }

    private double operator(double x, double x0) {
        return (1.0 / Math.sqrt(2.0 * Math.PI * sig2)) * Math.exp(-Math.pow(x - x0, 2.0) / (2.0 * sig2));
    }
}
