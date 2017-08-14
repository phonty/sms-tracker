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

import ij.process.ImageProcessor;

public class Blob_Detector {

    private static double sigma;
    private double sig2, sig4;
    private static int radius;

    public Blob_Detector() {

    }

    public Blob_Detector(double sigma, int radius) {
        this();
        this.sigma = sigma / Math.sqrt(2.0);
        this.radius = radius;
        sig2 = Math.pow(this.sigma, 2.0);
        sig4 = Math.pow(this.sigma, 4.0);
    }

    public ImageProcessor laplacianOfGaussian(ImageProcessor input) {
        ImageProcessor output = input.convertToFloat();
        float[] inputPix = (float[]) output.getPixels();
        int height = input.getHeight();
        int width = input.getWidth();
        double[] kernel = createKernel();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                output.putPixelValue(x, y, laplacianOfGaussianNeighbourhood(inputPix, x, y, kernel, new int[]{width, height}));
            }
        }
        return output;
    }

    private double laplacianOfGaussianNeighbourhood(float[] pixels, int x, int y, double[] kernel, int[] dims) {
        double output = 0.0;
        int x0 = x - radius;
        int y0 = y - radius;
        int kernelSize = 2 * radius + 1;
        for (int j = y0 < 0 ? 0 : y0; j <= y + radius && j < dims[1]; j++) {
            int iOffset = j * dims[0];
            int kOffset = (j - y0) * kernelSize - x0;
            for (int i = x0 < 0 ? 0 : x0; i <= x + radius && i < dims[0]; i++) {
                output += pixels[iOffset + i] * kernel[i + kOffset];
            }
        }
        return output;
    }

    public double[] createKernel() {
        int size = 2 * radius + 1;
        int size2 = size * size;
        double[] kernel = new double[size2];
        int x0 = radius, y0 = radius;
        double sum = 0.0;
        for (int y = 0; y < size; y++) {
            double y2 = Math.pow(y - y0, 2.0);
            int offset = y * size;
            for (int x = 0; x < size; x++) {
                double x2 = Math.pow(x - x0, 2.0);
                kernel[x + offset] = lOGOperator(x2, y2);
                sum += Math.abs(kernel[x + offset]);
            }
        }
        for (int i = 0; i < size2; i++) {
            kernel[i] /= sum;
        }
        return kernel;
    }

    private double lOGOperator(double x, double y) {
        return (-1.0 / (Math.PI * sig4))
                * (1.0 - ((x + y) / (2.0 * sig2)))
                * Math.exp(-(x + y) / (2.0 * sig2));
    }
}
