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

    public ImageProcessor laplacianOfGaussian(ImageProcessor input, String directory) {
        ImageProcessor output = input.convertToFloat();
        int height = input.getHeight();
        int width = input.getWidth();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                output.putPixelValue(x, y, laplacianOfGaussianNeighbourhood(input, x, y));
            }
        }
        IJ.saveAs(new ImagePlus("", output), "TIFF", String.format("%s%slog_%1.3f_%d", directory, File.separator, sigma, radius));
        return output;
    }

    private double laplacianOfGaussianNeighbourhood(ImageProcessor input, int x0, int y0) {
        double output = 0.0;
        for (int y = y0 - radius; y <= y0 + radius; y++) {
            double y2 = Math.pow(y - y0, 2.0);
            for (int x = x0 - radius; x <= x0 + radius; x++) {
                double x2 = Math.pow(x - x0, 2.0);
                output += input.getPixelValue(x, y) * lOGOperator(x2, y2);
            }
        }
        return output;
    }

    public void outputKernel() {
        FloatProcessor output = new FloatProcessor(2 * radius + 1, 2 * radius + 1);
        int x0 = radius, y0 = radius;
        for (int y = y0 - radius; y <= y0 + radius; y++) {
            double y2 = Math.pow(y - y0, 2.0);
            for (int x = x0 - radius; x <= x0 + radius; x++) {
                double x2 = Math.pow(x - x0, 2.0);
                output.putPixelValue(x, y, lOGOperator(x2, y2));
            }
        }
        IJ.saveAs(new ImagePlus("", output), "TIFF", String.format("C:/Users/barryd/particle_mapper_debug/kernel_%1.3f_%d", sigma, radius));
    }

    private double lOGOperator(double x, double y) {
        return (-1.0 / (Math.PI * sig4))
                * (1.0 - ((x + y) / (2.0 * sig2)))
                * Math.exp(-(x + y) / (2.0 * sig2));
    }
}
