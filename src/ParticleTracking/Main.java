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

import Detection.Filament_Detector;
import ij.IJ;
import ij.process.ImageProcessor;

public class Main {

//    public static void main(String args[]) {
//        Particle_Mapper instance = new Particle_Mapper();
//        instance.run(null);
//        System.exit(0);
//    }
//    public static void main(String args[]) {
//        TestGenerator tg = new TestGenerator();
//        tg.generateMulti(10, 512, 512, 100);
////        tg.generateMulti(40, 10, 512, 512, tg.generateNuclei(10, 512, 512, 24, 36));
//        System.exit(0);
//    }
    public static void main(String args[]) {
        GPUAnalyse instance = new GPUAnalyse();
        instance.run(null);
        System.exit(0);
    }
//    public static void main(String args[]) {
//        int radius = 5;
//        double[] sigmas = new double[10];
//        for (int i = 0; i < sigmas.length; i++) {
//            sigmas[i] = 5.0 / (i + 1);
//        }
//        ImageProcessor ip = IJ.openImage().getProcessor();
//        for (int sIndex = 0; sIndex < sigmas.length; sIndex++) {
//            Filament_Detector instance = new Filament_Detector(sigmas[sIndex], radius);
//            instance.convolve(ip.duplicate(), "C://Users/barryd/filament_detector_debug");
////            instance.outputKernel();
////            instance.laplacianOfGaussian(ip.duplicate(), "C://Users/barryd/blob_detector_debug");
//        }
//        System.exit(0);
//    }

//        public static void main(String args[]) {
//        Particle_Colocaliser instance = new Particle_Colocaliser();
//        instance.run(null);
//        System.exit(0);
//    }
}
