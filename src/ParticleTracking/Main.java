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

import ParticleMapping.Particle_Mapper;
import ij.IJ;

/**
 *
 * @author Dave Barry <david.barry at crick.ac.uk>
 */
public class Main {

//    public static void main(String args[]) {
//        Particle_Mapper instance = new Particle_Mapper();
//        instance.run(null);
//        System.exit(0);
//    }
    public static void main(String args[]) {
        TestGenerator tg = new TestGenerator();
        tg.generateMulti(100,30, 512, 512, tg.generateNuclei(20, 512, 512, 24, 36));
        System.exit(0);
    }
}
