package ini.cx3d;

import ini.cx3d.localBiology.factory.NeuriteElementFactory;
import ini.cx3d.localBiology.factory.SomaElementFactory;
import ini.cx3d.physics.factory.*;
import ini.cx3d.physics.interfaces.PhysicalBond;
import ini.cx3d.simulations.ECM;
import ini.cx3d.spatialOrganization.NewDelaunayTest;
import ini.cx3d.spatialOrganization.factory.OpenTriangleOrganizerFactory;
import ini.cx3d.spatialOrganization.interfaces.OpenTriangleOrganizer;
import ini.cx3d.swig.simulation.JavaUtilT_PhysicalNode;
import ini.cx3d.synapses.factory.PhysicalBoutonFactory;
import ini.cx3d.synapses.factory.PhysicalSpineFactory;
import ini.cx3d.utilities.Matrix;

import java.awt.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import static ini.cx3d.utilities.Matrix.concat;
import static ini.cx3d.utilities.Matrix.randomNoise;

/**
 * provides functionality that has not been implemented yet in C++
 * especially static methods as they can't be handled by SWIG directors
 */
public class JavaUtil2 extends ini.cx3d.swig.simulation.JavaUtil2 {

   public double[] matrixRandomNoise3(double k){
        return Matrix.randomNoise(k, 3);
    }

    public double exp(double d) { return Math.exp(d); }
    public double sqrt(double d) { return Math.sqrt(d); }
    public double cos(double d) { return Math.cos(d); }
    public double sin(double d) { return Math.sin(d); }
    public double asin(double d) { return Math.asin(d); }
    public double acos(double d) { return Math.acos(d); }
    public double atan2(double d, double d1) { return Math.atan2(d, d1); }
    public double cbrt(double d) { return Math.cbrt(d); }
    public ini.cx3d.physics.interfaces.PhysicalCylinder newPhysicalCylinder() {return PhysicalCylinderFactory.create();}
    public double getRandomDouble1(){
        return getRandomDouble();
    }
    public double matrixNextRandomDouble(){
        return Matrix.getRandomDouble();
    }
    public ini.cx3d.physics.interfaces.PhysicalNode newPhysicalNode() {return PhysicalNodeFactory.create();}
    public ini.cx3d.spatialOrganization.SpatialOrganizationNodeMovementListener newPhysicalNodeMovementListener() {return PhysicalNodeMovementListenerFactory.create();}
    public ini.cx3d.physics.interfaces.PhysicalSphere newPhysicalSphere() {return PhysicalSphereFactory.create();}
    public ini.cx3d.localBiology.interfaces.SomaElement newSomaElement() {return SomaElementFactory.create();}
    public ini.cx3d.localBiology.interfaces.NeuriteElement newNeuriteElement() {return NeuriteElementFactory.create();}
    public ini.cx3d.synapses.interfaces.PhysicalSpine newPhysicalSpine(ini.cx3d.physics.interfaces.PhysicalObject po, double[] origin, double length) {
        return PhysicalSpineFactory.create(po, origin, length);
    }
    public ini.cx3d.synapses.interfaces.PhysicalBouton newPhysicalBouton(ini.cx3d.physics.interfaces.PhysicalObject po, double[] origin, double length) {
        return PhysicalBoutonFactory.create(po, origin, length);
    }
    public PhysicalBond newPhysicalBond(ini.cx3d.physics.interfaces.PhysicalObject a, double[] positionOnA, ini.cx3d.physics.interfaces.PhysicalObject b , double[] positionOnB, double restingLength, double springConstant) {
        return PhysicalBondFactory.create(a, positionOnA, b, positionOnB, restingLength, springConstant);
    }

    public Color getRandomColor(){
        Color c = new Color((float) getRandomDouble(),(float) getRandomDouble(),(float) getRandomDouble(),0.7f);
        return c;
    }

    public void initPhysicalNodeMovementListener(){
        ini.cx3d.swig.simulation.PhysicalNodeMovementListener.setMovementOperationId((int) (10000 * getRandomDouble()));
    }


    // **************************************************************************
    // Random Number
    // **************************************************************************
    static Random random = new Random();

    /**
     * @return a random number between, from uniform probability 0 and 1;
     */
    public static double getRandomDouble(){
        return random.nextDouble();
    }

    /**
     * returns a random number from gaussian distribution
     * @param mean
     * @param standardDeviation
     * @return
     */
    public double getGaussianDouble(double mean, double standardDeviation){
        return mean + standardDeviation*random.nextGaussian();
    }



    /**
     * Initialises the random number generator.
     * @param seed
     */
    public static void setRandomSeed(long seed){
        random = new Random(seed);
        Matrix.setRandomSeedTo(seed);
    }

    public void setRandomSeed1(int seed){
        setRandomSeed(seed);
    }

    public double[] foo() {
        return concat(randomNoise(600,2), randomNoise(100,1));
    }
}
