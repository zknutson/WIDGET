package widget;

import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

public class Star {
    //Acceleration, velocity, and displacement values in ordered pairs (x,y)
    private Vector2D a,v,d;    
    //Mass in solar masses
    private double mass, radius;
    
    //Constructor
    public Star(Vector2D velocity, Vector2D displacement, double mass) {
        a = new Vector2D(0,0);
        v = velocity;
        d = displacement;
        this.mass = mass;
        this.radius = getRadius();
    }
    
    //Main timestepping methods
    public void updateVelocity(Vector2D Fg, double tStep) {
        //Find the acceleration in each coordinate component direction, based on the gravitational force and the star's distance in each direction
        this.a = Fg.scalarMultiply(1/mass);
        //System.out.println("GForce: " + Fg);
        //System.out.println("Mass: " + mass);
        //System.out.println("Before: " + v);
        this.v = v.add(tStep,a);
        //System.out.println("After: " + v);
        //System.out.println(tStep);
        //System.out.println("Acceleration:" + a);
        //System.out.println("Velocity:" + v);
    }
    public void updatePosition(double tStep) {
        //System.out.println("Before: " + d);
        d = d.add(tStep,v);
        //System.out.println("Velocity: " + v);
        //System.out.println("Middle: " + d);
        d = d.subtract(0.5 * tStep * tStep,a);
        //System.out.println("Acc: " + a);
        //System.out.println("Position: " + d);
    }

    //Getters and Setters
    public Vector2D getDis() {
        return d;
    }

    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }
    public double getRadius() {
        if (this.mass > 1) {
            return Math.pow(this.mass, 0.57);
        }
        return Math.pow(this.mass, 0.88);
    }
}