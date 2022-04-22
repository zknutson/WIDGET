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
        this.v = v.add(tStep,a);
    }
    public void updatePosition(double tStep) {
        d = d.add(tStep,v);
        d = d.subtract(0.5 * tStep * tStep,a);
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