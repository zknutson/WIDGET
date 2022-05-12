package widget;

import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

public class Star {
    //Acceleration, velocity, and displacement values in ordered pairs (x,y)
    private Vector2D a,v,d;    
    //Mass in solar masses
    protected double mass;
    protected double radius;
    
    //Constructor
    public Star(Vector2D velocity, Vector2D displacement, double mass) {
        a = new Vector2D(0,0);
        v = velocity;
        d = displacement;
        this.mass = mass;
        this.radius = computeRadius();
    }
    
    public double computeRadius() {
        if (this.mass > 1) {
            return Math.pow(this.mass, 0.57);
        }
        return Math.pow(this.mass, 0.88);
    }
    
    //Main timestepping methods
    public void updateVelocity(Vector2D Fg, double tStep) {
        this.a = Fg.scalarMultiply(1/mass);
        this.v = v.add(tStep,a);
    }
    public void updatePosition(double tStep) {
        d = d.add(tStep,v.scalarMultiply(1/696000.0));
        d = d.subtract(0.5 * tStep * tStep,a.scalarMultiply(1/696000.0));
    }

    //Getters and Setters
    public Vector2D getDis() {
        return d;
    }
    public Vector2D getVel() {
        return v;
    }
    public double getMass() {
        return mass;
    }
    public double getRadius() {
        return this.radius;
    }
    public void setMass(double mass) {
        this.mass = mass;
        this.radius = computeRadius();
    }
}