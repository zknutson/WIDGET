package widget;

public class Star {
    //Acceleration, velocity, and displacement values in ordered pairs (x,y)
    private double[] a = new double[2],v = new double[2],d = new double[2];
    //Mass in solar masses
    private double mass;
    
    //Constructor
    public Star(double[] velocity, double[] displacement, double mass) {
        v = velocity;
        d = displacement;
        this.mass = mass;
    }
    
    //Main timestepping methods
    public void updateVelocity(double[] Fg, double tStep) {
        //Find the acceleration in each coordinate component direction, based on the gravitational force and the star's distance in each direction
        this.a[0] = Fg[0] / mass;
        this.a[1] = Fg[1] / mass;
        v[0] += a[0] * tStep;
        v[1] += a[1] * tStep;
    }
    public void updatePosition(double tStep) {
        d[0] += v[0] * tStep - 0.5 * a[0] * tStep * tStep;
        d[1] += v[1] * tStep - 0.5 * a[1] * tStep * tStep; 
    }

    //Getters and Setters
    public double[] getDis() {
        return d;
    }

    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }
    
}