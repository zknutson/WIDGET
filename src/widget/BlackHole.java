package widget;

import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

public class BlackHole extends Star {
    public BlackHole(Vector2D velocity, Vector2D displacement, double mass) {
        super(velocity,displacement,mass);
        this.radius = 0.25;
    }
    @Override
    public void setMass(double mass) {
        this.mass = mass;
    }
}
