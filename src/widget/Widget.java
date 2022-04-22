package widget;
import java.awt.*;
import javax.swing.JFrame;
import javax.swing.JPanel;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

public class Widget extends JPanel{ 
    //Gravity constant in Solar-Radii, Solar Masses, and Kilometers/second
    static final double G = 1.90809 * 100000;
    //Time in seconds
    static double time = 0;
    static double tStep = 1;
    static double tFinal = 3600 * 24 * 365.0 * 1000;
    //Coordinates used for rendering output
    static double[] pos1,pos2;
    
    public static void main(String[] args) {
        //Initialize render output
        JFrame frame = new JFrame();
        frame.getContentPane().add(new Widget());
        frame.setSize(1000,1000);
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);      
        
        //Velocity and displacement from the origin
        Vector2D v,d;
        
        //Declare star 1
        v = new Vector2D(0,29.78);
        d = new Vector2D(214.93946938362,0);
        Star s1 = new Star(v,d,1/332946.0);
        
        //Declare star 2        
        v = new Vector2D(0,0);
        d = new Vector2D(0,0);
        BlackHole s2 = new BlackHole(v,d,1);
        
        //Orbit/Accuracy counting trackers
        double tmp = 0;
        boolean off = true;
        int count = 0;
        
        //Main update loop
        while (time <= tFinal) {
            
            //Saves the distance in the x direction and y direction for later, and computes the hypotenuse (total distance)
            double distSq = s2.getDis().distanceSq(s1.getDis());
            //Defines the timestep, based on the star's distance^2
            tStep = distSq * 0.0000000001;
            //Find the force of gravity on the star
            double GForce = calcGForce(s1, s2, distSq);
            //Determine the coordinate component forces for each star (newton's third law, saves cpu cycles so not computing twice)
            Vector2D Fg1 = getAcceleration(s1.getDis(),s2.getDis(),GForce);
            Vector2D Fg2 = Fg1.negate();
            //Update the velocity of each star given the gravitational force in coordinate component directions (will compute acceleration first)
            s1.updateVelocity(Fg1, tStep);
            s2.updateVelocity(Fg2, tStep);
            //Update the position of each star based off the velocity previously calculated
            s1.updatePosition(tStep);
            s2.updatePosition(tStep);

            //Iterate timestep
            time += tStep;
            //Define rendered coordinates
            pos1 = new double[] {s1.getDis().getX(),s1.getDis().getY()};
            pos2 = new double[] {s2.getDis().getX(),s2.getDis().getY()};
            //Accuracy tracking
            double x1a = s1.getDis().getX();
            if (x1a > tmp) {
                off = true;
            }
            else if (off) {
                System.out.println("Orbit #: " + count);
                System.out.println(tmp / 214.93946938362);
                off = false;
                count++;
            }
            tmp = x1a;
        }
    }
    
    //Distance calculation (stays squared to save expensive Math.sqrt)
    public static double calcDistanceSqrd(double xdis, double ydis) {
        return xdis * xdis + ydis * ydis;
    }
    
    //Universal gravitational equation based on G, masses, and distance squared found earlier
    public static double calcGForce(Star s1, Star s2, double distSqrd) {
        return (G * s1.getMass() * s2.getMass()) / distSqrd;
    }
    
    //Converts the star's distance in each direction and the total gravitational force between stars to compute force in coordinate component form
    public static Vector2D getAcceleration (Vector2D v1, Vector2D v2, double magnitude) {
        return v2.subtract(v1).normalize().scalarMultiply(magnitude);
    }

    public Widget() {
    }
    
    //Graphics
    @Override
    public void paint(Graphics g) {
        g.clearRect(0, 0, 1700, 1000);
        int height = 50;
        int width = 50;
        g.fillArc((int)pos1[0] - (width / 2) + 500, (int)pos1[1] - (height / 2) + 500, width, height, 0, 360);
        g.fillArc((int)pos2[0] - (width / 20) + 500, (int)pos2[1] - (height / 20) + 500, width / 10, height / 10, 0, 360);
        repaint();
    } 
    
}
