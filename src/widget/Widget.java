package widget;
import java.awt.*;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class Widget extends JPanel{ 
    //Gravity constant in (astronomical units * (astronomical units/s)^2) / Solar masses
    static final double G = 3.964913 * Math.pow(10, -14);
    //Time in seconds
    static double time = 0;
    static double tStep = 1;
    static double tFinal = 3600 * 24 * 365.0 * 1000;
    //Coordinates used for rendering output
    static int x1,x2,y1,y2;

    public static void main(String[] args) {
        
        //Initialize render output
        JFrame frame = new JFrame();
        frame.getContentPane().add(new Widget());
        frame.setSize(1700,1000);
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);      
        
        //Velocity and displacement from the origin, ordered pairs (x,y)
        double[] v,d;
        
        //Declare star 1
        v = new double[] {0,1.9906 * Math.pow(10, -7)};
        d = new double[] {1,0};
        Star s1 = new Star(v,d,1/332946.0);
        
        //Declare star 2        
        v = new double[] {0,0};
        d = new double[] {0,0};
        BlackHole s2 = new BlackHole(v,d,1);
        
        //Orbit/Accuracy counting trackers
        double tmp = 0;
        boolean off = true;
        int count = 0;
        
        //Main update loop
        while (time <= tFinal) {
            
            //Saves the distance in the x direction and y direction for later, and computes the hypotenuse (total distance)
            double xdist = s2.getDis()[0] - s1.getDis()[0];
            double ydist = s2.getDis()[1] - s1.getDis()[1];
            double distSqrd = calcDistanceSqrd(xdist, ydist);
            
            //Defines the timestep, based on the star's distance^2
            tStep = distSqrd * 40;
            //Find the force of gravity on the stars
            double GForce = calcGForce(s1, s2, distSqrd);

            //Determine the coordinate component forces for each star (newton's third law, saves cpu cycles so not computing twice)
            double [] Fg1 = getCoordCompAcceleration(xdist,ydist,GForce);
            double [] Fg2 = new double[] {Fg1[0] * -1, Fg1[1] * -1};

            //Update the velocity of each star given the gravitational force in coordinate component directions (will compute acceleration first)
            s1.updateVelocity(Fg1, tStep);
            s2.updateVelocity(Fg2, tStep);

            //Update the position of each star based off the velocity previously calculated
            s1.updatePosition(tStep);
            s2.updatePosition(tStep);

            //Iterate timestep
            time += tStep;

            //Define rendered coordinates
            x1 = (int)(s1.getDis()[0] * 300);
            x2 = (int)(s2.getDis()[0] * 300);
            y1 = (int)(s1.getDis()[1] * 300);
            y2 = (int)(s2.getDis()[1] * 300);

            //Accuracy tracking
            double x1a = s1.getDis()[0];
            if (x1a > tmp) {
                off = true;
            }
            else if (off) {
                System.out.println("Orbit #: " + count);
                System.out.println(tmp);
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
    public static double[] getCoordCompAcceleration (double i, double j, double magnitude) {
        double x = j/i;
        double y = 1 + x * x;
        double cosarctan = Math.sqrt(y) / (y);
        double sinarctan = cosarctan * x;
        if (i < 0) {
            cosarctan *= -1;
            sinarctan *= -1;
        }
        return new double[] {cosarctan * magnitude, sinarctan * magnitude};
    }
    
    //Graphics
    @Override
    public void paint(Graphics g) {
        g.clearRect(0, 0, 1700, 1000);
        int height = 50;
        int width = 50;
        g.fillArc(x1 - (width / 2) + 1000, y1 - (height / 2) + 500, width, height, 0, 360);
        g.fillArc(x2 - (width / 20) + 1000, y2 - (height / 20) + 500, width / 10, height / 10, 0, 360);
        repaint();
    } 
}
