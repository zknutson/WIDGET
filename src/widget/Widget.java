//If Volume is larger than the roche lobe transfer mass
//Figure out max accretion rate then apply buffer

package widget;
import java.awt.*;
import javax.swing.JFrame;
import javax.swing.JPanel;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

public class Widget extends JPanel{ 
//Gravity constant in Solar-Radii, Solar Masses, and Kilometers/second
    static final double G = 1.90809 * 100000;
    //Speedscalar controls the overall "fedelity" of the simulation Default: 1
    static final double SPEEDSCALAR = 1;
    //Length of the simulation in seconds * 3600 hrs * 24 days * 365 years;
    static final double TFINAL = 3600 * 24 * 365.0 * 1000;
    //Time in seconds
    static double time = 0;
    static double tStep = 1;
    //Coordinates used for rendering output
    static double[] pos1,pos2;
    static double r1,r2;
    static double rlobe;
    //ORBIT CHARACTERISTICS
    static double a,m1,m2,e,d1,d2;
    // a = Start distance between the two objects (perihelion) in solar radii
    // m1 and m2 = the start masses for the two objects (star and blackhole respectively) in solar masses
    //e = eccentricity of the orbit (%)
    //d1 and d2 = distances from the barycenter for objects 1 and 2
    
    public static void main(String[] args) {
        //System.out.println(sepinskyCorrectionFunction(Acurr,qcurr));  
        //System.out.println(eggletonVolumeEquivelantRocheLobe(1,1));
        //Initialize render output
        JFrame frame = new JFrame();
        InitRenderOutput(frame);     
        
        a = 50;
        m1 = 5;
        m2 = 5;
        e = 0;
        d1 = calcDisplacementFromBarycenter(a,m1,m2);
        d2 = a - d1;
        
        //Declare stars using the orbital perameters
        Star s1 = InitStar(a,m1,m2,e,d1,false);
        Star s2 = InitStar(a,m2,m1,e,d2,true);
        
        run(s1,s2);
    }
    
    //Completes one full simulation scenario
    public static void run(Star s1 , Star s2) {
        //Orbit/Accuracy counting trackers
        //double tmp = 0;
        //boolean off = true;
        //int count = 0;
        double count = 0;
        double mass = 0;
        while (time <= TFINAL) {
            //Saves the distance in the x direction and y direction for later, and computes the norm
            Vector2D dis1 = s1.getDis();
            Vector2D dis2 = s2.getDis();
            double distSq = dis2.distanceSq(dis1);
            //Defines the timestep, based on the star's distance^2
            tStep =  distSq * 0.00000000001 * SPEEDSCALAR;
            //Find the force of gravity on the star
            double GForce = calcGForce(s1, s2, distSq);
            //Determine the coordinate component forces for each star (newton's third law, saves cpu cycles so not computing twice)
            Vector2D Fg1 = calcAcceleration(dis1,dis2,GForce);
            Vector2D Fg2 = Fg1.negate();
            //Update the velocity of each star given the gravitational force in coordinate component directions (will compute acceleration first)
            s1.updateVelocity(Fg1, tStep);
            s2.updateVelocity(Fg2, tStep);
            //Update the position of each star based off the velocity previously calculated
            s1.updatePosition(tStep);
            s2.updatePosition(tStep);
            rlobe = eggletonVolumeEquivelantRocheLobe(Math.sqrt(distSq),m1/m2);
            //Iterate timestep
            time += tStep;
            //Define rendered coordinates
            Vector2D offset = findBarycenter(s1,s2);
            pos1 = new double[] {dis1.getX() - offset.getX() ,dis1.getY() - offset.getY()};
            pos2 = new double[] {dis2.getX() - offset.getX(),dis2.getY() - offset.getY()};
            r1 = s1.getRadius() * 2 + 5;
            r2 = s2.getRadius() * 10;
            //Accuracy tracking
            /*double x1a = s1.getDis().getX();
            if (x1a < tmp) {
                off = true;
            }
            else if (off) {
                System.out.println("Orbit #: " + count);
                System.out.println(tmp / 214.93946938362);
                off = false;
                count++;
            }
            tmp = x1a;*/
        }
    }
    
    //Initializes a star with given starting conditions
    private static Star InitStar(double a, double m1, double m2, double e, double displacement, boolean b) {
        if (b) {
            return new BlackHole(new Vector2D(0,calcStartVelocity(a,displacement,m2,e)),new Vector2D(displacement,0),m1);
        }
        return new Star(new Vector2D(0,-1 * calcStartVelocity(a,displacement,m2,e)),new Vector2D(-1 * displacement,0),m1);
    }
    
    //Universal gravitational equation based on G, masses, and distance squared found earlier
    public static double calcGForce(Star s1, Star s2, double distSqrd) {
        return (G * s1.getMass() * s2.getMass()) / distSqrd;
    }
    //Converts the star's distance in each direction and the total gravitational force between stars to compute force in coordinate component form
    public static Vector2D calcAcceleration (Vector2D v1, Vector2D v2, double magnitude) {
        return v2.subtract(v1).normalize().scalarMultiply(magnitude);
    }
    //Using a simple "center of mass" equation, find the distance each point is from the COM (barycenter)
    private static double calcDisplacementFromBarycenter(double a, double m1, double m2) {
        return (a) / (1 + (m1/m2));
    }
    //Given relevant perameters, calculate the required velocity at perihelion to establish the desired starting orbit
    public static double calcStartVelocity(double a, double radiusFromBarycenter, double mOther, double eccentricity) {
        return 1.2 * Math.sqrt((G * mOther * radiusFromBarycenter) / (a * a));
    }
    
    //Graphics
    private static void InitRenderOutput(JFrame frame) {
        frame.getContentPane().add(new Widget());
        frame.setSize(1000,1000);
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);      
    } 
    @Override
    public void paint(Graphics g) {
        g.clearRect(0, 0, 1000, 1000);
        g.setColor(Color.GRAY);
        g.fillRect(0, 0, 1000, 1000);
        if (m1 > 40) {
            g.setColor(Color.BLUE);
        }
        else if (m1 > 10) {
            g.setColor(Color.CYAN);
        }
        else if (m1 > 1.5) {
            g.setColor(Color.WHITE);
        }
        else if (m1 > 1) {
            g.setColor(Color.YELLOW);
        }
        else if (m1 > 0.5) {
            g.setColor(Color.ORANGE);
        }
        else {
            g.setColor(Color.RED);
        }
        g.fillArc((int)pos1[0] - (int)(r1) + 500, (int)pos1[1] - (int)(r1) + 500, (int)(r1 * 2), (int)(r1 * 2), 0, 360);
        g.setColor(Color.pink);
        g.drawArc((int)pos1[0] - (int)(rlobe) + 500, (int)pos1[1] - (int)(rlobe) + 500,(int)rlobe * 2,(int)rlobe * 2,0,360);
        g.setColor(Color.black);
        g.fillArc((int)pos2[0] - (int)(r2) + 500, (int)pos2[1] - (int)(r2) + 500, (int)(r2 * 2), (int)(r2 * 2), 0, 360);
        repaint();
    }
    //VALIDATED
    public static double sepinskyCorrectionFunction(double A, double q) {
        double in = Math.log10(A);
        return j0(in) + (j1(in) * Math.exp(-j2(in) * Math.pow(Math.log10(q),j3(in))));
    }
    //VALIDATED
    public static double j0(double A) {
        return (1.895 * Math.pow(A,0.837))/(Math.exp(1.636 * Math.pow(A,0.789)) - 1);
    }
    //VALIDATED
    public static double j1(double A) {
        return (4.3 * Math.pow(A, 0.98))/(Math.exp(2.5 * Math.pow(A,0.66)) + 4.7);
    }
    //VALIDATED
    public static double j2(double A) {
        return (1.0)/((8.8 * Math.exp(-2.95 * Math.pow(A,0.76))) + (1.64 * Math.exp(-0.03 * Math.pow(A,0.76))));
    }
    //VALIDATED
    public static double j3(double A) {
        return 0.256 * Math.exp(-1.33 * Math.pow(A, 2.9)) * (5.5 * Math.exp(1.33 * Math.pow(A,2.9)) + 1);
    }
    //VALIDATED
    public static double eggletonVolumeEquivelantRocheLobe(double seperation, double q) {
        return (seperation * 0.49 * Math.pow(q,2.0/3.0))/(0.6 * Math.pow(q,2.0/3.0) + Math.log(1 + Math.pow(q, 1.0/3.0)));
    }
    public static double alphaFunction (double eccentricity, double cosV) {
        return Math.pow(1 + eccentricity, 4) / Math.pow(1 + eccentricity * cosV, 3);
    }
    public static Vector2D findBarycenter(Star s1, Star s2) {
        double xBar = (s1.getMass() * s1.getDis().getX() + s2.getMass() * s2.getDis().getX()) / (s1.getMass() + s2.getMass());
        double yBar = (s1.getMass() * s1.getDis().getY() + s2.getMass() * s2.getDis().getY()) / (s1.getMass() + s2.getMass());
        return new Vector2D(xBar,yBar);
    }
    public static double apoPeriastronChecker(Star s1,Star s2) {
        //System.out.println("VELOCITY: " + s1.getVel());
        //System.out.println("POSITION: " + s1.getDis());
        //System.out.println("BARYCENTER POS: " + findBarycenter(s1,s2));
        return s1.getDis().subtract(findBarycenter(s1,s2)).dotProduct(s1.getVel());
    }
}