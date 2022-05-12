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
    static final double SPEEDSCALAR = 10;
    //Length of the simulation in seconds * 3600 hrs * 24 days * 365 years;
    static final double TFINAL = 3600 * 24 * 365.0 * 1000;
    //Time in seconds
    static double time = 0;
    static double tStep = 1;
    //Coordinates used for rendering output
    static double[] pos1 = {0,0} ,pos2 = {0,0};
    static double r1,r2;
    static double rlobe;
    static final double ZOOM = 20;
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
        
        a = 1/0.2054 ;
        m1 = 1;
        m2 = 10;
        e = 0.6;
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
        int count = 0;
        double last = 0;
        Vector2D periastronpos;
        while (time <= TFINAL) {
            //Saves the distance in the x direction and y direction for later, and computes the norm
            Vector2D dis1 = s1.getDis();
            Vector2D dis2 = s2.getDis();
            double distSq = dis2.distanceSq(dis1);
            //Defines the timestep, based on the star's distance^2
            tStep =  distSq * 0.0000000001 * SPEEDSCALAR;
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
            rlobe = eggletonVolumeEquivelantRocheLobe(Math.sqrt(distSq),m1/m2) * sepinskyCorrectionFunction(alphaFunction(getEccentricityVector(s1,s2).getNorm(),getCosTrueAnomaly(s1,s2)),m1/m2);
            Vector2D offset = findBarycenter(s1,s2);
            //double curr = apoPeriastronChecker(s1,s2);
            //if (curr >= 0 && last <= 0) {
            //    periastronpos = offset.subtract(s1.getDis());
            //    System.out.println("Periastron vector from Barycenter: " + periastronpos);
            //}
            //last = curr;
            //if(count % 100000 == 0) {
            //System.out.println("Correction FUnction: " + sepinskyCorrectionFunction(alphaFunction(getEccentricityVector(s1,s2).getNorm(),getCosTrueAnomaly(s1,s2)),m1/m2));
            //System.out.println("Eccentricity Vector: " + getEccentricityVector(s1,s2));
            //System.out.println("True Anomaly: " + getCosTrueAnomaly(s1,s2));
            //System.out.println("A: " + alphaFunction(getEccentricityVector(s1,s2).getNorm(),getCosTrueAnomaly(s1,s2)));
            //}
            //count++;
            //Iterate timestep
            time += tStep;
            //Define rendered coordinates
            pos1 = new double[] {dis1.getX() - offset.getX(),dis1.getY() - offset.getY()};
            pos2 = new double[] {dis2.getX() - offset.getX(),dis2.getY() - offset.getY()};
            r1 = s1.getRadius() * 1;
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
            return new BlackHole(new Vector2D(0, 0 * calcStartVelocity(a,m1,m2,e)),new Vector2D(displacement,0),m1);
        }
        return new Star(new Vector2D(0,-1 * calcStartVelocity(a,m1,m2,e)),new Vector2D(-1 * displacement,0),m1);
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
    public static double calcStartVelocity(double a, double m, double mOther, double eccentricity) {
        return  Math.sqrt(((eccentricity + 1) * (G * (m + mOther))) / a); // * Math.sqrt((G * mOther * radiusFromBarycenter) / (a * a)); 
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
        g.fillArc((int)(pos1[0] * ZOOM) - (int)(r1 * ZOOM) + 500, (int)(pos1[1] * ZOOM) - (int)(r1 * ZOOM) + 500, (int)(r1 * ZOOM) * 2, (int)(r1 * ZOOM) * 2, 0, 360);
        g.setColor(Color.pink);
        g.drawArc((int)(pos1[0] * ZOOM) - (int)(rlobe * ZOOM) + 500, (int)(pos1[1] * ZOOM) - (int)(rlobe * ZOOM) + 500,(int)(rlobe * ZOOM) * 2,(int)(rlobe * ZOOM) * 2,0,360);
        g.setColor(Color.black);
        g.fillArc((int)(pos2[0] * ZOOM) - (int)(r2) + 500, (int)(pos2[1] * ZOOM) - (int)(r2) + 500, (int)(r2) * 2, (int)(r2) * 2, 0, 360);
        //System.out.println(rlobe);
        //System.out.println(pos1[0] + " " + pos1[1]);
        //System.out.println(pos2[0] + " " + pos2[1]);
        repaint();
    }
    /*//VALIDATED
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
    }*/
    //VALIDATED
    public static double sepinskyCorrectionFunction(double A, double q) {
        double in = Math.log10(A);
        if (in >= 0.2) {
        return i0(in) + i1(in) * Math.exp(-i2(in) * Math.pow(Math.log10(q) + i3(in),2));
        }
        return g0(in) + g1(in) * Math.log10(q) + g2(in) * Math.pow(Math.log10(q), 2);
    }
    //VALIDATED
    public static double i0(double A) {
        return (6.3014 * Math.pow(A, 1.3643))/(Math.exp(2.3644 * Math.pow(A,0.70748)) - (1.4413 * Math.exp(-0.0000184 * Math.pow(A, -4.5693))));
    }
    //VALIDATED
    public static double i1(double A) {
        return A / (0.0015 * Math.exp(8.84 * Math.pow(A,0.282)) + 15.78);
    }
    //VALIDATED
    public static double i2(double A) {
        return (1 + 0.036 * Math.exp(8.01 * Math.pow(A, 0.879))) / (0.105 * Math.exp(7.91 * Math.pow(A, 0.879)));
    }
    //VALIDATED
    public static double i3(double A) {
        return 0.991 / (1.38 * Math.exp(-0.035 * Math.pow(A, 0.76)) + 23.0 * Math.exp(-2.89 * Math.pow(A,0.76)));
    }
    public static double g0(double A) {
        return 0.9978 - 0.1229 * A - 0.1273 * A * A;
    }
    public static double g1(double A) {
        return 0.001 + 0.02556 * A;
    }
    public static double g2(double A) {
        return 0.0004 + 0.0021 * A;
    }
    //VALIDATED
    public static double eggletonVolumeEquivelantRocheLobe(double seperation, double q) {
        return (seperation * 0.49 * Math.pow(q,2.0/3.0))/(0.6 * Math.pow(q,2.0/3.0) + Math.log(1 + Math.pow(q, 1.0/3.0)));
    }
    //VALIDATED
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
    //VALIDATED
    public static Vector2D getEccentricityVector(Star s1, Star s2) {
        Vector2D v = s1.getVel().subtract(s2.getVel());
        Vector2D r = s1.getDis().subtract(s2.getDis());
        double meu = G * (s1.getMass() + s2.getMass());
        return r.scalarMultiply((v.dotProduct(v)/meu) - (1/r.getNorm())).subtract(v.scalarMultiply(r.dotProduct(v) / meu));
    }
    //VALIDATED
    public static double getCosTrueAnomaly(Star s1, Star s2) {
        Vector2D r = s1.getDis().subtract(s2.getDis());
        Vector2D e = getEccentricityVector(s1,s2);
        return e.dotProduct(r) / (e.getNorm() * r.getNorm());
    }
}