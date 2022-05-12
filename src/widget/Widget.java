//If Volume is larger than the roche lobe transfer mass
//Figure out max accretion rate then apply buffer
package widget;
import java.awt.*;
import javax.swing.JFrame;
import javax.swing.JPanel;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

public class Widget extends JPanel{ 
//Gravity constant in Solar-Radii, Solar Masses, and Kilometers/second ^ 2
    static final double G = 1.90809 * 100000;
    //Speedscalar controls the overall "fedelity" of the simulation Default: 1
    static final double SPEEDSCALAR = 100;
    //Length of the simulation in seconds * 3600 s/hr * 24 hr/day * 365 days/year;
    static final double TFINAL = 3600 * 24 * 365;
    //Time in seconds
    static double time = 0;
    static double tStep = 1;
    //Coordinates/ perameters used for rendering output
    static final double ZOOM = 30;
    static double[] renderedPositionStar1, renderedPositionStar2;
    static double renderedRadiusStar1,renderedRadiusStar2;
    static double renderedRocheLobe;
    static double renderedMassStar1;
    //ORBIT CHARACTERISTICS
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
        
        double a = 5;//214.93946938362;
        double m1 = 1;//1/332954.355179;
        double m2 = 10;
        double e = 0.6;
        double d1 = calcDisplacementFromBarycenter(a,m1,m2);
        double d2 = a - d1;
        
        //Declare stars using the orbital perameters
        Star s1 = InitStar(a,m1,m2,e,d1,false);
        Star s2 = InitStar(a,m2,m1,e,d2,true);
        
        boolean[] result = run(s1,s2,m1);
        System.out.println("The star lost mass to the black hole: " + result[0]);
        System.out.println("The star lost sufficient mass that neuclear fusion is no longer possible: " + result[1]);
    }
    //Completes one full simulation scenario
    public static boolean[] run(Star s1 , Star s2, double startMassStar1) {
        while (time <= TFINAL) {
            
            //Saves the distance in the x direction and y direction for later, and computes the norm
            Vector2D dis1 = s1.getDis();
            Vector2D dis2 = s2.getDis();
            double distSq = dis2.distanceSq(dis1);
            
            //Defines the timestep, based on the star's distance^2
            tStep =  distSq * 0.0001 * SPEEDSCALAR;
            
            //Find the force of gravity on the star
            double GForce = calcGForce(s1, s2, distSq);
            
            //Determine the coordinate component forces for each star (newton's third law, saves cpu cycles so not computing twice)
            Vector2D Fg1 = calcAcceleration(dis1,dis2,GForce);
            Vector2D Fg2 = Fg1.negate();
            
            //Update the velocity of each star given the gravitational force by first computing acceleration
            s1.updateVelocity(Fg1, tStep);
            s2.updateVelocity(Fg2, tStep);
            
            //Update the position of each star based off the velocity previously calculated
            s1.updatePosition(tStep);
            s2.updatePosition(tStep);
            
            //Determines the size of the roche lobe based off of orbital perameters
            double q = s1.getMass()/s2.getMass();
            double cosV = getCosTrueAnomaly(s1,s2);
            double e = getEccentricityVector(s1,s2).getNorm();
            double A = alphaFunction(e,cosV);
            double rlobe = eggletonVolumeEquivelantRocheLobe(Math.sqrt(distSq),q) * sepinskyCorrectionFunction(A,q);
            
            //Define rendered coordinates / perameters
            Vector2D offset = findBarycenter(s1,s2);
            renderedPositionStar1 = new double[] {dis1.getX() - offset.getX(),dis1.getY() - offset.getY()};
            renderedPositionStar2 = new double[] {dis2.getX() - offset.getX(),dis2.getY() - offset.getY()};
            renderedRadiusStar1 = s1.getRadius() * 1;
            renderedRadiusStar2 = s2.getRadius() * 10;
            renderedMassStar1 = s1.getMass();
            renderedRocheLobe = rlobe;
            //Iterate timestep
            time += tStep;
            
            /*qa    DEBUG TRACERS GO HERE
            if (Math.floor(time) % 3600 == 0) {
            System.out.println("Eccentricity: " + getEccentricityVector(s1,s2));
            System.out.println("Alpha value: " + A);
            System.out.println("Roche Lobe Size: " + rlobe);
            System.out.println("True anomaly (cos): " + cosV);
            }
            //END DEBUG TRACERS */
            
            //Check if star mass is less than what is required for nuclear fusion to continue
            if (s1.getMass() < 0.00000001 )
                return new boolean[] {true,true};
        }
        if (s1.getMass() < startMassStar1)
            return new boolean[] {true,false};
        return new boolean[] {false, false};
    }
    
    //Initializes a star with given starting conditions
    private static Star InitStar(double a, double m1, double m2, double e, double displacement, boolean b) {
        if (b)
            return new BlackHole(new Vector2D(0, 0 * calcStartVelocity(a,m1,m2,e)),new Vector2D(displacement,0),m1);
        return new Star(new Vector2D(0,-1 * calcStartVelocity(a,m1,m2,e)),new Vector2D(-1 * displacement,0),m1);
    }
    
    //Universal gravitational equation based on G, masses, and distance squared found earlier
    public static double calcGForce(Star s1, Star s2, double distSqrd) {
        return (G * s1.getMass() * s2.getMass()) / (distSqrd * 696000.0);
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
        if (renderedMassStar1 > 40)
            g.setColor(Color.BLUE);
        else if (renderedMassStar1 > 10)
            g.setColor(Color.CYAN);
        else if (renderedMassStar1 > 1.5)
            g.setColor(Color.WHITE);
        else if (renderedMassStar1 > 1)
            g.setColor(Color.YELLOW);
        else if (renderedMassStar1 > 0.5)
            g.setColor(Color.ORANGE);
        else
            g.setColor(Color.RED);
        g.fillArc((int)(renderedPositionStar1[0] * ZOOM) - (int)(renderedRadiusStar1 * ZOOM) + 500, (int)(renderedPositionStar1[1] * ZOOM) - (int)(renderedRadiusStar1 * ZOOM) + 500, (int)(renderedRadiusStar1 * ZOOM) * 2, (int)(renderedRadiusStar1 * ZOOM) * 2, 0, 360);
        g.setColor(Color.pink);
        g.drawArc((int)(renderedPositionStar1[0] * ZOOM) - (int)(renderedRocheLobe * ZOOM) + 500, (int)(renderedPositionStar1[1] * ZOOM) - (int)(renderedRocheLobe * ZOOM) + 500,(int)(renderedRocheLobe * ZOOM) * 2,(int)(renderedRocheLobe * ZOOM) * 2,0,360);
        g.setColor(Color.black);
        g.fillArc((int)(renderedPositionStar2[0] * ZOOM) - (int)(renderedRadiusStar2) + 500, (int)(renderedPositionStar2[1] * ZOOM) - (int)(renderedRadiusStar2) + 500, (int)(renderedRadiusStar2) * 2, (int)(renderedRadiusStar2) * 2, 0, 360);
        repaint();
    }
    public static Vector2D findBarycenter(Star s1, Star s2) {
        double xBar = (s1.getMass() * s1.getDis().getX() + s2.getMass() * s2.getDis().getX()) / (s1.getMass() + s2.getMass());
        double yBar = (s1.getMass() * s1.getDis().getY() + s2.getMass() * s2.getDis().getY()) / (s1.getMass() + s2.getMass());
        return new Vector2D(xBar,yBar);
    }
    
    public static double eggletonVolumeEquivelantRocheLobe(double seperation, double q) {
        return (seperation * 0.49 * Math.pow(q,2.0/3.0))/(0.6 * Math.pow(q,2.0/3.0) + Math.log(1 + Math.pow(q, 1.0/3.0)));
    }
    public static Vector2D getEccentricityVector(Star s1, Star s2) {
        Vector2D v = s1.getVel().subtract(s2.getVel());
        Vector2D r = s1.getDis().subtract(s2.getDis());
        double meu = G * (s1.getMass() + s2.getMass());
        return r.scalarMultiply((v.dotProduct(v)/meu) - (1/r.getNorm())).subtract(v.scalarMultiply(r.dotProduct(v) / meu));
    }
    public static double getCosTrueAnomaly(Star s1, Star s2) {
        Vector2D r = s1.getDis().subtract(s2.getDis());
        Vector2D e = getEccentricityVector(s1,s2);
        return e.dotProduct(r) / (e.getNorm() * r.getNorm());
    }
    public static double alphaFunction (double eccentricity, double cosV) {
        return Math.pow(1 + eccentricity, 4) / Math.pow(1 + eccentricity * cosV, 3);
    }
    public static double sepinskyCorrectionFunction(double A, double q) {
        double in = Math.log10(A);
        if (in >= 0.2)
            return i0(in) + i1(in) * Math.exp(-i2(in) * Math.pow(Math.log10(q) + i3(in),2));
        return g0(in) + g1(in) * Math.log10(q) + g2(in) * Math.pow(Math.log10(q), 2);
    }
    public static double i0(double A) {
        return (6.3014 * Math.pow(A, 1.3643))/(Math.exp(2.3644 * Math.pow(A,0.70748)) - (1.4413 * Math.exp(-0.0000184 * Math.pow(A, -4.5693))));
    }
    public static double i1(double A) {
        return A / (0.0015 * Math.exp(8.84 * Math.pow(A,0.282)) + 15.78);
    }
    public static double i2(double A) {
        return (1 + 0.036 * Math.exp(8.01 * Math.pow(A, 0.879))) / (0.105 * Math.exp(7.91 * Math.pow(A, 0.879)));
    }
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
}