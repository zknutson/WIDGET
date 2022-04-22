/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package widget;

import org.apache.commons.math3.geometry.Vector;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

/**
 *
 * @author Majesticc
 */
public class BlackHole extends Star {
    public BlackHole(Vector2D velocity, Vector2D displacement, double mass) {
        super(velocity,displacement,mass);
    }
}
