

import java.awt.Color;
import java.awt.Graphics;
import java.awt.MouseInfo;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import javax.swing.JPanel;

public class Panel extends JPanel{
	
	public Fluid fluid;
	public int prevX=0, prevY=0;

	public int N;
	public int SCALE = 5;

	public Panel() {

		this.setBackground(new Color(255, 255, 255));

		// initialize new Loop() as a variable if you want access ever
		new Loop(this);

		this.addMouseMotionListener(new MouseMotionListener() {
			@Override
			public void mouseDragged(MouseEvent e) {
				if(fluid!=null) {
					float amtX = e.getX()/SCALE - prevX;
					float amtY = e.getY()/SCALE - prevY;
					fluid.addDye(e.getX()/SCALE, e.getY()/SCALE, 255f);
					fluid.addVelocity(e.getX()/SCALE, e.getY()/SCALE, amtX, amtY);
					
					prevX = e.getX()/SCALE;
					prevY = e.getY()/SCALE;
				}
			}
			// this function is not being used
			@Override public void mouseMoved(MouseEvent e) {}
		});
	}
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		if(fluid!=null) {
			fluid.step();
			fluid.drawDens(g);
		}
	}
}
