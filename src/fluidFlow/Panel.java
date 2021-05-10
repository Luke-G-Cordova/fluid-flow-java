

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
	
	public int[][] colors = {
//			{255, 255, 255},
			{253,254,2},
			{11,255,1},
			{254,0,246},
			{255, 0, 0}
//			{255, 0  , 255},
//			{0  , 255, 255}
	};//255, 255, 255, 0, 0, 0
	public int[] col = colors[3];
	public Panel() {
		this.setBackground(new Color(255, 255, 255));
		Loop loop = new Loop(this);
		this.addMouseListener(new MouseListener() {
			@Override
			public void mouseClicked(MouseEvent e) {
				col = colors[(int) (Math.random()*colors.length)];
//				fluid.addDye(e.getX()/SCALE, e.getY()/SCALE, 255f);
			}
			@Override
			public void mousePressed(MouseEvent e) {
//				col = colors[(int) (Math.random()*colors.length)];
				prevX = e.getX()/SCALE;
				prevY = e.getY()/SCALE;
			}
			@Override
			public void mouseReleased(MouseEvent e) {}

			@Override
			public void mouseEntered(MouseEvent e) {}

			@Override
			public void mouseExited(MouseEvent e) {}
			
		});
		this.addMouseMotionListener(new MouseMotionListener() {

			@Override
			public void mouseDragged(MouseEvent e) {
				if(fluid!=null) {
					float amtX = e.getX()/SCALE - prevX;
					float amtY = e.getY()/SCALE - prevY;
//					, col[0], col[1], col[2]
					fluid.addDye(e.getX()/SCALE, e.getY()/SCALE, 255f);
//					fluid.addDye(e.getX()/SCALE, e.getY()/SCALE, 255f);
//					System.out.println("r: "+color[0]+" g: "+color[1]+" b: "+color[2]);
					fluid.addVelocity(e.getX()/SCALE, e.getY()/SCALE, amtX, amtY);
					
					prevX = e.getX()/SCALE;
					prevY = e.getY()/SCALE;
				}
//				if(fluid!=null) {
//					
//					fluid.addDye(e.getX()/SCALE, e.getY()/SCALE, 255);
//					
//					float amtX = e.getX()/SCALE - prevX;
//					float amtY = e.getY()/SCALE - prevY;
//					
//					fluid.addVelocity(e.getX()/SCALE, e.getY()/SCALE, amtX, amtY);
//					
//					prevX = e.getX()/SCALE;
//					prevY = e.getY()/SCALE;
//				}
			}
			@Override
			public void mouseMoved(MouseEvent e) {
//				if(fluid!=null) {
//					int[] color = {255, 255, 255};
//					fluid.addDye(e.getX()/SCALE, e.getY()/SCALE, 255, color);
//					
//					float amtX = e.getX()/SCALE - prevX;
//					float amtY = e.getY()/SCALE - prevY;
//					
//					fluid.addVelocity(e.getX()/SCALE, e.getY()/SCALE, amtX, amtY);
//					
//					prevX = e.getX()/SCALE;
//					prevY = e.getY()/SCALE;
//				}
			}
		});
	}
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		if(fluid!=null) {
			fluid.step();
			fluid.drawDens(g);
//			fluid.drawVel(g);
			
		}
	}

}
