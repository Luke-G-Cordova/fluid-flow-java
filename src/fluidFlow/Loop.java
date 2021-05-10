

public class Loop extends Thread{
	private Panel panel;

	public Loop(Panel panel) {
		this.panel = panel;
		this.start();
	}

	public void run() {
		// wait for panel to size
		while(panel.getWidth()==0);{}
		
		// this only works for a square space right now, N is the length of one side of the draw space
		panel.N = panel.getWidth();
		
		// create new Fluid()
		panel.fluid = new Fluid(panel.N, panel.SCALE, .01f, 0f, 0f);

		while(true) {
			panel.repaint();
		}
	}
}
