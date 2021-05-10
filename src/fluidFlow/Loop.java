

public class Loop extends Thread{
	private Panel panel;
	public Loop(Panel panel) {
		this.panel = panel;
		this.start();
	}
	public void run() {
		while(panel.getWidth()==0);{}
		panel.N = panel.getWidth();
		System.out.println(panel.N);
		panel.fluid = new Fluid(panel.N, panel.SCALE, .01f, 0f, 0f);
		while(true) {
			panel.repaint();
		}
	}
}
