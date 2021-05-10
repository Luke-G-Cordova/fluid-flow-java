

import javax.swing.JFrame;

public class Index {
	public static void main(String[] args) {
		JFrame jf = new JFrame("Fluid Flow");
		jf.setSize(500, 500);
		Panel panel = new Panel();
		jf.add(panel);
		jf.setResizable(false);
		jf.setLocationRelativeTo(null);
		jf.setDefaultCloseOperation(jf.EXIT_ON_CLOSE);
		jf.setVisible(true);
	}
}
