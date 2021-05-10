
import java.awt.Color;
import java.awt.Graphics;

public class Fluid {
    private int N;
	private int SCALE;
	private int iter = 1;
	public float[] Vx, Vy, Vx0, Vy0, density, density0;
	private float dt;
	private float[] s;
	private float visc, diff;
	private float r=0, gr=0, b=0;
	private float r1=0, g1=0, b1=0;
	public Fluid(int N, int SCALE, float dt, float visc, float diff) {
		this.N = N/SCALE;
		this.SCALE = SCALE;
		this.dt = dt;
		this.visc = visc;
		this.diff = diff;
		
		this.s = new float[N*N];
		this.Vx = new float[N*N];
		this.Vy = new float[N*N];
		this.Vx0 = new float[N*N];
		this.Vy0 = new float[N*N];
		this.density = new float[N*N];
		this.density0 = new float[N*N];
	}
	public void step(){
	    float visc     = this.visc;
	    float diff     = this.diff;
	    float dt       = this.dt;
	    float[] Vx      = this.Vx;
	    float[] Vy      = this.Vy;
	    float[] Vx0     = this.Vx0;
	    float[] Vy0     = this.Vy0;
	    float[] s       = this.s;
	    float[] density = this.density;
	    
	    diffuse(1, Vx0, Vx, visc, dt);
	    diffuse(2, Vy0, Vy, visc, dt);
	    
	    project(Vx0, Vy0, Vx, Vy);
	    
	    advect(1, Vx, Vx0, Vx0, Vy0, dt);
	    advect(2, Vy, Vy0, Vx0, Vy0, dt);
	    
	    project(Vx, Vy, Vx0, Vy0);
	    
	    diffuse(0, s, density, diff, dt);
	    advect(0, density, s, Vx, Vy, dt);
	}
	public void diffuse(int b, float[] x, float[] x0, float diff, float dt) {
		float a = dt * diff * (N-2)*(N-2);
		for(int k = 0;k<iter;k++) {
			for(int i = 1;i<=N;i++) {
				for(int j = 1;j<=N;j++) {
					x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+
							 x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+4*a);
				}
			}
			setBnd(b, x);
		}
	}
	public void advect(int b, float[] d, float[] d0, float[] Vx, float[] Vy, float dt) {
		int i, j, i0, j0, i1, j1;
		float x, y, s0, t0, s1, t1, dt0;
		dt0 = dt*N;
		for ( i=1 ; i<=N ; i++ ) {
			for ( j=1 ; j<=N ; j++ ) {
				x = i-dt0*Vx[IX(i,j)]; 
				y = j-dt0*Vy[IX(i,j)];
				if (x<0.5) x=0.5f; 
				if (x>N+0.5) x=N+ 0.5f; 
				i0=(int)x; 
				i1=i0+1;
				if (y<0.5) y=0.5f; 
				if (y>N+0.5) y=N+ 0.5f; 
				j0=(int)y; 
				j1=j0+1;
				s1 = x-i0; 
				s0 = 1-s1; 
				t1 = y-j0; 
				t0 = 1-t1;
				d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
			}
		}
		setBnd (b, d);
	}
	public void advectdeleteme(int b, float[] d, float[] d0, float[] Vx, float[] Vy, int[][] colors, int[][] colors0, float dt) {
		int i, j, trueXindex, trueYindex, trueXindexP1, trueYindexP1;
		float trueX, trueY, trueXdeciPtSubFrom1, trueYdeciPtSubFrom1, trueXdeciPt, trueYdeciPt, pixToCoverInTstep;
		pixToCoverInTstep = dt*N;
		for ( i=1 ; i<=N ; i++ ) {
			for ( j=1 ; j<=N ; j++ ) {
				trueX = i-pixToCoverInTstep*Vx[IX(i,j)]; 
				trueY = j-pixToCoverInTstep*Vy[IX(i,j)];
				if (trueX<0.5) trueX=0.5f; //trueX must be .5 or higher
				if (trueX>N+0.5) trueX=N+ 0.5f; //trueX must be N+.5 or lower
				trueXindex=(int)trueX; 
				trueXindexP1=trueXindex+1;
				if (trueY<0.5) trueY=0.5f; //trueY must be .5 or higher
				if (trueY>N+0.5) trueY=N+ 0.5f; //trueY must be N+.5 or lower
				trueYindex=(int)trueY; 
				trueYindexP1=trueYindex+1;
				trueXdeciPt = trueX-trueXindex; 
				trueXdeciPtSubFrom1 = 1-trueXdeciPt; 
				trueYdeciPt = trueY-trueYindex; 
				trueYdeciPtSubFrom1 = 1-trueYdeciPt;
				
				d[IX(i,j)] = 
						trueXdeciPtSubFrom1*(trueYdeciPtSubFrom1*d0[IX(trueXindex,trueYindex)]+trueYdeciPt*d0[IX(trueXindex,trueYindexP1)])+
						trueXdeciPt*(trueYdeciPtSubFrom1*d0[IX(trueXindexP1,trueYindex)]+trueYdeciPt*d0[IX(trueXindexP1,trueYindexP1)]);

			}
		}
		setBnd (b, d);
	}
	public void project(float[] Vx, float[] Vy, float[] p, float[] div) {
		int i, j, k;
		float h;
		h = (float) (1.0/N);
		for ( i=1 ; i<=N ; i++ ) {
			for ( j=1 ; j<=N ; j++ ) {
				div[IX(i,j)] = (float) (-0.5*h*(Vx[IX(i+1,j)]-Vx[IX(i-1,j)]+
						Vy[IX(i,j+1)]-Vy[IX(i,j-1)]));
				p[IX(i,j)] = 0;
			}
		}
		setBnd(0, div); 
		setBnd(0, p);
		for ( k=0 ; k<20 ; k++ ) {
			for ( i=1 ; i<=N ; i++ ) {
				for ( j=1 ; j<=N ; j++ ) {
					p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+
							p[IX(i,j-1)]+p[IX(i,j+1)])/4;
				}
			}
			setBnd( 0, p);
		}
		for ( i=1 ; i<=N ; i++ ) {
			for ( j=1 ; j<=N ; j++ ) {
				Vx[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
				Vy[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
			}
		}
		setBnd (1, Vx); 
		setBnd (2, Vy);

	}
	public void setBnd(int b, float[] x) {

	    for(int k = 1; k < N - 1; k++) {
	        for(int i = 1; i < N - 1; i++) {
	            x[IX(i, 0  )] = b == 2 ? -x[IX(i, 1  )] : x[IX(i, 1  )];
	            x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
	        }
	    }
	    for(int k = 1; k < N - 1; k++) {
	        for(int j = 1; j < N - 1; j++) {
	            x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
	            x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
	        }
	    }
	    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	    x[IX(0, N-1)] = 0.5f * (x[IX(1, N-1)] + x[IX(0, N-2)]);
	    x[IX(N-1, 0)] = 0.5f * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
	    x[IX(N-1, N-1)] = 0.5f * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);
	}
	
	public void addDye(int x, int y, float amt) {
		int num = 3;
		for(int i = x-num;i<=x+num;i++) {
			for(int j = y-num;j<=y+num;j++) {
				if(density[IX(i, j)]+amt>255) {
					density[IX(i, j)] = 255;
				}else {
					density[IX(i, j)] += amt;
				}
			}
		}
	}
	public void addVelocity(int x, int y, float amtX, float amtY) {
		int num = 3;
		for(int i = x-num;i<=x+num;i++) {
			for(int j = y-num;j<=y+num;j++) {
				if(Vx[IX(i, j)]+amtX>255) {
					Vx[IX(i, j)] = 255;
				}else {
					Vx[IX(i, j)] += amtX;
				}
				if(Vy[IX(i, j)]+amtY>255) {
					Vy[IX(i, j)] = 255;
				}else {
					Vy[IX(i, j)] += amtY;
				}
			}
		}
	}
	public void drawDens(Graphics g) {
		r+=r1;
		gr+=g1;
		b+=b1;
		if(r==0&&gr==0&&b==0) {
			r1=3;
		}else if(r==255) {
			r1=-3f;
			g1=3f;
			b1=0;
		}else if(gr==255) {
			g1=-3f;
			b1=3f;
			r1=0;
		}else if(b==255) {
			b1=-3f;
			r1=3f;
			g1=0;
		}
		for(int i = 0;i<this.N;i++) {
			for(int j = 0;j<this.N;j++) {
				int d = (int) this.density[IX(i, j)];
				
				g.setColor(new Color((int)r, (int)gr,(int) b, (int)d));
				g.fillRect(i*SCALE, j*SCALE, SCALE, SCALE);
			}
		}
		
	}
	int yeet = 0;
	public void drawVel(Graphics g) {
		for(int i = 0;i<this.N;i++) {
			for(int j = 0;j<this.N;j++) {
				int v = (Math.max((int)Math.abs(this.Vx[IX(i, j)]), (int)Math.abs(this.Vy[IX(i, j)])));
				if(v!=0) {
					v=150;
				}
				g.setColor(new Color(238,210,67, v));
				g.fillRect(i*SCALE, j*SCALE, SCALE, SCALE);
			}
		}
	}
	public int IX(int x, int y) {
		if(x+this.N*y>(1256*this.N)/10) {
			return (x+this.N*y)/((1256*this.N)/10);
		}else if(x+this.N*y<0) {
			return (x+this.N*y)+((1256*this.N)/10);
		}
		return x + this.N*y;
	}
}