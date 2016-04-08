package weiss.ressources;



/**
 * Create an object with a 3D array to manipulate it. You can use the gradient and the FFT
 * @author Guillaume
 *
 */
public class Double3DReal {

	private double [][][] array  ; 
	private int rows, columns, depth ; 
	private boolean complexExist ;
	private Double3DComplex complex ;  
	
	/**
	 * Construct the Double3DReal with the dimensions 
	 * All of the values of the terms will be zero
	 * @param d
	 * @param r
	 * @param c
	 */
	public Double3DReal(int d , int r, int c ) {
		rows=r ; 
		columns=c ;
		depth=d ; 
		array = new double[depth][rows][columns] ;
		complexExist=false ; 
	}
	
	/**
	 * Construct the Double3DReal with an array in 3D
	 * @param tab
	 */
	public Double3DReal(double[][][] tab) {
		depth=tab.length ;
		rows=tab[0].length ; 
		columns=tab[0][0].length ;
		array = new double[depth][rows][columns] ;
		
		for (int k=0 ; k<depth ; k++) {
			for (int i=0 ; i<rows ; i++) {
				for (int j=0 ; j<columns ; j++) {
					array[k][i][j]=tab[k][i][j] ; 
				}
			}
		}
		complexExist=false ; 
	}
	
	/**
	 * Set the value at the coordinates d , r  , c 
	 * @param d
	 * @param r
	 * @param c
	 * @param value
	 */
	public void setValue(int d, int r, int c,  double value) {
		array[d][r][c]=value ; 
	}
	
	/**
	 * Get the value at the coordinates d, r , c
	 * @param d
	 * @param r
	 * @param c
	 * @return a real
	 */
	public double getValue(int d, int r, int c ) {
		return array[d][r][c] ; 
	}
	
	/**
	 * Get the depth of the Double3DReal
	 * @return
	 */
	public int getSizeK () {
		return depth ; 
	}
	
	/**
	 * get the number of rows of the Double3DReal
	 * @return
	 */
	public int getSizeI () {
		return rows ; 
	}
	
	/**
	 * Get the number of columns of the Double3DReal
	 * @return
	 */
	public int getSizeJ () {
		return columns ; 
	}
	/**
	 * Get the array of the object
	 * @return
	 */
	public double[][][] getArray () {
		return array ; 
	}
	
	
	
	
	
	//********************GRADIENT*********************//
	//
	//
	//*************************************************//
	/**
	 * Create the filter h1 to calculate the gradient
	 * @return an array in 2D
	 */
	public double [][][] createH1Filter() {   
		double [][][] h1 = new double[3][3][3] ; 
		
		h1[0][0][0] =0.0153 ;   h1[0][1][0] =0.0568 ;   h1[0][2][0] =0.0153 ;
		h1[1][0][0] =0 ;	    h1[1][1][0] =0 ;        h1[1][2][0] =0 ;
		h1[2][0][0] = -0.0153 ; h1[2][1][0] = -0.0568 ; h1[2][2][0] = -0.0153 ;
		
		h1[0][0][1] =0.0568 ;   h1[0][1][1] =0.2117 ;   h1[0][2][1] =0.0568 ;
		h1[1][0][1] =0 ;        h1[1][1][1] =0 ;        h1[1][2][1] =0 ;
		h1[2][0][1] = -0.0568 ; h1[2][1][1] = -0.2117 ; h1[2][2][1] = -0.0568 ;
		
		h1[0][0][2] =0.0153 ;   h1[0][1][2] =0.0568 ;   h1[0][2][2] =0.0153 ;
		h1[1][0][2] =0 ;        h1[1][1][2] =0 ;        h1[1][2][2] =0 ;
		h1[2][0][2] = -0.0153 ; h1[2][1][2] = -0.0568 ; h1[2][2][2] = -0.0153 ;
		
		return h1 ; 
	}
	
	/**
	 * Create the filter h2 to calculate the gradient
	 * @return an array in 2D
	 */
	public double [][][] createH2Filter() { 
		double [][][] h2 = new double[3][3][3] ; 
		
		h2[0][0][0]=0.0153 ;   h2[0][0][1]=0.0568 ;   h2[0][0][2]= 0.0153 ; 
		h2[0][1][0]=0 ;        h2[0][1][1]=0 ;        h2[0][1][2]= 0 ;
		h2[0][2][0]= -0.0153 ; h2[0][2][1]= -0.0568 ; h2[0][2][2]= -0.0153 ;
		
		h2[1][0][0]=0.0568 ;   h2[1][0][1]=0.2117 ;   h2[1][0][2]= 0.0568 ; 
		h2[1][1][0]=0 ;        h2[1][1][1]=0 ;        h2[1][1][2]= 0 ; 
		h2[1][2][0]= -0.0568 ; h2[1][2][1]= -0.2117 ; h2[1][2][2]= -0.0568 ;
		
		h2[2][0][0]=0.0153   ; h2[2][0][1]=0.0568   ; h2[2][0][2]= 0.0153 ; 
		h2[2][1][0]=0        ; h2[2][1][1]=0        ; h2[2][1][2]= 0 ; 
		h2[2][2][0]= -0.0153 ; h2[2][2][1]= -0.0568 ; h2[2][2][2]= -0.0153 ; 
		
		return h2 ; 
	}
	
	/**
	 * Create the filter h3 to calculate the gradient
	 * @return an array in 2D
	 */
	public double[][][] createH3Filter() { 
		double[][][] h3 = new double[3][3][3] ; 
		
		h3[0][0][0] =0.0153 ; h3[0][0][1] =0 ; h3[0][0][2] = -0.0153 ;
		h3[1][0][0] =0.0568 ; h3[1][0][1] =0 ; h3[1][0][2] = -0.0568 ;
		h3[2][0][0] =0.0153 ; h3[2][0][1] =0 ; h3[2][0][2] = -0.0153 ;
		
		h3[0][1][0] =0.0568 ; h3[0][1][1] =0 ; h3[0][1][2] = -0.0568 ;
		h3[1][1][0] =0.2117 ; h3[1][1][1] =0 ; h3[1][1][2] = -0.2117 ;
		h3[2][1][0] =0.0568 ; h3[2][1][1] =0 ; h3[2][1][2] = -0.0568 ;
		
		h3[0][2][0] =0.0153 ; h3[0][2][1] =0 ; h3[0][2][2] = -0.0153 ;
		h3[1][2][0] =0.0568 ; h3[1][2][1] =0 ; h3[1][2][2] = -0.0568 ;
		h3[2][2][0] =0.0153 ; h3[2][2][1] =0 ; h3[2][2][2] = -0.0153 ;
		
		return h3 ; 
		
	}
	
	/**
	 * Calculate the gradient in 3D
	 * @return
	 */
	public double[][][][] calculGradient() {
		double [][][][] gradient = new double[3][depth] [rows] [columns]  ; 
		double d1=0 , d2=0 , d3=0 ; 
		int zTmp,iTmp,jTmp ; 
		boolean out = false ; 
		
		double[][][] h1 = createH1Filter() ,h2 = createH2Filter() ,h3 = createH3Filter() ; 
		
		for (int k=0 ; k<depth ; k++) {
			for (int i=0 ; i<rows ; i++) {
				for (int j=0 ; j<columns ; j++) {
					
					for (int kP=0 ; kP<3 && !out; kP++) {
						for (int iP=0 ; iP<3 && !out; iP++) {
							for (int jP=0 ; jP<3 && !out; jP++) {
								zTmp=k+(kP-1) ; 
								iTmp=i+(iP-1) ; 
								jTmp=j+(jP-1) ; 
								
								if (zTmp<0 || zTmp>=depth || iTmp<0 || iTmp>=rows || jTmp<0 || jTmp>=columns) 
									out=true ; 
								else {
									d1+= ( getValue(zTmp, iTmp, jTmp) * h1[kP][iP][jP] ) ; 
									d2+= ( getValue(zTmp, iTmp, jTmp) * h2[kP][iP][jP] ) ; 
									d3+= ( getValue(zTmp, iTmp, jTmp) * h3[kP][iP][jP] ) ; 
								}
								
							}
						}
					}
					if (out) {
						gradient[0][k][i][j] = 0 ; 
						gradient[1][k][i][j] = 0 ; 
						gradient[2][k][i][j] = 0 ; 
					}
					else 
					{
						gradient[0][k][i][j] = d1 ; 
						gradient[1][k][i][j] = d2 ; 
						gradient[2][k][i][j] = d3 ; 
					}
					d1=0 ; 
					d2=0 ; 
					d3=0 ; 
					out=false ; 
					
					
				}
			}
		}
		
		
		return gradient ; 
	}
	
	
	
	
	
	
	
	
	
	
	//**************FOURIER**************//
	//
	//
	//***********************************//
	/**
	 * Initialization of the complex matrix
	 */
	public void initFourier () {
		double[][][] tmp = transformInComplex() ;
		complex = new Double3DComplex(tmp) ;
		complexExist = true ; 
	}
	
	/**
	 * Transform a real matrix in complex matrix
	 * @return tmp
	 */
	private double[][][] transformInComplex () {
		double [][][] tmp = new double[depth][rows][columns*2] ;
		for (int k=0 ; k<array.length ; k++) {
			for (int i = 0 ; i<array[0].length ; i++) {
				for (int j=0 ; j<array[0][0].length ; j++) {
					tmp[k][i][j*2]=array[k][i][j] ; 
				}
			}
		}
		return tmp ;
	}
	
	/**
	 * Calculate the FFT 
	 * @return fftn
	 */
	public Double3DComplex getFFTn () {
		if (!complexExist) 
			initFourier() ; 
		return   complex.getFFTn() ; 
	}
	
	/**
	 * Calculate the inverse of the fft
	 * @return ifftn
	 */
	public Double3DComplex getIFFTn() {
		if (!complexExist) 
			initFourier() ; 
		return   complex.getIFFTn(true) ; 
	}
	
	
	/**
	 * Calculate the convolution between u and h 
	 * @param u
	 * @param h
	 * @return Double3DReal
	 */
	public static Double3DReal convolution(Double3DReal u , Double3DReal h) {
		return new Double3DReal(  (u.getFFTn().multiply(h.getFFTn())).getIFFTn(true).transformInReal()     ) ; 
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	//**************PRINT**************//
	//
	//
	//*********************************//
	/**
	 * Print the value of the Double3DReal
	 * @return a string
	 */
	public String printArray() {
		String s="" ; 
		for (int k=0 ; k<depth ; k++) {
			for (int i=0 ; i<rows ; i++) {
				for (int j=0 ; j<columns ; j++) {
					s+="["+array[k][i][j]+"] " ; 
				}
				s+="\n" ; 
			}
			s+="\n" ; 
		}
		
		return s ; 
	}
	
	
	
	
	
	
	
	//**************TEST**************//
	//
	//
	//********************************//
	/**
	 * Create a random matrix
	 * Print a debug of the matrix (print the value and the matlab command line)
	 * @param nk
	 * @param ni
	 * @param nj
	 * @return an array in 3D 
	 */
	public static double [][][] createMatrixRand3x3x3Debug(int nk , int ni , int nj) {
		double[][][] u = new double [nk][ni][nj] ; 
		String s = "", m = "u=[" ; 
		
		for (int k=0 ; k<nk ; k++) {
			for (int i=0 ; i<ni ; i++) {
				for (int j=0 ; j<nj ; j++) {
					u[k][i][j]=(int)(Math.random()*100) ; 
				}
			}
		}
		
		//PRINT
		for (int k=0 ; k<nk ; k++) {
			m+="[" ; 
			for (int i=0; i<ni ; i++) {
				for (int j=0 ; j<nj ; j++) {
					s+="["+u[k][i][j]+"] " ; 
					m+=u[k][i][j] ; 
					if (j!=nj-1)
						m+="," ; 
				}
				s+="\n" ;
				if (i!=ni-1)
					m+=";" ; 
			}
			s+="\n" ;
			if (k!=nk-1)
				m+="]," ; 
			else
				m+="]" ; 
		}
		m+="]  ; u=reshape(u,"+ni+","+nj+","+nk+")" ; 
		
		
		System.out.println(s+"\n"+m);
		
		return u ; 
	}
	
}
