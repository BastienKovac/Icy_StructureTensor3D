package weiss.ressources;

import weiss.fourier.res.DoubleFFT_3D;

/**
 * Create a 3D array with complex to manipulate the FFT more easily
 * @author Guillaume
 *
 */
public class Double3DComplex {
	
	private double[][][] array ; 
	private int depth , rows, columns ;
	private DoubleFFT_3D fft ; 
	
	/**
	 * Create a Double3DComplex with the dimensions d , r , c 
	 * @param d
	 * @param r
	 * @param c
	 */
	public Double3DComplex(int d , int r , int c) {
		depth = d ; 
		rows = r ; 
		columns = c ; 
		array = new double[d][r][c*2] ; 
	}
	
	/**
	 * Create a Double3DComplex with another array 3D 
	 * @param tab
	 */
	public Double3DComplex(double[][][] tab ) {
		depth = tab.length ; 
		rows = tab[0].length ; 
		columns = tab[0][0].length/2 ; 
		array = new double[depth][rows][columns*2] ; 
		for (int k=0 ; k<depth ; k++) {
			for (int i=0 ; i<rows ; i++) {
				for (int j=0 ; j<columns*2 ; j++) {
					array[k][i][j]=tab[k][i][j] ; 
				}
			}
		}
		fft = new DoubleFFT_3D(depth, rows, columns) ; 
		
	}
	
	/**
	 * Get the array of the Double3DComplex
	 * @return
	 */
	public double[][][] getArray() {
		return array ; 
	}
	
	/**
	 * Get the value at the coordinates k , i  j
	 * @param k
	 * @param i
	 * @param j
	 * @param complex
	 * @return a real
	 * @throws IllegalArgumentException
	 */
	public double getValue(int k , int i, int j, boolean complex) throws IllegalArgumentException{
		
		if(k < depth && i < rows && j < columns)
		{
			if(complex)
				return array[k][i][2*j+1];
			else 
				return array[k][i][2*j];
		} else {
			throw new IllegalArgumentException("Row and colomuns wrong") ; 
		}
	}
	
	/**
	 * Set the value at the coordinates k, i , j  
	 * @param value
	 * @param k
	 * @param i
	 * @param j
	 * @param complex
	 * @throws IllegalArgumentException
	 */
	public void setValue(double value, int k, int i, int j, boolean complex) throws IllegalArgumentException{
		
		if(k < depth && i < rows && j < columns)
		{
			if(complex)
				array[k][i][2*j+1] = value;
			else
				array[k][i][2*j] = value;
		} else {
			throw new IllegalArgumentException("rows and columns wrong") ; 
		}
		
	}
	
	
	
	
	
	
	
	
	
	//*************FOURIER*************//
	//								
	//
	//*********************************//
	/**
	 * initialize the FFT 
	 */
	public void initFFT() {
		if (fft==null) 
			fft = new DoubleFFT_3D(depth, rows, columns) ; 
	}
	
	/**
	 * Transform complex in real
	 * @return tmp
	 */
	public double [][][] transformInReal () {
		double [][][] tmp = new double[depth][rows][columns] ; 
		for (int k=0 ; k<array.length ; k++) {
			for (int i=0 ; i<array[0].length ; i++) {
				for (int j=0 ; j<(array[0][0].length)/2 ; j++) {
					tmp[k][i][j] = array[k][i][j*2] ; 
				}
			}
		}
		return tmp ; 
	}
	
	/** 
	 * Calculate FFT 
	 * @return res
	 */
	public Double3DComplex getFFTn()
	{
		initFFT();
		double[][][] backup = new double[depth][rows][columns*2];
		for (int k=0 ; k<depth ; k++) {
			for(int i= 0; i < rows; i++) {
				for(int j= 0; j< columns*2; j++) { 
					backup[k][i][j]= array[k][i][j];
				}
			}
		}
		fft.complexForward(this.array);	// call the extern library
		Double3DComplex res = new Double3DComplex(this.array);
		array=backup ; 
		return res;
	}
	
	
	
	/**
	 * Calculate the inverse of FFT 
	 * @param scaling
	 * @return res
	 */
	public Double3DComplex getIFFTn(boolean scaling) {
		this.initFFT();
		// lets backup array
		double[][][] backup = new double[depth][rows][columns*2];
		for (int k=0 ; k<depth ; k++) {
			for(int i= 0; i < rows; i++) {
				for(int j= 0; j< columns*2; j++) { 
					backup[k][i][j]= array[k][i][j];
				}
			}
		}
		fft.complexInverse(this.array, scaling);  // call the extern library
		Double3DComplex res = new Double3DComplex(this.array);
		array=backup ; 
		return res;
	}
	
	
	
	
	//**************Operations**************//
	//								
	//			
	//
	//**************************************//
	
	/**
	 * Multiplication of two complex matrix
	 * @param d
	 * @return Double3DCompelex
	 */
	public Double3DComplex multiply(Double3DComplex d)
	{
		Double3DComplex res = new Double3DComplex(depth,rows,columns);
		double[] val = new double[2];
		for (int k=0 ; k<depth ; k++) {
			for(int i= 0; i< rows; i++) { 
				for(int j= 0; j < columns; j++) {
					val = multiplyComplex(this.getValue(k, i, j, false), this.getValue(k, i, j, true), d.getValue(k, i, j, false), d.getValue(k, i, j, true));
					res.setValue(val[0],k, i, j, false);
					res.setValue(val[1],k, i, j, true);
				}
			}
		}
		return res;
	}
	
	

	/**
	 * Multiplication of two complex
	 * @param real1
	 * @param im1
	 * @param real2
	 * @param im2
	 * @return an array in 1D (including the real and imaginary)
	 */
	public static double[] multiplyComplex(double real1, double im1, double real2, double im2) {
		double r = real1 * real2 - im1 * im2;
		double imag = real1 * im2 + im1 * real2;
		double[] res =  {r, imag};
		return res;
	}
	
	
	
	
	
	
	
	
	

}
