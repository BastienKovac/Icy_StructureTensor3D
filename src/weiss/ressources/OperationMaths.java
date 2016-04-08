package weiss.ressources;

import java.util.LinkedList;
import java.util.List;

/**
 * This clas includes all of the function which can be use in the structure
 * Tensor algorithm All of the matrix will be arrays in 2 or 3 dimensions
 * 
 * @author Guillaume
 *
 */
public class OperationMaths {

	// TO BE THREAD
	/**
	 * Create the equivalence of the function meshgrid in matlab
	 * 
	 * @param vmin1
	 * @param pas1
	 * @param vmax1
	 * @param vmin2
	 * @param pas2
	 * @param vmax2
	 * @param vmin3
	 * @param pas3
	 * @param vmax3
	 * @return A list of 3 array in 3D (tab1 for x, tab2 for y, tab3 for z)
	 */
	public static List<double[][][]> meshgrid(double vmin1, double pas1, double vmax1, // X
			double vmin2, double pas2, double vmax2, // Y
			double vmin3, double pas3, double vmax3) { // Z

		double[][][] tab1, tab2, tab3;
		int t1, t2, t3;
		t1 = (int) (((vmax1 - vmin1) / pas1) + 1);
		t2 = (int) (((vmax2 - vmin2) / pas2) + 1);
		t3 = (int) (((vmax3 - vmin3) / pas3) + 1);
		tab1 = new double[t3][t1][t2];
		tab2 = new double[t3][t1][t2];
		tab3 = new double[t3][t1][t2];

		double vtmpI, vtmpK;

		for (int k = 0; k < t3; k++) {
			vtmpK = vmin3 + k * pas3;
			for (int i = 0; i < t1; i++) {
				vtmpI = vmin2 + i * pas2;
				for (int j = 0; j < t2; j++) {
					tab1[k][i][j] = vmin1 + j * pas1;
					tab2[k][i][j] = vtmpI;
					tab3[k][i][j] = vtmpK;
				}
			}
		}

		List<double[][][]> l = new LinkedList<double[][][]>();
		l.add(tab1);
		l.add(tab2);
		l.add(tab3);

		return l;

	}

	/**
	 * Calculate the equivalence of the function fftShit in matlab
	 * 
	 * @param d
	 * @return an array in 3D including the fftshift of the matrix d
	 */
	public static double[][][] fftShiftMat(double[][][] d) {
		double[][][] tmp = new double[d.length][d[0].length][d[0][0].length];
		int n1 = d[0].length;
		int n2 = d[0][0].length;
		int n3 = d.length;
		for (int k = 0; k < tmp.length; k++) {
			for (int i = 0; i < tmp[0].length; i++) {
				for (int j = 0; j < tmp[0][0].length; j++) {
					tmp[(k + (n3 / 2)) % n3][(i + (n1 / 2)) % n1][(j + (n2 / 2)) % n2] = d[k][i][j];
				}
			}
		}
		return tmp;
	}

	// ----------Matrix---------//
	// Operations //
	// //
	// //
	// //
	// matrix square
	/**
	 * Make the square root of a matrix in 2D
	 * 
	 * @param d
	 */
	public static void matSquare2D(double[][] d) {
		for (int i = 0; i < d.length; i++) {
			for (int j = 0; j < d[0].length; j++) {
				d[i][j] = (d[i][j] * d[i][j]);
			}
		}
	}

	/**
	 * Calculate the square of each term of the matrix Use it only when all the
	 * z dimensions are equal
	 * 
	 * @param d
	 */
	public static void matSquare3DSpecific(double[][][] d) {
		double[][] dTmp = d[0];
		matSquare2D(dTmp);
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					d[k][i][j] = dTmp[i][j];
				}
			}
		}
	}

	/**
	 * Calculate the square of each term of the matrix Use it by default (less
	 * effiency)
	 * 
	 * @param d
	 */
	public static void matSquare3D(double[][][] d) {
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					d[k][i][j] = d[k][i][j] * d[k][i][j];
				}
			}
		}
	}

	// TO BE THREAD IF LONG
	// matrix add
	/**
	 * Make the addition of two matrix 2D
	 * 
	 * @param d
	 * @param c
	 * @return a matrix in 2D
	 */
	public static double[][] matAdd2D(double[][] d, double[][] c) {
		double[][] tmp = new double[d.length][d[0].length];
		for (int i = 0; i < tmp.length; i++) {
			for (int j = 0; j < tmp[0].length; j++) {
				tmp[i][j] = (d[i][j] + c[i][j]);
			}
		}
		return tmp;
	}

	/**
	 * Add each term of the matrix d with the matrix c Use it only when the z
	 * dimensions are equal
	 * 
	 * @param d
	 * @param c
	 * @return a matrix in 3D
	 */
	public static double[][][] matAdd3DSpecific(double[][][] d, double[][][] c) {
		double[][][] tmp = new double[d.length][d[0].length][d[0][0].length];
		double[][] firstAdd = matAdd2D(d[0], c[0]);
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					tmp[k][i][j] = firstAdd[i][j];
				}
			}
		}
		return tmp;
	}

	/**
	 * Add each term of the matrix d with the matrix c Use it by default (less
	 * efficency)
	 * 
	 * @param d
	 * @param c
	 * @return a matrix in 3D
	 */
	public static double[][][] matAdd3D(double[][][] d, double[][][] c) {
		double[][][] tmp = new double[d.length][d[0].length][d[0][0].length];
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					tmp[k][i][j] = d[k][i][j] + c[k][i][j];
				}
			}
		}
		return tmp;
	}

	/**
	 * Divide each term of the matrix by a real
	 * 
	 * @param d
	 * @param r
	 */
	public static void divideByReal(double[][][] d, double r) {
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					d[k][i][j] = d[k][i][j] / r;
				}
			}
		}
	}

	/**
	 * Make the exponential function of each term of the matrix
	 * 
	 * @param d
	 */
	public static void matExp(double[][][] d) {
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					d[k][i][j] = Math.exp(d[k][i][j]);
				}
			}
		}
	}

	/**
	 * Make the sum of all of the term of the matrix
	 * 
	 * @param d
	 * @return a real (sum of the term)
	 */
	public static double matSum(double[][][] d) {
		double tmp = 0;
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					tmp += d[k][i][j];
				}
			}
		}
		return tmp;
	}

	/**
	 * Multiply each term of the matrix by a real
	 * 
	 * @param d
	 * @param r
	 */
	public static void matMultByReal(double[][][] d, double r) {
		for (int k = 0; k < d.length; k++) {
			for (int i = 0; i < d[0].length; i++) {
				for (int j = 0; j < d[0][0].length; j++) {
					d[k][i][j] *= r;
				}
			}
		}
	}

	// **********PRINT**********//
	//
	//
	// *************************//
	/**
	 * Print the value of the matrix (3D) to debug
	 * 
	 * @param tab
	 * @return a string
	 */
	public static String printArray3D(double[][][] tab) {
		String s = "";
		for (int k = 0; k < tab.length; k++) {
			for (int i = 0; i < tab[0].length; i++) {
				for (int j = 0; j < tab[0][0].length; j++) {
					s += "[" + tab[k][i][j] + "] ";
				}
				s += "\n";
			}
			s += "\n";
		}

		return s;
	}

	/**
	 * Print the value of the matrix (2D) to debug
	 * 
	 * @param tab
	 * @return a string
	 */
	public static String printArray2D(double[][] tab) {
		String s = "";
		for (int i = 0; i < tab.length; i++) {
			for (int j = 0; j < tab[0].length; j++) {
				s += "[" + tab[i][j] + "] ";
			}
			s += "\n";
		}

		return s;
	}

	/**
	 * Print the command line to write in matlab to have the same matrix
	 * 
	 * @param d
	 * @return a string
	 */
	public static String debugMatLab(double[][][] d) {
		String m = "";
		int nk = d.length, ni = d[0].length, nj = d[0][0].length;

		for (int k = 0; k < nk; k++) {
			m += "[";
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					m += d[k][i][j];
					if (j != nj - 1)
						m += ",";
				}
				if (i != ni - 1)
					m += ";";
			}
			if (k != nk - 1)
				m += "],";
			else
				m += "]";
		}
		m += "]  ; u=reshape(u," + ni + "," + nj + "," + nk + ")";

		return m;
	}

}
