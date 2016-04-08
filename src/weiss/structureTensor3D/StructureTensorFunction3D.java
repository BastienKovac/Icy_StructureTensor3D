package weiss.structureTensor3D;

import java.util.LinkedList;
import java.util.List;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import icy.painter.Overlay;
import plugins.weiss.StructureTensor3D.StructureTensor3D;
import vtk.vtkPoints;
import weiss.ressources.Double3DReal;
import weiss.ressources.OperationMaths;
import weiss.vtkressources.EllipsoidGlyph3D;

public class StructureTensorFunction3D {

	/**
	 * Calculate the tensor algorithm
	 * 
	 * @param u
	 * @param voxelTresh
	 * @param sigma
	 * @param delta
	 * @return Overlay of ellipsoids
	 */
	public static Overlay structureTensor(Double3DReal u, double[][][] voxelTresh, int sigma, int delta,
			int resolution) {
		long t0 = System.currentTimeMillis();
		int nk = u.getSizeK(), ni = u.getSizeI(), nj = u.getSizeJ();
		double[][][][] gradient = u.calculGradient();
		double[][][] d1 = gradient[0], d2 = gradient[1], d3 = gradient[2];
		Double3DReal J11, J12, J13, J22, J23, J33;
		J11 = new Double3DReal(nk, ni, nj);
		J12 = new Double3DReal(nk, ni, nj);
		J13 = new Double3DReal(nk, ni, nj);
		J22 = new Double3DReal(nk, ni, nj);
		J23 = new Double3DReal(nk, ni, nj);
		J33 = new Double3DReal(nk, ni, nj);

		for (int k = 0; k < nk; k++) {
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					J11.setValue(k, i, j, d1[k][i][j] * d1[k][i][j]);
					J12.setValue(k, i, j, d1[k][i][j] * d2[k][i][j]);
					J13.setValue(k, i, j, d1[k][i][j] * d3[k][i][j]);
					J22.setValue(k, i, j, d2[k][i][j] * d2[k][i][j]);
					J23.setValue(k, i, j, d2[k][i][j] * d3[k][i][j]);
					J33.setValue(k, i, j, d3[k][i][j] * d3[k][i][j]);
				}
			}
		}
		d1 = null;
		d2 = null;
		d3 = null;
		gradient = null;

		long t1 = System.currentTimeMillis();
		Double3DReal h = new Double3DReal(createH(ni, nj, nk, sigma));
		long t2 = System.currentTimeMillis();
		System.out.println("Time : " + (t2 - t1));

		Double3DReal J11d, J12d, J13d, J22d, J23d, J33d;
		J11d = Double3DReal.convolution(J11, h);
		J11 = null;
		J12d = Double3DReal.convolution(J12, h);
		J12 = null;
		J13d = Double3DReal.convolution(J13, h);
		J13 = null;
		J22d = Double3DReal.convolution(J22, h);
		J22 = null;
		J23d = Double3DReal.convolution(J23, h);
		J23 = null;
		J33d = Double3DReal.convolution(J33, h);
		J33 = null;
		h = null;

		double[][] tabVal; // array tmp to create SR

		// array of the SR matrix
		List<Matrix> tabSR = new LinkedList<Matrix>();

		// VTK
		vtkPoints points = new vtkPoints();
		int nbTuples = 0;

		double[][] st = new double[3][3];
		double gamma = 0.7;

		// Create the ellipsoid
		for (int k = 0; k < nk; k += delta) { // +=delta
			for (int i = 0; i < ni; i += delta) {
				for (int j = 0; j < nj; j += delta) {

					if (voxelTresh[k][i][j] > 0) {

						double pointJ = j;
						double pointI = i;
						double pointK = k;

						if (StructureTensor3D.randomizeCenterOfEllipsis.getValue()) {
							// add the center to the points (randomized)
							pointJ = j + (Math.random() - 0.5) * 2 * delta * (1 - gamma);
							pointI = i + (Math.random() - 0.5) * 2 * delta * (1 - gamma);
							pointK = k + (Math.random() - 0.5) * 2 * delta * (1 - gamma);
						}

						points.InsertNextPoint(pointJ, pointI, pointK);

						tabVal = new double[3][3];
						tabVal[0][0] = J11d.getValue(k, i, j);
						tabVal[1][0] = J12d.getValue(k, i, j);
						tabVal[2][0] = J13d.getValue(k, i, j);
						tabVal[0][1] = J12d.getValue(k, i, j);
						tabVal[1][1] = J22d.getValue(k, i, j);
						tabVal[2][1] = J23d.getValue(k, i, j);
						tabVal[0][2] = J13d.getValue(k, i, j);
						tabVal[1][2] = J23d.getValue(k, i, j);
						tabVal[2][2] = J33d.getValue(k, i, j);

						// create S and R then SR (SU here)
						Matrix tmp = new Matrix(tabVal);
						SingularValueDecomposition svd = tmp.svd();
						Matrix S = svd.getS();
						Matrix U = svd.getU();

						double[][] test = U.getArray(), test2 = new double[3][3];
						for (int iT = 0; iT < test.length; iT++) {
							for (int jT = 0; jT < test[0].length; jT++) {
								test2[(2 * iT + 2) % 3][jT] = test[iT][jT];
							}
						}
						U = new Matrix(test2);

						// normalize S
						st[0][0] = Math.sqrt(S.get(2, 2) / S.get(0, 0)) * delta * gamma;
						st[1][1] = Math.sqrt(S.get(2, 2) / S.get(1, 1)) * delta * gamma;
						st[2][2] = Math.sqrt(S.get(2, 2) / S.get(2, 2)) * delta * gamma;
						S = new Matrix(st);

						tabSR.add(S.times(U.transpose()));

						nbTuples++;

					}

				}
			}
		}
		long t0b = System.currentTimeMillis();
		System.out.println("Time function : " + (t0b - t0));
		return new EllipsoidGlyph3D(points, nbTuples, tabSR, resolution);

	}

	/**
	 * Create the gaussian kernel
	 * 
	 * @param nX
	 * @param nY
	 * @param nZ
	 * @param sigma
	 * @return gaussian kernel in a 3D array
	 */
	public static double[][][] createH(int nX, int nY, int nZ, int sigma) {
		// TODO : Math.floor() ?
		List<double[][][]> l = OperationMaths.meshgrid(-nX / 2 + 1, 1, nX / 2, -nY / 2 + 1, 1, nY / 2, -nZ / 2 + 1, 1,
				nZ / 2);
		double[][][] X = l.get(0), Y = l.get(1), Z = l.get(2);
		OperationMaths.matSquare3DSpecific(X);
		OperationMaths.matSquare3DSpecific(Y);
		OperationMaths.matSquare3D(Z);
		double[][][] tmp = OperationMaths.matAdd3DSpecific(X, Y);
		X = null;
		Y = null;
		double[][][] h = OperationMaths.matAdd3D(tmp, Z);
		tmp = null;
		Z = null;

		OperationMaths.matMultByReal(h, -1);

		OperationMaths.divideByReal(h, (2 * sigma * sigma));
		OperationMaths.matExp(h);

		OperationMaths.divideByReal(h, OperationMaths.matSum(h));
		double[][][] hFinal = OperationMaths.fftShiftMat(h);

		return hFinal;
	}

}
