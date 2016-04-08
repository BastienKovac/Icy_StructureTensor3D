package weiss.vtkressources;

import java.util.List;

import Jama.Matrix;
import icy.painter.Overlay;
import icy.painter.VtkPainter;
import plugins.weiss.StructureTensor3D.StructureTensor3D;
import vtk.vtkActor;
import vtk.vtkDoubleArray;
import vtk.vtkPoints;
import vtk.vtkPolyData;
import vtk.vtkPolyDataMapper;
import vtk.vtkProp;
import vtk.vtkSphereSource;
import vtk.vtkTensorGlyph;

/**
 * Create a vtkActor to represent some ellipsoids on the vtk Viewer on Icy You
 * have to create the object EllipsoidGlyph3D with the points of the center of
 * the ellipsoids, the number of the ellipsoids and the matrix of paramters (SR)
 * of each ellipsoid (put it on an array of Matrix)
 * 
 * @author Guillaume
 *
 */
public class EllipsoidGlyph3D extends Overlay implements VtkPainter {

	private vtkActor actor;
	private vtkPoints points;
	private int nbTuples;
	private int resolution;
	private List<Matrix> tabSR;

	public EllipsoidGlyph3D(vtkPoints points, int nbTuples, List<Matrix> tabSR, int resolution) {
		super("Ellipsoid");
		this.points = points;
		this.nbTuples = nbTuples;
		this.tabSR = tabSR;
		this.resolution = resolution;
		init();
	}

	/**
	 * Initialize the vtk actor to create the ellipsoids
	 */
	private void init() {

		vtkPolyData polyData = new vtkPolyData();
		polyData.SetPoints(points);

		vtkDoubleArray tensors = new vtkDoubleArray();

		tensors.SetNumberOfTuples(nbTuples);
		tensors.SetNumberOfComponents(9);

		double[][] tabTmp;
		for (int i = 0; i < tabSR.size(); i++) {
			tabTmp = tabSR.get(i).getArray();
			tensors.InsertTuple9(i, tabTmp[0][0], tabTmp[0][1], tabTmp[0][2], tabTmp[1][0], tabTmp[1][1], tabTmp[1][2],
					tabTmp[2][0], tabTmp[2][1], tabTmp[2][2]);
		}

		// apply tensor to the polydata
		polyData.GetPointData().SetTensors(tensors);

		// create the sphere
		vtkSphereSource sphereSource = new vtkSphereSource();
		sphereSource.SetPhiResolution(resolution);
		sphereSource.SetThetaResolution(resolution);
		sphereSource.Update();

		// create the tensor
		vtkTensorGlyph tensorGlyph = new vtkTensorGlyph();

		// set the Data
		tensorGlyph.SetInputData(polyData);

		// connect the tensorGlyph with the sphere
		tensorGlyph.SetSourceConnection(sphereSource.GetOutputPort());

		// define the tensor
		tensorGlyph.ColorGlyphsOff();
		tensorGlyph.ThreeGlyphsOff();
		tensorGlyph.ExtractEigenvaluesOff();
		tensorGlyph.Update();

		// mapper
		vtkPolyDataMapper mapper = new vtkPolyDataMapper();
		mapper.SetInputData(tensorGlyph.GetOutput());

		// actor
		actor = new vtkActor();
		actor.SetMapper(mapper);
		actor.SetScale(StructureTensor3D.se.getPixelSizeX(), StructureTensor3D.se.getPixelSizeY(),
				StructureTensor3D.se.getPixelSizeZ());
		actor.GetProperty().SetColor(1, 0, 0);

		// actor.GetProperty().SetOpacity(0.7);

	}

	@Override
	public vtkProp[] getProps() {
		return new vtkProp[] { actor };
	}

}
