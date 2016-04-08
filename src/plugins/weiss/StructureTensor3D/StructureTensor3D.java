package plugins.weiss.StructureTensor3D;

import java.util.LinkedList;
import java.util.List;

import weiss.listener.DeleteOverlay;
import weiss.listener.OpenAbout;
import weiss.ressources.Double3DReal;
import weiss.structureTensor3D.StructureTensorFunction3D;
import icy.painter.Overlay;
import icy.sequence.Sequence;
import icy.type.collection.array.Array2DUtil;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;

/**
 * <p>
 * Apply the structure tensor algorithm in 3D of the PRIMO team of the ITAV
 * (Toulouse, FRANCE) This is a plugin for the software Icy To use it, you have
 * to get a 8-bits image in 3D You have to fix the parameters sigma and delta Do
 * not hesitate to do some different test with different parameters
 * </p>
 * <p>
 * In case you use this algorithm, please cite : WEISS Pierre, FEHRENBACH
 * Jerome, DE BRITO Guillaume - ITAV (Toulouse, France)
 * </p>
 * 
 * @see <a href=
 *      "http://www.math.univ-toulouse.fr/~weiss/Publis/Journals/2014/Structure_Tensor_Cell_Organization_Zhang_Weiss_2014.pdf">
 *      http://www.math.univ-toulouse.fr/~weiss/Publis/Journals/2014/
 *      Structure_Tensor_Cell_Organization_Zhang_Weiss_2014.pdf</a>
 * @author DE BRITO Guillaume with the help of WEISS Pierre and FEHRENBACH
 *         Jer�me
 *
 */
public class StructureTensor3D extends EzPlug implements Block {

	private EzVarSequence varCurrentSeq = new EzVarSequence("Sequence");
	private EzVarSequence varThresholderSeq = new EzVarSequence("Threshold");
	private EzVarInteger varSigma = new EzVarInteger("Sigma :");
	private EzVarInteger varDelta = new EzVarInteger("Delta : ");
	private EzVarInteger varResolution = new EzVarInteger("Resolution of the ellipsoids : ");
	private EzButton button;
	private EzButton buttonAbout;
	private int width, height, depth;

	public static Sequence se;
	public static Sequence threshold;
	public static Overlay o;
	public static EzVarBoolean randomizeCenterOfEllipsis;

	@Override
	protected void initialize() {
		varCurrentSeq.setToolTipText("Select the sequence to use");
		super.addEzComponent(varCurrentSeq);

		varThresholderSeq.setToolTipText("Select the sequence including the threshold to use");
		super.addEzComponent(varThresholderSeq);

		varSigma.setToolTipText("Set the value of sigma, characteristic size of the ellipses to analyze");
		varSigma.setValue(15);
		super.addEzComponent(varSigma);

		varDelta.setToolTipText("Set the value of delta, step size between segments/ellipses");
		varDelta.setValue(25);
		super.addEzComponent(varDelta);

		varResolution.setToolTipText("Set the value of the resolution, reduce it to gain more performances ");
		varResolution.setValue(15);
		super.addEzComponent(varResolution);

		randomizeCenterOfEllipsis = new EzVarBoolean("Randomize", true);
		randomizeCenterOfEllipsis
				.setToolTipText("Randomize the position of the ellipsises within the blocks determined with Delta");
		super.addEzComponent(randomizeCenterOfEllipsis);

		button = new EzButton("Delete overlay", new DeleteOverlay());
		button.setVisible(true);
		super.addEzComponent(button);

		buttonAbout = new EzButton("About", new OpenAbout());
		super.addEzComponent(buttonAbout);

		o = null;
		se = null;
		threshold = null;

	}

	@Override
	protected void execute() {
		se = varCurrentSeq.getValue();
		threshold = varThresholderSeq.getValue();
		int sigma = varSigma.getValue();
		int delta = varDelta.getValue();
		int resolution = varResolution.getValue();
		width = se.getWidth();
		height = se.getHeight();
		depth = se.getSizeZ();

		List<double[][][]> l = getVoxelsValues();
		double[][][] voxel = l.get(0), voxelTresh = l.get(1);

		Double3DReal u = new Double3DReal(voxel);
		o = StructureTensorFunction3D.structureTensor(u, voxelTresh, sigma, delta, resolution);
		se.addOverlay(o);

		se.overlayChanged(o);

	}

	@Override
	public void clean() {
		o = null;
		se = null;
	}

	/**
	 * Get the voxels values of an image
	 * 
	 * @return voxel array in 3D
	 */
	private List<double[][][]> getVoxelsValues() {

		double[][] tmpVoxel = Array2DUtil.arrayToDoubleArray(se.getDataXYZ(0, 0), se.isSignedDataType());
		double[][] tmpThreshVoxel = Array2DUtil.arrayToDoubleArray(threshold.getDataXYZ(0, 0),
				threshold.isSignedDataType());
		double[][][] voxel = new double[depth][height][width];
		double[][][] threshVoxel = new double[depth][height][width];
		for (int k = 0; k < depth; k++) {
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					voxel[k][i][j] = tmpVoxel[k][i * width + j];
					threshVoxel[k][i][j] = tmpThreshVoxel[k][i * width + j];
				}
			}
		}
		List<double[][][]> l = new LinkedList<double[][][]>();
		l.add(voxel);
		l.add(threshVoxel);
		return l;
	}

	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add("Sequence", varCurrentSeq.getVariable());
		inputMap.add("Threshold", varThresholderSeq.getVariable());
		inputMap.add("Var sigma", varSigma.getVariable());
		inputMap.add("Var delta", varDelta.getVariable());
		inputMap.add("Var resolution", varResolution.getVariable());

	}

	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add("Sequence", varCurrentSeq.getVariable());
	}

}
