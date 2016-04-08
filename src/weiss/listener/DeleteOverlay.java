package weiss.listener;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import icy.painter.Overlay;
import plugins.weiss.StructureTensor3D.StructureTensor3D;

/**
 * Create an actionListener to delete the overlay on the vtkViewer in Icy 
 * @author Guillaume
 *
 */
public class DeleteOverlay implements ActionListener {

	@Override
	public void actionPerformed(ActionEvent arg0) {
		Overlay o = StructureTensor3D.o ; 
		if(o!=null)
			StructureTensor3D.se.removeOverlay(o) ; 
		
	}

}
