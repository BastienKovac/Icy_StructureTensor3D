package weiss.listener;

import java.awt.BorderLayout;
import java.awt.Desktop;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import org.jdesktop.swingx.JXHyperlink;

import icy.gui.frame.IcyFrame;




public class OpenAbout implements ActionListener{

	public static final String url = "http://www.math.univ-toulouse.fr/~weiss/" ; 
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		IcyFrame j = new IcyFrame() ; 
		j.setVisible(true);
		j.externalize();
		j.setSize(600,300);
		FlowLayout b = new FlowLayout() ;  

		JPanel ip = new JPanel() ; 
		ip.setLayout(b);
		
		JLabel text = new JLabel("<html><br><br>In case you use this algorithm, please cite  :"
				+ "<br> Structure tensor based analysis of cells and nuclei organization in tissues."
				+ "<br>W. Zhang, J. Fehrenbach, A. Desmaison, V. Lobjois,"
				+ "<br> B. Ducommun, P. Weiss, IEEE TMI (2015) "
				+ "<br><br>Further information there:<br></html>",SwingConstants.CENTER) ;
		Font f = new Font("test", 0, 14) ;  
		text.setFont(f) ; 
		ip.add(text,BorderLayout.CENTER);
		JXHyperlink link = new JXHyperlink(new SampleAction()) ; 
		link.setFont(f);
		ip.add(link) ; 
		
		j.add(ip);
		
		
		
	}
	
	private static void open(URI uri) {
		if (Desktop.isDesktopSupported()) {
			try {
				Desktop.getDesktop().browse(uri);
			}
			catch(IOException e) {
				//ERROR
			}
		}
		else {
			//ERROR
		}
	}
	
	private static class SampleAction extends AbstractAction {
		private static final long serialVersionUID = 1L;

		public SampleAction() {
			super.putValue(Action.NAME, "http://www.math.univ-toulouse.fr/~weiss/PagePublications.html");
			super.putValue(Action.SHORT_DESCRIPTION, "http://www.math.univ-toulouse.fr/~weiss/PagePublications.html");
		}

		@Override
		public void actionPerformed(ActionEvent arg0) {
			URI uri;
			try {
				uri = new URI(OpenAbout.url);
				open(uri);
			} catch (URISyntaxException e1) {
				e1.printStackTrace();
			} 
		}
	}
	

}
