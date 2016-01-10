package sprouting;

/**
 * Sprout segmentation plugin for Fiji/ImageJ
 * 
 * @version 0.2.0
 * 
 * (C) 2012-2013 Jan Eglinger
 * Heinrich-Heine University Düsseldorf
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;

import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;

import ij.measure.Calibration;
import ij.measure.ResultsTable;

import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.PlugInFilterRunner;
import ij.plugin.frame.Recorder;

import ij.process.AutoThresholder;
import ij.process.FloodFiller;
import ij.process.ImageProcessor;
import ij.process.LUT;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Font;
import java.awt.Label;

import skeleton_analysis.AnalyzeSkeleton_;
import skeleton_analysis.SkeletonResult;

/**
 * Sprout segmentation plugin
 * 
 * @author Jan Eglinger (jan.eglinger at gmail.com)
 */
public class Sprout_Analyzer implements ExtendedPlugInFilter, DialogListener {
	/*
	 * Private variables
	 */
	private final boolean PERICYTES = false;
	private final String PREF_KEY = "sprout_analyzer.";
	private ImagePlus imp;
	private int flags = DOES_8G | DOES_16 | DOES_32 | NO_CHANGES;
	private int nPasses = 0;
	private int dialog = 0;
	private Label messageArea;
	private ImagePlus bead_imp;

	/*  Image-dependent variables */
	private double pixel_size;
	private Calibration cal;

	/* 
	 * Parameters
	 */

	/* Channel configuration */
	private int ch_bead, ch_nuc, ch_sprout;
	private boolean use_bead_mask = false, use_sprout_mask = false, use_nuc_mask = false;

	/* Bead recognition */
	private String thr_bead;
	private double blur_bead, bead_radius, bead_radius_multiplier;

	/* Sprout recognition */
	private String thr_sprout;
	private double blur_sprout, min_plexus_area, min_sprout_area, max_hole_area;
	private boolean do_recover = false;
	private boolean do_exclude_borders = false;

	/* Nuclei */
	private double min_nuc_area;

	/* Pericytes */
	private int ch_endo, ch_peri;

	/* Output */
	private boolean do_n_beads, do_n_sprouts, do_sprout_area, do_total_length, do_n_cells;
	private boolean do_avg_length, do_avg_width, do_cell_density;
	private boolean do_pericyte;

	/*  Results   */
	private int num_beads, num_sprouts, num_nuc, num_peri;
	private double sprout_area, avg_sprout_length, avg_sprout_width, cell_density, totalLength, peri_area;

	/* Overridden functions */

	/**
	 * setup
	 */
	@Override
	public int setup (String arg, ImagePlus imp) {
		this.imp = imp;
		cal = imp.getCalibration();
		pixel_size = cal.getX(1.0);
		return flags;
	}

	/**
	 * setNPasses
	 */
	@Override
	public void setNPasses (int nPasses) {
		this.nPasses = nPasses;
	}

	/**
	 * showDialog
	 */
	@Override
	public int showDialog (ImagePlus imp, String cmd, PlugInFilterRunner pfr) {
		Font bold = new Font("", Font.BOLD, 12);

		int nChannels = imp.getNChannels();
		String[] channels = new String[nChannels];
		for (int i=0; i<nChannels; i++) {
			channels[i] = "Channel " + (i+1);
		}

		// read prefs
		readPrefs(nChannels);
		
		// Dialog #1
		GenericDialog gd1 = new GenericDialog("Sprout Analyzer (1/4)");
		gd1.setInsets(5, 0, 10);
		gd1.addMessage("Please choose the image channels\ncontaining the required information", bold);
		gd1.setInsets(5, 10, 0);
		gd1.addChoice("Beads", channels, channels[ch_bead - 1]);
		gd1.setInsets(0, 60, 15);
		gd1.addCheckbox("Use_bead_mask as specified", false);
		gd1.setInsets(0, 10, 0);
		gd1.addChoice("Sprouts", channels, channels[ch_sprout - 1]);
		gd1.setInsets(0, 60, 15);
		gd1.addCheckbox("Use_sprout_mask as specified", false);
		gd1.setInsets(0, 10, 0);
		gd1.addChoice("Nuclei", channels, channels[ch_nuc - 1]);
		gd1.setInsets(0, 60, 15);
		gd1.addCheckbox("Use_nucleus_mask as specified", false);
		gd1.addDialogListener(this);

		dialog = 1;
		gd1.showDialog();

		if (gd1.wasCanceled())
			return DONE;

		if (!use_bead_mask) {
			// Dialog #2
			GenericDialog gd2 = new GenericDialog("Sprout Analyzer (2/4)");
			//gd2.setInsets(0, 0, 0);
			//gd2.addMessage("Bead detection______________________________________________________________", bold);
			gd2.addChoice("Bead_threshold method", AutoThresholder.getMethods(), thr_bead);
			gd2.addSlider("Blur_radius_for_bead detection (" + cal.getUnits() + "):", 0.05, 5.0, blur_bead); // slider granularity = 0.05 only if max-min <= 5.0
			gd2.addSlider("Minimum_bead_radius (" + cal.getUnits() + "):", 0, 150, bead_radius);
			gd2.addSlider("Dilate_beads by factor:", 1, 5.5, bead_radius_multiplier);
			gd2.addPreviewCheckbox(pfr, "Preview bead detection");
			gd2.addMessage(" ");
			messageArea = (Label)gd2.getMessage();
			gd2.addDialogListener(this);

			dialog = 2;
			gd2.showDialog();

			if (gd2.wasCanceled())
				return DONE;
		}

		if (!use_sprout_mask) {
			// Dialog #3
			GenericDialog gd3 = new GenericDialog("Sprout Analyzer (3/4)");
			//gd3.setInsets(0, 0, 0);
			//gd3.addMessage("Sprout detection____________________________________________________________", bold);
			gd3.addChoice("Sprout_threshold method", AutoThresholder.getMethods(), thr_sprout);
			gd3.addSlider("Blur_radius_for_sprout detection (" + cal.getUnits() + "):", 0.05, 5.0, blur_sprout);
			// TODO add do_recover and do_exclude_borders options
			gd3.addPreviewCheckbox(pfr, "Preview sprout detection");
			gd3.addDialogListener(this);

			dialog = 3;
			gd3.showDialog();

			if (gd3.wasCanceled())
				return DONE;
		}

		// Dialog #4
		GenericDialog gd4 = new GenericDialog("Sprout Analyzer (4/4)");
		gd4.setInsets(5, 0, 10);
		gd4.addMessage("Please choose the parameters to be quantified", bold);
		String[] labels = {
			"Number_of_beads",		"Average_sprout_length",
			"Number_of_sprouts",		"Average_sprout_width",
			"Total_sprout_area",		"Cell_density",
			"Total_network_length",		"Number_of_cells"
		};
		boolean[] defaults = {
			do_n_beads,			do_avg_length,
			do_n_sprouts,			do_avg_width,
			do_sprout_area,			do_cell_density,
			do_total_length,		do_n_cells
		};
		//gd4.setInsets(10, 50, 0);
		gd4.addCheckboxGroup(4, 2, labels, defaults);
		gd4.addDialogListener(this);

		dialog = 4;
		gd4.showDialog();

		if (gd4.wasCanceled())
			return DONE;


		// write prefs
		writePrefs();

		dialog = 0;
		return flags;
	}

	/**
	 * dialogItemChanged
	 */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		switch (dialog) {
			case 1:	// read parameters from dialog 1
				IJ.log("ItemChanged 1");
				ch_bead = gd.getNextChoiceIndex() + 1;
				ch_sprout = gd.getNextChoiceIndex() + 1;
				ch_nuc = gd.getNextChoiceIndex() + 1;
				use_bead_mask = gd.getNextBoolean();
				use_sprout_mask = gd.getNextBoolean();
				use_nuc_mask = gd.getNextBoolean();
				break;
			case 2: // read parameters from dialog 2
				IJ.log("ItemChanged 2");
				thr_bead = gd.getNextChoice();
				blur_bead = gd.getNextNumber();
				bead_radius = gd.getNextNumber();
				bead_radius_multiplier = gd.getNextNumber();
				if (!gd.getPreviewCheckbox().getState()) {
					messageArea.setText(""); // clear "x beads found" message
					imp.setOverlay(null);
				}
				break;
			case 3: // read parameters from dialog 3
				IJ.log("ItemChanged 3");
				thr_sprout = gd.getNextChoice();
				blur_sprout = gd.getNextNumber();
				break;
			case 4: // read parameters from dialog 4
				IJ.log("ItemChanged 4");
				do_n_beads = gd.getNextBoolean();
				do_avg_length = gd.getNextBoolean();
				do_n_sprouts = gd.getNextBoolean();
				do_avg_width = gd.getNextBoolean();
				do_sprout_area = gd.getNextBoolean();
				do_cell_density = gd.getNextBoolean();
				do_total_length = gd.getNextBoolean();
				// do_pericyte = gd.getNextBoolean();
				do_n_cells = gd.getNextBoolean();
				break;
		}
		return (!gd.invalidNumber());
	}

	/**
	 * run
	 */
	@Override
	public void run (ImageProcessor ip) {
		if (dialog == 2) { // bead preview
			bead_imp = findBeads(imp, ch_bead); // find beads (takes time)
			ImageProcessor bead_ip = bead_imp.getProcessor();
			bead_ip.setLut(LUT.createLutFromColor(Color.WHITE));
			ImageRoi roi = new ImageRoi(0, 0, bead_ip);
			roi.setZeroTransparent(true);
			roi.setOpacity(0.3);
			Overlay ovl = new Overlay(roi);
			imp.setOverlay(ovl); // TODO save current overlay before overwriting
			messageArea.setText(count(bead_imp) + " bead(s) found");
		}
		if (dialog == 3) { // sprout preview
			if (bead_imp == null) // only if no preview was run on dialog 2
				bead_imp = findBeads(imp, ch_bead, use_bead_mask);
			ImagePlus sprout_imp = findSprouts(imp, ch_sprout, bead_imp);
			ImageProcessor sprout_ip = sprout_imp.getProcessor();
			sprout_ip.setLut(LUT.createLutFromColor(Color.WHITE));
			ImageRoi roi = new ImageRoi(0, 0, sprout_ip);
			roi.setZeroTransparent(true);
			roi.setOpacity(0.3);
			Overlay ovl = new Overlay(roi);
			imp.setOverlay(ovl);
			
		}
		if (dialog == 0) { // full processing
			processAndShow();
			// findBeads
			// findSprouts
			// skeletonize
			// analyze
			// output
			imp.setOverlay(null);
		}
		IJ.log("dialog = " + dialog);
	}

	/**
	 * Read parameters from ImageJ Prefs
	 *
	 * @param nChannels
	 */
	private void readPrefs (int nChannels) {
		/* Channels */
		ch_bead = 			(int)Prefs.get(PREF_KEY + "beads", 1);
		if (ch_bead > nChannels)
			ch_bead = 1;
		ch_sprout =			(int)Prefs.get(PREF_KEY + "sprouts", 2);
		if (ch_sprout > nChannels)
			ch_sprout = 1;
		ch_nuc =			(int)Prefs.get(PREF_KEY + "nucleus_marker", 1);
		if (ch_nuc > nChannels)
			ch_nuc = 1;
		use_bead_mask =			Prefs.get(PREF_KEY + "bead_mask", false);
		use_sprout_mask =		Prefs.get(PREF_KEY + "sprout_mask", false);
		use_nuc_mask =			Prefs.get(PREF_KEY + "nuc_mask", false);

		/* Bead recognition */
		thr_bead =			Prefs.get(PREF_KEY + "bead_threshold", "Triangle");
		blur_bead =			Prefs.get(PREF_KEY + "blur_radius_for_bead", 2.0);
		bead_radius =			Prefs.get(PREF_KEY + "minimum_bead_radius", 60.0);
		bead_radius_multiplier =	Prefs.get(PREF_KEY + "dilate_beads", 1.2);

		/* Sprout recognition */
		thr_sprout =			Prefs.get(PREF_KEY + "threshold", "Li");
		blur_sprout =			Prefs.get(PREF_KEY + "blur_radius_for_sprout", 4.0);

		/* Output */
		do_n_beads =			Prefs.get(PREF_KEY + "number_of_beads", true);
		do_n_sprouts =			Prefs.get(PREF_KEY + "number_of_sprouts", true);
		do_sprout_area =		Prefs.get(PREF_KEY + "total_sprout_area", true);
		do_total_length =		Prefs.get(PREF_KEY + "total_network_length", true);
		do_n_cells =			Prefs.get(PREF_KEY + "number_of_cells", true);
		do_avg_length =			Prefs.get(PREF_KEY + "average_sprout_length", true);
		do_avg_width =			Prefs.get(PREF_KEY + "average_sprout_width", true);
		do_cell_density =		Prefs.get(PREF_KEY + "cell_density", true);
		// do_pericyte =		Prefs.get(PREF_KEY + "pericyte_coverage", true);
	}

	/**
	 * Save parameters to ImageJ Prefs
	 */
	private void writePrefs () {
		/* Save settings to Prefs */
		Prefs.set(PREF_KEY + "beads", ch_bead);
		Prefs.set(PREF_KEY + "predefined_bead_mask", use_bead_mask);
		Prefs.set(PREF_KEY + "minimum_bead_radius", bead_radius);
		Prefs.set(PREF_KEY + "blur_radius_for_bead", blur_bead);
		Prefs.set(PREF_KEY + "sprouts", ch_sprout);
		Prefs.set(PREF_KEY + "predefined_sprout_mask", use_sprout_mask);
		Prefs.set(PREF_KEY + "recover_interrupted_structures", do_recover);
		Prefs.set(PREF_KEY + "bead_threshold", thr_bead);
		Prefs.set(PREF_KEY + "threshold", thr_sprout);
		Prefs.set(PREF_KEY + "blur_radius_for_sprout", blur_sprout);
		Prefs.set(PREF_KEY + "exclude_well-border_artefacts", do_exclude_borders);
		Prefs.set(PREF_KEY + "dilate_beads", bead_radius_multiplier);
		Prefs.set(PREF_KEY + "nucleus_marker", ch_nuc);
		//Prefs.set(PREF_KEY + "classify_endothelial", classify_ec);
		//Prefs.set(PREF_KEY + "endothelial_cell_nuclei", ch_endo);
		//Prefs.set(PREF_KEY + "quantify_pericyte_area", pericyte_area);
		//Prefs.set(PREF_KEY + "pericytes", ch_peri);
		Prefs.set(PREF_KEY + "number_of_beads", do_n_beads);
		Prefs.set(PREF_KEY + "number_of_sprouts", do_n_sprouts);
		Prefs.set(PREF_KEY + "total_sprout_area", do_sprout_area);
		Prefs.set(PREF_KEY + "total_network_length", do_total_length);
		Prefs.set(PREF_KEY + "number_of_cells", do_n_cells);
		Prefs.set(PREF_KEY + "average_sprout_length", do_avg_length);
		Prefs.set(PREF_KEY + "average_sprout_width", do_avg_width);
		Prefs.set(PREF_KEY + "cell_density", do_cell_density);
		//Prefs.set(PREF_KEY + "pericyte_coverage", do_pericyte);
	}

	/**
	 * Do the actual processing.
	 */
	private void processAndShow() {
		/* Private intermediate images */
		ImagePlus sprout_imp, ssp_imp, skel_imp, nuc_imp, endo_imp, peri_imp;
		ImageStack result_stack;
		/* Segmentation */
		bead_imp = findBeads(imp, ch_bead, use_bead_mask);

		if (use_sprout_mask) {
			sprout_imp = new Duplicator().run(imp, ch_sprout, ch_sprout, 1, 1, 1, 1); 
			IJ.run(sprout_imp, "Convert to Mask", ""); 
		} else {
			sprout_imp = findSprouts(imp, ch_sprout, bead_imp);
		}
		/* Morphometrical Analysis */
		// --- Number of beads ---      <= bead_imp
		num_beads = count(bead_imp);
		ImageCalculator ic = new ImageCalculator();
		ssp_imp = ic.run("OR create", sprout_imp, bead_imp);
		// --- Total sprout area ---      <= sprout_imp
		ResultsTable rt = new ResultsTable();
		IJ.setThreshold(sprout_imp, 1, 255);
		Analyzer an = new Analyzer(sprout_imp, Analyzer.AREA + Analyzer.LIMIT, rt);  
		an.measure();
		sprout_area = rt.getValueAsDouble(rt.getLastColumn(), rt.getCounter()-1);
		// --- Number of sprouts --- and --- Total length ---
		skel_imp = getCleanSkeleton(ssp_imp, bead_imp);
		if (analyzeSproutSkeleton(skel_imp, bead_imp)) {
			IJ.showStatus("Finished analyzing sprout skeletons");
		}
		// --- Number of cells --- and --- Pericyte coverage ---
		nuc_imp = getNucleusMask(imp, sprout_imp, ch_nuc);
		/*
		endo_imp = null;
		if (classify_ec) {
			endo_imp = classifyEC(imp, nuc_imp, ch_endo);
			num_nuc = count(endo_imp, 1);
			num_peri = count(endo_imp, 2);
			num_nuc += num_peri;
		} else {
			num_nuc = count(nuc_imp);
		}
		peri_imp = null;
		if (pericyte_area) {
			peri_imp = quantPericytes(imp, sprout_imp, ch_peri);
		}
		*/
		num_nuc = count(nuc_imp);

		/* Show the results and display result images */
		ResultsTable result = ResultsTable.getResultsTable();
		result.incrementCounter();
		result.setPrecision(5);
		result.addLabel(imp.getTitle());
		if (do_n_beads) result.addValue("n(beads)", num_beads);
		if (do_n_sprouts) result.addValue("n(sprouts)", num_sprouts);
		if (do_n_cells) result.addValue("n(cells)", num_nuc);
		if (do_sprout_area) result.addValue("Total sprout area (" + cal.getUnits() + "\u00B2)", sprout_area);
		if (do_total_length) result.addValue("Total network length (" + cal.getUnits() + ")", totalLength);
		if (do_avg_length) result.addValue("Average sprout length (" + cal.getUnits() + ")", avg_sprout_length);
		if (do_avg_width) result.addValue("Average sprout width (" + cal.getUnits() + ")", sprout_area / totalLength);
		if (do_cell_density) result.addValue("Cell density (1/" + cal.getUnits() + "\u00B2)", num_nuc / sprout_area);
		//if (do_pericyte & classify_ec) result.addValue("Pericytes per total cells", (double) num_peri / num_nuc);
		//if (do_pericyte & pericyte_area) result.addValue("Pericyte area fraction", peri_area / sprout_area);
		result.show("Results");

		/* Subtract beads from skeleton   */
		ic.run("Subtract", skel_imp, bead_imp);
		/* Show results stack */
		result_stack = bead_imp.getStack();
		result_stack.addSlice(ssp_imp.getProcessor());
		result_stack.addSlice(sprout_imp.getProcessor());
		result_stack.addSlice(skel_imp.getProcessor());
		result_stack.addSlice(nuc_imp.getProcessor()); // disabled for screencast
		/*
		if (classify_ec) {
			endo_imp.setSlice(1);
			result_stack.addSlice(endo_imp.getProcessor().duplicate());
			endo_imp.setSlice(2);
			result_stack.addSlice(endo_imp.getProcessor().duplicate());
		}
		if (pericyte_area) {
			result_stack.addSlice(peri_imp.getProcessor());
		}
		*/

		/*
		Overlay overlay = new Overlay();
		overlay.add(makeTransparentRoi(skel_imp, new Color(255, 255, 255, 128)));
		bead_imp.setOverlay(overlay);
		*/
		bead_imp.setTitle("ResultImage");
		// bead_imp.setSlice(2); // for Screencast
		bead_imp.show();

		/* Show result image with several overlays */
		// TODO: 
		/*
			float opacity = 30.0;

			roi = new ImageRoi(x, y, overlayimp.getProcessor());
			IndexColorModel cm = new IndexColorModel(8, 256, fi.reds, fi.greens, fi.blues);
		        ip.setColorModel(cm);
		        //roi.setName(overlay.getShortTitle());
		        ((ImageRoi)roi).setOpacity(opacity/100.0);
		            Overlay overlayList = imp.getOverlay();
		            if (overlayList==null) overlayList = new Overlay();
		            overlayList.add(roi);
		            imp.setOverlay(overlayList);
		            Undo.setup(Undo.OVERLAY_ADDITION, imp);
		 */

		/*
		bead_imp.setTitle("Bead mask");
		//bead_imp.show(); //debug
		sprout_imp.setTitle("Sprout mask");
		sprout_imp.show(); // debug
		ssp_imp.setTitle("Sprouts and bead mask");
		//ssp_imp.show(); // debug
		skel_imp.setTitle("Sprout skeleton");
		skel_imp.show(); // debug
		nuc_imp.setTitle("Nucleus mask");
		//nuc_imp.show(); //debug
		if (classify_ec) {
			endo_imp.setTitle("Classified nuclei");
			IJ.run(endo_imp, "Make Composite", "display=Composite"); //debug (this displays the image)
		}
		if (pericyte_area) {
			peri_imp.setTitle("Pericyte area");
			peri_imp.show(); //debug
		}
		*/		
	}

	/**
	 * Find beads in a given channel of an ImagePlus.
	 * 
	 * @param imp
	 * @param channel
	 */
	private ImagePlus findBeads(ImagePlus imp, int channel) {
		return findBeads(imp, channel, false);
	}

	/**
	 * Find beads in a given channel of an ImagePlus, or simply return the correct channel.
	 * 
	 * @param imp
	 * @param channel
	 * @param specified Use specified mask without processing
	 */
	private ImagePlus findBeads(ImagePlus imp, int channel, boolean specified) {
	 	ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
		if (!specified) {
			IJ.run(output, "Gaussian Blur...", "sigma=" + blur_bead + " scaled");
			IJ.setAutoThreshold(output, thr_bead + " dark");
			IJ.run(output, "Convert to Mask", "");
			ResultsTable rt = new ResultsTable(); // necessary to avoid interference with standard ResultsTable
			ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.IN_SITU_SHOW + ParticleAnalyzer.INCLUDE_HOLES, 0, rt, 10, 50000);
			pa.analyze(output);
			IJ.showStatus("Finding beads...");
			IJ.run(output, "Minimum...", "radius=" + IJ.d2s(bead_radius / pixel_size));
			IJ.run(output, "Maximum...", "radius=" + IJ.d2s(bead_radius_multiplier * bead_radius / pixel_size));
		}
		else {
			IJ.run(output, "Convert to Mask", "");
		}
		return output;
	 }

	/**
	 * Find sprouts in a given channel of an ImagePlus with a bead mask.
	 * 
	 * @param imp
	 * @param channel
	 * @param beads
	 */
	private ImagePlus findSprouts(ImagePlus imp, int channel, ImagePlus beads) {
		ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
		IJ.run(output, "Gaussian Blur...", "sigma=" + blur_sprout + " scaled");
		IJ.setAutoThreshold(output, thr_sprout + " dark"); // Use combined threshold here??
		ResultsTable rt = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 200, Double.POSITIVE_INFINITY);
		pa.analyze(output);
		IJ.showStatus("Finding sprouts...");
		if (do_recover) {
			/* Dilate and Erode with different radii */
			IJ.run(output, "Maximum...", "radius=" + IJ.d2s(10 / pixel_size));
			IJ.showStatus("Finding sprouts....");
			IJ.run(output, "Minimum...", "radius=" + IJ.d2s(8 / pixel_size));
			IJ.showStatus("Finding sprouts.....");
		}

		/* Find regions connected do beads */
		ImageCalculator ic = new ImageCalculator();
		ic.run("OR", output, beads);
		ResultsTable rt2 = new ResultsTable();
		ParticleAnalyzer pa2 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE | ParticleAnalyzer.RECORD_STARTS, ParticleAnalyzer.CENTER_OF_MASS, rt2, 0, Double.POSITIVE_INFINITY);
		pa2.analyze(beads);
		int nResults = rt2.getCounter();
		ImageProcessor ip = output.getProcessor();
		ip.setValue(128);
		FloodFiller filler = new FloodFiller(ip);
		for (int i = 0; i < nResults; i++) {
			filler.fill((int)rt2.getValue("XStart", i), (int)rt2.getValue("YStart", i));
		}
		IJ.setThreshold(output, 127, 129);
		IJ.run(output, "Convert to Mask", "");
		ic.run("XOR", output, beads);
		if (do_exclude_borders) {
			// TODO: discard big area at the corners
			ImagePlus artefact_mask = new Duplicator().run(output);
			// maybe watershed first? what about wholes?
			// analyze particles with exclude_edges
			ParticleAnalyzer pa4 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 0, Double.POSITIVE_INFINITY);
			pa4.analyze(artefact_mask);
			// xor the result with output
			ic.run("XOR", artefact_mask, output);
			// select particles with minimum size = min_artefact_size
			ParticleAnalyzer pa5 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 20000, Double.POSITIVE_INFINITY);
			pa5.analyze(artefact_mask);
			// XOR(output, artefact_mask)
			ic.run("XOR",output, artefact_mask);			
		}
		/* Discard sprouts smaller than min_sprout_area */
		ParticleAnalyzer pa3 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 5000, Double.POSITIVE_INFINITY);
		pa3.analyze(output);
	 	return output;
	}


	/**
	 * Count the number of objects in a segmented binary image.
	 * 
	 * @param imp Segmented binary image
	 */
	private int count(ImagePlus imp) {
		ResultsTable rt = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE, 0, rt, 0, Double.POSITIVE_INFINITY);
		pa.analyze(imp);
		return rt.getCounter();
	}

	/**
	 * Skeletonize a given sprout image and removes unimportant branches.
	 * 
	 * @param sprouts binary image containing sprout segmentation, the image to be skeletonized
	 * @param beads binary image containing bead segmentation
	 */
	private ImagePlus getCleanSkeleton(ImagePlus sprouts, ImagePlus beads) {
		ImagePlus output = new Duplicator().run(sprouts);
		IJ.run(output, "Skeletonize (2D/3D)", "");
		// TODO: remove "short" branches (pruning algorithm?)
		return output;
	}

	/**
	 * Populate result variables with values from skeleton analysis.
	 * 
	 * @param skeleton
	 * @param beads
	 * 
	 *  Determine the number of sprouts (num_sprouts) from the intersections
	 *  of the skeleton with the bead frames.
	 *  Determine the average sprout length (avg_sprout_length) by summing up
	 *  the lengths of all branches, and dividing by the determined number of
	 *  sprouts.
	 *  
	 *  Fill the following variables:
	 *   int	num_sprouts
	 *   double	avg_sprout_length
	 */
	private boolean analyzeSproutSkeleton(ImagePlus skeleton, ImagePlus beads) {
		ImagePlus temp = new Duplicator().run(beads);
		IJ.run(temp, "Dilate", "");
		ImageCalculator ic = new ImageCalculator();
		ic.run("XOR", temp, beads);
		ic.run("AND", temp, skeleton);
		/* Count the number of sprouts */
		ResultsTable rt = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE, 0, rt, 0, Double.POSITIVE_INFINITY);
		pa.analyze(temp);
		num_sprouts = rt.getCounter();

		ImagePlus sprout_skel = new Duplicator().run(skeleton);
		ic.run("Subtract", sprout_skel, beads);
		/* determine average network length per sprout */
		AnalyzeSkeleton_ skel = new AnalyzeSkeleton_();
		skel.calculateShortestPath = true;
		skel.setup("", sprout_skel);
		SkeletonResult sr = skel.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
		double[] branchLengths = sr.getAverageBranchLength();
		int[] branchNumbers = sr.getBranches();
		totalLength = 0;
		if (branchNumbers != null) {
			for (int i = 0; i < branchNumbers.length; i++) {
				totalLength += branchNumbers[i] * branchLengths[i];
			}
		}
		/* maybe count the longest_shortest_paths here, instead of total network length */
		avg_sprout_length = totalLength / num_sprouts;
		return true; // TODO: some error capturing here -> return false
	}

	/**
	 * Simple Segmentation of Nuclei in the given image.
	 * 
	 * @param imp
	 * @param sprouts
	 * @param channel
	 */
	private ImagePlus getNucleusMask(ImagePlus imp, ImagePlus sprouts, int channel) {
	 	ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
	 	IJ.run(output, "Subtract Background...", "rolling=50");
		IJ.run(output, "Gaussian Blur...", "sigma=2");
		IJ.setAutoThreshold(output, "Li dark");
		ResultsTable rt = new ResultsTable();		
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 20, Double.POSITIVE_INFINITY);
		pa.analyze(output);
		IJ.run(output, "Watershed", "");
		ImageCalculator ic = new ImageCalculator();
		ic.run("AND", output, sprouts);
		return output;
	}
}
