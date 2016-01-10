package sprouting;

/**
 * Sprout segmentation plugin for Fiji/ImageJ
 * 
 * @version 0.3.0
 * 
 * (C) 2012-2014 Jan Eglinger
 * Heinrich-Heine University Düsseldorf
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;

import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;

import ij.measure.Calibration;
import ij.measure.ResultsTable;

import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.plugin.Thresholder;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.PlugInFilterRunner;
import ij.plugin.filter.RankFilters;
import ij.plugin.frame.Recorder;

import ij.process.AutoThresholder;
import ij.process.FloodFiller;
import ij.process.ImageProcessor;
import ij.process.LUT;

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Font;
import java.awt.Label;

import java.util.Vector;

import morphology.BinaryReconstruct_;

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
	private static final boolean PERICYTES = false;
	private static final String PREF_KEY = "sprout_analyzer.";
	private static final double OVERLAY_OPACITY = 0.5;
	private static final int NO_DIALOG = 0, CHANNEL_DIALOG = 1, BEAD_DIALOG = 2, SPROUT_DIALOG = 3, NUCLEUS_DIALOG = 4, PERICYTE_DIALOG = 5, PERICYTE_AREA_DIALOG = 6;
	private static final int NUM_BEADS = 0, NUM_SPROUTS = 2, NUM_CELLS = 4, TOT_AREA = 6, TOT_LENGTH = 8, AVG_LENGTH = 1, AVG_WIDTH = 3, AVG_DENSITY = 5, NUM_EC = 7, PERI_AREA = 9; // custom order for param dialog
	private ImagePlus imp;
	private boolean is16Bit;
	private int flags = DOES_8G | DOES_16 | DOES_32 | NO_CHANGES;
	private int nPasses = 0;
	private int dialog = NO_DIALOG;
	private Label messageArea;
	private ImagePlus bead_imp, sprout_imp, nuc_imp;
	private GaussianBlur gb;
	private Thresholder thr;
	private RankFilters rf;
	private boolean userHasBlackBackground;

	/*  Image-dependent variables */
	private double pixel_size;
	private Calibration cal;

	/* 
	 * Parameters
	 */

	/* Channel configuration */
	private int ch_bead, ch_nuc, ch_sprout, ch_endo, ch_peri;
	private boolean use_bead_mask = false, use_sprout_mask = false, use_nuc_mask = false;

	/* Bead recognition */
	private String thr_bead;
	private double blur_bead, bead_radius, bead_radius_multiplier;

	/* Sprout recognition */
	private String thr_sprout;
	private double blur_sprout, min_plexus_area, min_sprout_area, max_hole_area, min_cluster_size;
	private boolean do_recover = false;
	private boolean do_exclude_borders;

	/* Nuclei */
	private String thr_nuc;
	private double min_nuc_area, blur_nuc, max_tolerance;

	/* EC classification */
	private String thr_ec;
	private double min_ec_area;
	private boolean quant_cell_numbers, quant_cell_fraction;

	/* Pericyte area */
	private String thr_peri;
	private boolean quant_peri_area, quant_peri_fraction;

	/* Output */
	private boolean quantify[];

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
		if (imp == null) return flags;
		is16Bit = (imp.getType() == ImagePlus.GRAY16);
		cal = imp.getCalibration();
		pixel_size = cal.getX(1.0);
		thr = new Thresholder();
		gb = new GaussianBlur();
		rf = new RankFilters();
		userHasBlackBackground = Prefs.blackBackground; // get user-set value
		Prefs.blackBackground = true;// set blackbackground to true for this plugin 
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

		// TODO warn if imp is uncalibrated

		// read prefs
		readPrefs(nChannels);
		
		// Prepare CheckboxGroup
		String[] labels = {
			"Number_of_beads",			"Average_sprout_length",
			"Number_of_sprouts",		"Average_sprout_width",
			"Number_of_cells",			"Cell_density",
			"Total_sprout_area",		"Numbers of ECs/pericytes",
			"Total_network_length",		"Pericyte_coverage (%area)"
		};

		// Dialog #1
		dialog = CHANNEL_DIALOG;
		GenericDialog gd1 = new GenericDialog("Sprout Analyzer - Configuration");
		gd1.setInsets(5, 0, 0);
		gd1.addMessage("Parameters", bold);
		gd1.addCheckboxGroup(5, 2, labels, quantify);
		
		gd1.setInsets(10, 0, 0);
		gd1.addMessage("Channels", bold);
		gd1.addChoice("Beads", channels, channels[ch_bead - 1]);
		gd1.addChoice("Sprouts", channels, channels[ch_sprout - 1]);
		gd1.addChoice("Nuclei", channels, channels[ch_nuc - 1]);
		gd1.addChoice("Endothelial_cell_marker", channels, channels[ch_endo - 1]);
		gd1.addChoice("Pericytes", channels, channels[ch_peri - 1]);

		gd1.addDialogListener(this);

		dialogItemChanged(gd1, null);
		gd1.showDialog();

		if (gd1.wasCanceled())
			return DONE;

		if (!use_bead_mask) {
			// Dialog #2
			dialog = BEAD_DIALOG;
			GenericDialog gd2 = new GenericDialog("Sprout Analyzer - Bead detection");
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

			gd2.showDialog();

			if (gd2.wasCanceled())
				return DONE;
		}

		if (!use_sprout_mask) {
			// Dialog #3
			dialog = SPROUT_DIALOG;
			GenericDialog gd3 = new GenericDialog("Sprout Analyzer - Sprout detection");
			//gd3.setInsets(0, 0, 0);
			//gd3.addMessage("Sprout detection____________________________________________________________", bold);
			gd3.addChoice("Sprout_threshold method", AutoThresholder.getMethods(), thr_sprout);
			gd3.addSlider("Blur_radius_for_sprout detection (" + cal.getUnits() + "):", 0.05, 5.0, blur_sprout);
			gd3.addSlider("Minimal_plexus_area (" + cal.getUnits() + "\u00B2):", 1000, 50000, min_plexus_area);
			gd3.addSlider("Minimal_sprout_area (" + cal.getUnits() + "\u00B2):", 10, 5000, min_sprout_area);
			gd3.addCheckbox("Exclude_cell_clusters at borders", do_exclude_borders);
			gd3.addSlider("Cluster_size for exclusion (" + cal.getUnits() + "\u00B2):", 1000, 100000, min_cluster_size);
			// TODO add do_recover and do_exclude_borders options
			gd3.addPreviewCheckbox(pfr, "Preview sprout detection");
			gd3.addDialogListener(this);

			gd3.showDialog();

			if (gd3.wasCanceled())
				return DONE;
		}

		if (!use_nuc_mask && (quantify[NUM_CELLS] || quantify[AVG_DENSITY] || quantify[NUM_EC])) {
			// Dialog #4
			dialog = NUCLEUS_DIALOG;
			GenericDialog gd4 = new GenericDialog("Sprout Analyzer - Nucleus segmentation");
			gd4.addChoice("Nucleus_threshold method", AutoThresholder.getMethods(), thr_nuc);
			gd4.addSlider("Blur_radius_for_nucleus segmentation (" + cal.getUnits() + "):", 0.05, 5.0, blur_nuc);
			gd4.addSlider("Tolerance for nuclei separation:", 0, is16Bit ? 2000 : 20, max_tolerance);
			gd4.addSlider("Minimal_nucleus_area (" + cal.getUnits() + "\u00B2):", 0, 200, min_nuc_area);
			gd4.addPreviewCheckbox(pfr, "Preview nucleus detection");
			gd4.addDialogListener(this);

			gd4.showDialog();

			if (gd4.wasCanceled())
				return DONE;
		}

		if (quantify[NUM_EC]) {
			// Dialog #5
			dialog = PERICYTE_DIALOG;
			GenericDialog gd5 = new GenericDialog("Sprout Analyzer - Cell classification");
			gd5.addChoice("EC_threshold method", AutoThresholder.getMethods(), thr_ec);
			gd5.addSlider("Minimal_EC_area (" + cal.getUnits() + "\u00B2):", 0, 200, min_ec_area);
			gd5.addCheckbox("Measure_number_of_classified cells", quant_cell_numbers);
			gd5.addCheckbox("Measure_pericyte_cell_fraction", quant_cell_fraction);
			gd5.addPreviewCheckbox(pfr, "Preview cell classification"); // TODO: separate preview from parameter choice
			gd5.addDialogListener(this);

			gd5.showDialog();

			if (gd5.wasCanceled())
				return DONE;
		}

		if (quantify[PERI_AREA]) {
			// Dialog #6
			dialog = PERICYTE_AREA_DIALOG;
			GenericDialog gd6 = new GenericDialog("Sprout Analyzer - Pericyte coverage");
			gd6.addChoice("Pericyte_threshold method", AutoThresholder.getMethods(), thr_peri);
			gd6.addCheckbox("Measure_total_pericyte_area", quant_peri_area);
			gd6.addCheckbox("Measure_pericyte_area_fraction", quant_peri_fraction);
			gd6.addPreviewCheckbox(pfr, "Preview pericyte area");
			gd6.addDialogListener(this);

			gd6.showDialog();

			if (gd6.wasCanceled())
				return DONE;
		}

		// write prefs
		writePrefs();

		dialog = NO_DIALOG;
		return flags;
	}

	/**
	 * dialogItemChanged
	 */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		switch (dialog) {
			case CHANNEL_DIALOG:	// read parameters from dialog 1
				//IJ.log("ItemChanged 1");
				//use_bead_mask = gd.getNextBoolean();
				//use_sprout_mask = gd.getNextBoolean();
				quantify[NUM_BEADS] = gd.getNextBoolean();
				quantify[AVG_LENGTH] = gd.getNextBoolean();
				quantify[NUM_SPROUTS] = gd.getNextBoolean();
				quantify[AVG_WIDTH] = gd.getNextBoolean();
				quantify[NUM_CELLS] = gd.getNextBoolean();
				quantify[AVG_DENSITY] = gd.getNextBoolean();
				quantify[TOT_AREA] = gd.getNextBoolean();
				quantify[NUM_EC] = gd.getNextBoolean();
				quantify[TOT_LENGTH] = gd.getNextBoolean();
				quantify[PERI_AREA] = gd.getNextBoolean();
				//use_nuc_mask = gd.getNextBoolean();
				ch_bead = gd.getNextChoiceIndex() + 1;
				ch_sprout = gd.getNextChoiceIndex() + 1;
				ch_nuc = gd.getNextChoiceIndex() + 1;
				ch_endo = gd.getNextChoiceIndex() + 1;
				ch_peri = gd.getNextChoiceIndex() + 1;

				Vector<?> choices = gd.getChoices();
				Choice choice2 = (Choice) choices.get(2);
				Choice choice3 = (Choice) choices.get(3);
				Choice choice4 = (Choice) choices.get(4);
				choice2.setEnabled(quantify[NUM_CELLS] || quantify[AVG_DENSITY] || quantify[NUM_EC] ? true : false);
				choice3.setEnabled(quantify[NUM_EC] ? true : false);
				choice4.setEnabled(quantify[PERI_AREA] ? true : false);
				break;
			case BEAD_DIALOG: // read parameters from dialog 2
				//IJ.log("ItemChanged 2");
				thr_bead = gd.getNextChoice();
				blur_bead = gd.getNextNumber();
				bead_radius = gd.getNextNumber();
				bead_radius_multiplier = gd.getNextNumber();
				if (!gd.getPreviewCheckbox().getState()) {
					messageArea.setText(""); // clear "x beads found" message
					imp.setOverlay(null);
				}
				break;
			case SPROUT_DIALOG: // read parameters from dialog 3
				//IJ.log("ItemChanged 3");
				thr_sprout = gd.getNextChoice();
				blur_sprout = gd.getNextNumber();
				min_plexus_area = gd.getNextNumber();
				min_sprout_area = gd.getNextNumber();
				do_exclude_borders = gd.getNextBoolean();
				min_cluster_size = gd.getNextNumber();
				if (!gd.getPreviewCheckbox().getState())
					imp.setOverlay(null);
				break;
			case NUCLEUS_DIALOG: // read parameters from dialog 4
				//IJ.log("ItemChanged 4");
				thr_nuc = gd.getNextChoice();
				blur_nuc = gd.getNextNumber();
				max_tolerance = gd.getNextNumber();
				min_nuc_area = gd.getNextNumber();
				if (!gd.getPreviewCheckbox().getState())
					imp.setOverlay(null);
				break;
			case PERICYTE_DIALOG: // read parameters from dialog 5
				//IJ.log("ItemChanged 5");
				thr_ec = gd.getNextChoice();
				min_ec_area = gd.getNextNumber();
				quant_cell_numbers = gd.getNextBoolean();
				quant_cell_fraction = gd.getNextBoolean();
				if (!gd.getPreviewCheckbox().getState())
					imp.setOverlay(null);
				break;
			case PERICYTE_AREA_DIALOG: // read parameters from dialog 6
				//IJ.log("ItemChanged 6");
				thr_peri = gd.getNextChoice();
				quant_peri_area = gd.getNextBoolean();
				quant_peri_fraction = gd.getNextBoolean();				
				if (!gd.getPreviewCheckbox().getState())
					imp.setOverlay(null);
				break;
		}
		return (!gd.invalidNumber());
	}

	/**
	 * run
	 */
	@Override
	public void run (ImageProcessor ip) {
		//Recorder.recordInMacros = false;
		// TODO save current overlay before overwriting
		imp.setOverlay(null);
		if (dialog == BEAD_DIALOG) { // bead preview
			bead_imp = findBeads(imp, ch_bead); // find beads (takes time)
			imp.setOverlay(makeOverlay(bead_imp.getProcessor(), Color.WHITE, OVERLAY_OPACITY));
			messageArea.setText(count(bead_imp) + " bead(s) found");
		}
		if (dialog == SPROUT_DIALOG) { // sprout preview
			if (bead_imp == null) // only if no preview was run on dialog 2
				bead_imp = findBeads(imp, ch_bead, use_bead_mask);
			sprout_imp = findSprouts(imp, ch_sprout, bead_imp, use_sprout_mask);
			imp.setOverlay(makeOverlay(sprout_imp.getProcessor(), Color.WHITE, OVERLAY_OPACITY));
			
		}
		if (dialog == NUCLEUS_DIALOG) { // nuclei preview
			if (sprout_imp == null) {
				if (bead_imp == null)
					bead_imp = findBeads(imp, ch_bead, use_bead_mask);
				sprout_imp = findSprouts(imp, ch_sprout, bead_imp, use_sprout_mask);
			}
			nuc_imp = getNucleusMask(imp, sprout_imp, ch_nuc);
			if (null != nuc_imp) imp.setOverlay(makeOverlay(nuc_imp.getProcessor(), Color.WHITE, OVERLAY_OPACITY));
		}
		if (dialog == PERICYTE_DIALOG) { // cell classification preview
			if (nuc_imp == null) {
				if (sprout_imp == null) {
					if (bead_imp == null)
						bead_imp = findBeads(imp, ch_bead, use_bead_mask);
					sprout_imp = findSprouts(imp, ch_sprout, bead_imp, use_sprout_mask);
				}
				nuc_imp = getNucleusMask(imp, sprout_imp, ch_nuc);
			}
			ImagePlus endo_imp = classifyEC(imp, nuc_imp, ch_endo);
			if (null != endo_imp) {
				imp.setOverlay(makeDoubleOverlay(endo_imp.getStack().getProcessor(1), endo_imp.getStack().getProcessor(2), Color.YELLOW, Color.MAGENTA, OVERLAY_OPACITY));
			}
		}
		if (dialog == PERICYTE_AREA_DIALOG) { // pericyte area preview
			if (sprout_imp == null) {
				if (bead_imp == null)
					bead_imp = findBeads(imp, ch_bead, use_bead_mask);
				sprout_imp = findSprouts(imp, ch_sprout, bead_imp, use_sprout_mask);
			}
			ImagePlus peri_imp = getPericyteArea(imp, sprout_imp, ch_peri);
			if (null != peri_imp) imp.setOverlay(makeOverlay(peri_imp.getProcessor(), Color.WHITE, OVERLAY_OPACITY));			
		}
		if (dialog == NO_DIALOG) { // full processing
			processAndShow();
			Prefs.blackBackground = userHasBlackBackground;
		}
	}

	/**
	 * Read parameters from ImageJ Prefs
	 *
	 * @param nChannels
	 */
	private void readPrefs (int nChannels) {
		/* Channels */
		ch_bead = 			(int)Prefs.get(PREF_KEY + "beads", 1);
		if (ch_bead > nChannels || ch_bead < 1)
			ch_bead = 1;
		ch_sprout =			(int)Prefs.get(PREF_KEY + "sprouts", 2);
		if (ch_sprout > nChannels || ch_sprout < 1)
			ch_sprout = 1;
		ch_nuc =			(int)Prefs.get(PREF_KEY + "nucleus_marker", 1);
		if (ch_nuc > nChannels || ch_nuc < 1)
			ch_nuc = 1;
		ch_endo =			(int)Prefs.get(PREF_KEY + "endothelial_cell_nuclei", 3);
		if (ch_endo > nChannels || ch_endo < 1)
			ch_endo = 1;
		ch_peri =			(int)Prefs.get(PREF_KEY + "pericyte_marker", 4);
		if (ch_peri > nChannels || ch_peri < 1)
			ch_peri = 1;
		use_bead_mask =			Prefs.get(PREF_KEY + "bead_mask", false);
		use_sprout_mask =		Prefs.get(PREF_KEY + "sprout_mask", false);
		use_nuc_mask =			Prefs.get(PREF_KEY + "nuc_mask", false);

		/* Bead recognition */
		thr_bead =				Prefs.get(PREF_KEY + "bead_threshold", "Triangle");
		blur_bead =				Prefs.get(PREF_KEY + "blur_radius_for_bead", 2.0);
		bead_radius =			Prefs.get(PREF_KEY + "minimum_bead_radius", 60.0);
		bead_radius_multiplier =	Prefs.get(PREF_KEY + "dilate_beads", 1.2);

		/* Sprout recognition */
		thr_sprout =			Prefs.get(PREF_KEY + "sprout_threshold", "Li");
		blur_sprout =			Prefs.get(PREF_KEY + "blur_radius_for_sprout", 4.0);
		min_plexus_area =		Prefs.get(PREF_KEY + "minimum_plexus_area", 5000);
		min_sprout_area =		Prefs.get(PREF_KEY + "minimum_sprout_area", 100);
		do_exclude_borders =	Prefs.get(PREF_KEY + "exclude_cell_clusters", false);
		min_cluster_size =		Prefs.get(PREF_KEY + "minimum_cluster_for_exclusion", 20000);

		/* Nucleus segmentation */
		thr_nuc =				Prefs.get(PREF_KEY + "nucleus_threshold", "Minimum");
		blur_nuc =				Prefs.get(PREF_KEY + "blur_radius_for_nuclei", 1.0);
		max_tolerance =			Prefs.get(PREF_KEY + "nucleus_tolerance", 5);
		min_nuc_area =			Prefs.get(PREF_KEY + "minimum_nucleus_area", 10);

		/* EC classification */
		thr_ec =				Prefs.get(PREF_KEY + "ec_threshold", "Default");
		min_ec_area =			Prefs.get(PREF_KEY + "minimum_ec_area", 10);
		quant_cell_numbers =	Prefs.get(PREF_KEY + "quantify_cell_numbers", true);
		quant_cell_fraction =	Prefs.get(PREF_KEY + "quantify_cell_fraction", true);

		/* Pericyte area */
		thr_peri =				Prefs.get(PREF_KEY + "pericyte_threshold", "Default");
		quant_peri_area =		Prefs.get(PREF_KEY + "quantify_pericyte_area", false);
		quant_peri_fraction =	Prefs.get(PREF_KEY + "quantify_pericyte_fraction", true);

		/* Output */
		quantify = new boolean[10];
		quantify[NUM_BEADS] =	Prefs.get(PREF_KEY + "number_of_beads", true);
		quantify[NUM_SPROUTS] =	Prefs.get(PREF_KEY + "number_of_sprouts", true);
		quantify[TOT_AREA] =	Prefs.get(PREF_KEY + "total_sprout_area", true);
		quantify[TOT_LENGTH] =	Prefs.get(PREF_KEY + "total_network_length", true);
		quantify[NUM_CELLS] =	Prefs.get(PREF_KEY + "number_of_cells", true);
		quantify[AVG_LENGTH] =	Prefs.get(PREF_KEY + "average_sprout_length", true);
		quantify[AVG_WIDTH] =	Prefs.get(PREF_KEY + "average_sprout_width", true);
		quantify[AVG_DENSITY] =	Prefs.get(PREF_KEY + "cell_density", true);
		quantify[NUM_EC] =		Prefs.get(PREF_KEY + "ec_number", false);
		quantify[PERI_AREA] =	Prefs.get(PREF_KEY + "pericyte_coverage", false);
		/* for (int j = 0; j <= quantify.length; j++) {
			IJ.log("Pos: " + Integer.toString(j) + " " + Boolean.toString(quantify[j]));
		} */
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
		Prefs.set(PREF_KEY + "sprout_threshold", thr_sprout);
		Prefs.set(PREF_KEY + "nucleus_threshold", thr_nuc);
		Prefs.set(PREF_KEY + "blur_radius_for_sprout", blur_sprout);
		Prefs.set(PREF_KEY + "minimum_plexus_area", min_plexus_area);
		Prefs.set(PREF_KEY + "minimum_sprout_area", min_sprout_area);
		Prefs.set(PREF_KEY + "minimum_nucleus_area", min_nuc_area);
		Prefs.set(PREF_KEY + "exclude_cell_clusters", do_exclude_borders);
		Prefs.set(PREF_KEY + "minimum_cluster_for_exclusion", min_cluster_size);		
		Prefs.set(PREF_KEY + "dilate_beads", bead_radius_multiplier);
		Prefs.set(PREF_KEY + "nucleus_marker", ch_nuc);
		Prefs.set(PREF_KEY + "blur_radius_for_nuclei", blur_nuc);
		Prefs.set(PREF_KEY + "nucleus_tolerance", max_tolerance);
		Prefs.set(PREF_KEY + "endothelial_cell_nuclei", ch_endo);
		Prefs.set(PREF_KEY + "pericyte_marker", ch_peri);
		Prefs.set(PREF_KEY + "number_of_beads", quantify[NUM_BEADS]);
		Prefs.set(PREF_KEY + "number_of_sprouts", quantify[NUM_SPROUTS]);
		Prefs.set(PREF_KEY + "total_sprout_area", quantify[TOT_AREA]);
		Prefs.set(PREF_KEY + "total_network_length", quantify[TOT_LENGTH]);
		Prefs.set(PREF_KEY + "number_of_cells", quantify[NUM_CELLS]);
		Prefs.set(PREF_KEY + "average_sprout_length", quantify[AVG_LENGTH]);
		Prefs.set(PREF_KEY + "average_sprout_width", quantify[AVG_WIDTH]);
		Prefs.set(PREF_KEY + "cell_density", quantify[AVG_DENSITY]);
		Prefs.set(PREF_KEY + "ec_number", quantify[NUM_EC]);
		Prefs.set(PREF_KEY + "pericyte_coverage", quantify[PERI_AREA]);
		Prefs.set(PREF_KEY + "ec_threshold", thr_ec);
		Prefs.set(PREF_KEY + "minimum_ec_area", min_ec_area);
		Prefs.set(PREF_KEY + "pericyte_threshold", thr_peri);
		Prefs.set(PREF_KEY + "quantify_pericyte_area", quant_peri_area);
		Prefs.set(PREF_KEY + "quantify_pericyte_fraction", quant_peri_fraction);
		Prefs.set(PREF_KEY + "quantify_cell_numbers", quant_cell_numbers);
		Prefs.set(PREF_KEY + "quantify_cell_fraction", quant_cell_fraction);
	}

	/**
	 * Do the actual processing.
	 */
	private void processAndShow() {
		/* Private intermediate images */
		ImagePlus ssp_imp, skel_imp, endo_imp = null, peri_imp = null;
		ImageStack result_stack;
		/* Segmentation */
		bead_imp = findBeads(imp, ch_bead, use_bead_mask);
		sprout_imp = findSprouts(imp, ch_sprout, bead_imp, use_sprout_mask);
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
		if (quantify[NUM_CELLS] || quantify[AVG_DENSITY] || quantify[NUM_EC]) {
			nuc_imp = getNucleusMask(imp, sprout_imp, ch_nuc);
			if (quantify[NUM_EC]) {
				endo_imp = classifyEC(imp, nuc_imp, ch_endo);
				num_nuc = count(endo_imp, 1);
				num_peri = count(endo_imp, 2);
				num_nuc += num_peri;
			} else {
				num_nuc = count(nuc_imp);
			}
		}
		
		if (quantify[PERI_AREA]) {
			peri_imp = getPericyteArea(imp, sprout_imp, ch_peri);
		}

		/* Show the results and display result images */
		ResultsTable result = ResultsTable.getResultsTable();
		result.incrementCounter();
		result.setPrecision(5);
		result.addLabel(imp.getTitle());
		if (quantify[NUM_BEADS]) result.addValue("n(beads)", num_beads);
		if (quantify[NUM_SPROUTS]) result.addValue("n(sprouts)", num_sprouts);
		if (quantify[NUM_CELLS]) result.addValue("n(cells)", num_nuc);
		if (quantify[TOT_AREA]) result.addValue("Total sprout area (" + cal.getUnits() + "\u00B2)", sprout_area);
		if (quantify[TOT_LENGTH]) result.addValue("Total network length (" + cal.getUnits() + ")", totalLength);
		if (quantify[AVG_LENGTH]) result.addValue("Average sprout length (" + cal.getUnits() + ")", avg_sprout_length);
		if (quantify[AVG_WIDTH]) result.addValue("Average sprout width (" + cal.getUnits() + ")", sprout_area / totalLength);
		if (quantify[AVG_DENSITY]) result.addValue("Cell density (1/" + cal.getUnits() + "\u00B2)", num_nuc / sprout_area);
		if (quantify[NUM_EC]) {
			if (quant_cell_numbers) {
				result.addValue("Number of ECs", (int) num_nuc - num_peri);
				result.addValue("Number of Pericytes", (int) num_peri);
			}
			if (quant_cell_fraction) result.addValue("Pericytes per total cells", (double) num_peri / num_nuc);
		}
		// TODO: optionally include total numbers of EC and pericytes
		if (quantify[PERI_AREA]) {
			peri_area = measureArea(peri_imp);
			if (quant_peri_area) result.addValue("Total pericyte area (" + cal.getUnits() + "\u00B2)", peri_area);
			if (quant_peri_fraction) result.addValue("Pericyte area fraction", peri_area / sprout_area);
		}
		result.show("Results");

		/* Subtract beads from skeleton   */
		ic.run("Subtract", skel_imp, bead_imp);
		/* Show results stack */
		result_stack = bead_imp.getStack();
		result_stack.addSlice(ssp_imp.getProcessor());
		result_stack.addSlice(sprout_imp.getProcessor());
		result_stack.addSlice(skel_imp.getProcessor());
		if (quantify[NUM_CELLS] || quantify[AVG_DENSITY] || quantify[NUM_EC])
			result_stack.addSlice(nuc_imp.getProcessor()); // disabled for screencast
		if (quantify[NUM_EC]) {
			result_stack.addSlice(endo_imp.getStack().getProcessor(1));
			result_stack.addSlice(endo_imp.getStack().getProcessor(2));			
		}
		if (quantify[PERI_AREA]) {
			result_stack.addSlice(peri_imp.getProcessor());
		}

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
			//IJ.run(output, "Gaussian Blur...", "sigma=" + blur_bead + " scaled");
			gb.blurGaussian(output.getProcessor() ,blur_bead ,blur_bead , 0.02); // scaled!
			//IJ.setAutoThreshold(output, thr_bead + " dark");
			output.getProcessor().setAutoThreshold(thr_bead, true, 0);
			//IJ.run(output, "Convert to Mask", "");
			WindowManager.setTempCurrentImage(output);
			thr.run("mask");
			WindowManager.setTempCurrentImage(null);
			// ResultsTable rt = new ResultsTable(); // necessary to avoid interference with standard ResultsTable
			ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.IN_SITU_SHOW + ParticleAnalyzer.INCLUDE_HOLES, 0, null, 10, 50000);
			pa.analyze(output);
			IJ.showStatus("Finding beads...");
			//IJ.run(output, "Minimum...", "radius=" + IJ.d2s(bead_radius / pixel_size));
			rf.rank(output.getProcessor(), bead_radius / pixel_size, RankFilters.MIN);
			//IJ.run(output, "Maximum...", "radius=" + IJ.d2s(bead_radius_multiplier * bead_radius / pixel_size));
			rf.rank(output.getProcessor(), bead_radius_multiplier * bead_radius / pixel_size, RankFilters.MAX);
		} else {
			IJ.run(output, "Convert to Mask", ""); // (new Thresholder()).run("mask");
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
	private ImagePlus findSprouts(ImagePlus imp, int channel, ImagePlus beads, boolean specified) {
		ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
		if (!specified) {
			//IJ.run(output, "Gaussian Blur...", "sigma=" + blur_sprout + " scaled");
			gb.blurGaussian(output.getProcessor() ,blur_sprout ,blur_sprout , 0.02);
			//IJ.setAutoThreshold(output, thr_sprout + " dark"); // Use combined threshold here??
			output.getProcessor().setAutoThreshold(thr_sprout, true, 0);
			// ResultsTable rt = new ResultsTable();
			ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, null, min_plexus_area, Double.POSITIVE_INFINITY);
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
			IJ.setThreshold(output, 127, 129); // TODO: avoid that this gets recorded
			//IJ.run(output, "Convert to Mask", "");
			WindowManager.setTempCurrentImage(output);
			thr.run("mask");
			WindowManager.setTempCurrentImage(null);
			ic.run("XOR", output, beads);
			if (do_exclude_borders) {
				// TODO: discard big area at the corners
				ImagePlus artefact_mask = new Duplicator().run(output);
				// maybe watershed first? what about wholes?
				// analyze particles with exclude_edges
				ParticleAnalyzer pa4 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES + ParticleAnalyzer.IN_SITU_SHOW, 0, null, 0, Double.POSITIVE_INFINITY);
				pa4.analyze(artefact_mask);
				// xor the result with output
				ic.run("XOR", artefact_mask, output);
				// select particles with minimum size = min_artefact_size
				ParticleAnalyzer pa5 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, null, min_cluster_size, Double.POSITIVE_INFINITY);
				pa5.analyze(artefact_mask);
				// XOR(output, artefact_mask)
				ic.run("XOR",output, artefact_mask);			
			}
			/* Discard sprouts smaller than min_sprout_area */
			ParticleAnalyzer pa3 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, null, min_sprout_area, Double.POSITIVE_INFINITY);
			pa3.analyze(output);
		} else {
			IJ.run(output, "Convert to Mask", ""); // (new Thresholder()).run("mask");
		}
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
	 * Counts the number of objects in a given slice of binary image stack
	 * 
	 * @param imp
	 * @param slice
	 */
	public int count(ImagePlus imp, int slice) {
		imp.setSlice(slice);
		return count(imp);
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
		IJ.run(output, "Subtract Background...", "rolling=50"); // TODO: avoid use of IJ.run during preview
		IJ.run(output, "Gaussian Blur...", "sigma=" + blur_nuc + " scaled"); // scaling?
		IJ.setAutoThreshold(output, thr_nuc + " dark");
		ImageProcessor ip = output.getProcessor();
		ip = (new MaximumFinder()).findMaxima(ip, max_tolerance, ip.getMinThreshold(), MaximumFinder.SEGMENTED, false, false);
		if (null != ip) output.setProcessor(ip);// make sure this finds its way back to output
		/*
		IJ.run(output, "Subtract Background...", "rolling=50");
		IJ.run(output, "Gaussian Blur...", "sigma=2");
		IJ.setAutoThreshold(output, "Li dark"); // customize threshold
		*/
		// ResultsTable rt = new ResultsTable();		
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, null, min_nuc_area / (pixel_size * pixel_size), Double.POSITIVE_INFINITY);
		pa.analyze(output);
		//IJ.run(output, "Watershed", "");
		ImageCalculator ic = new ImageCalculator();
		ic.run("AND", output, sprouts);
		return output;
	}

	/**
	 * Classify endothelial cells based on EC-specific nuclear staining
	 * 
	 * @param imp
	 * @param nuclei
	 * @param channel
	 */
	private ImagePlus classifyEC(ImagePlus imp, ImagePlus nuclei, int channel) {
	 	/* Create EC-positive mask */
	 	ImagePlus temp = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
	 	IJ.run(temp, "Subtract Background...", "rolling=50");
		IJ.run(temp, "Gaussian Blur...", "sigma=2"); // TODO: make blur radius configurable
		ImageCalculator ic = new ImageCalculator();
		ImagePlus output = ic.run("Multiply create 32-bit", temp, nuclei);
		IJ.setAutoThreshold(output, thr_ec + " dark"); // TODO: avoid IJ during preview
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, null, min_ec_area / (pixel_size * pixel_size), Double.POSITIVE_INFINITY);
		pa.analyze(output);
		// use BinaryReconstruct output > nuclei
		BinaryReconstruct_ br = new BinaryReconstruct_();
		ImagePlus ec_imp = (ImagePlus)br.exec(nuclei, output, null, false, true, false)[1];
		// XOR ec X nuclei -> non-EC
		ImagePlus ec_neg = ic.run("XOR create", ec_imp, nuclei);
		ec_imp.getStack().addSlice(ec_neg.getProcessor());
		return ec_imp;
	}

	/**
	 * Quantifies the pericyte coverage (area fraction)
	 * 
	 * @param imp
	 * @param sprouts
	 * @param channel
	 */
	public ImagePlus getPericyteArea(ImagePlus imp, ImagePlus sprouts, int channel) {
		ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
		IJ.setAutoThreshold(output, thr_peri + " dark"); // TODO: avoid IJ
		IJ.run(output, "Convert to Mask", ""); // TODO: avoid IJ
		// mask with sprouts
		ImageCalculator ic = new ImageCalculator();
		ic.run("AND", output, sprouts);
		/* // this is just quantification, not needed during preview
		ResultsTable rt = new ResultsTable();
		IJ.setThreshold(output, 1, 255);
		Analyzer an = new Analyzer(output, Analyzer.AREA + Analyzer.LIMIT, rt);  
		an.measure();
		peri_area = rt.getValueAsDouble(rt.getLastColumn(), rt.getCounter()-1);
		*/
		return output;
	}

	/**
	 * Measure area of a binary image
	 * 
	 * @param binaryImp The ImagePlus binary image to be measured
	 */
	private Double measureArea(ImagePlus binaryImp) {
		ResultsTable rt = new ResultsTable();
		IJ.setThreshold(binaryImp, 1, 255);
		Analyzer an = new Analyzer(binaryImp, Analyzer.AREA + Analyzer.LIMIT, rt);  
		an.measure();
		return rt.getValueAsDouble(rt.getLastColumn(), rt.getCounter()-1);
		
	}

	/**
	 * Transform an ImageProcessor into an Overlay with a given color and opacity
	 *
	 * @param ip The ImageProcessor that should serve as overlay
	 * @param color The color to create the Overlay's LUT
	 * @param opacity The opacity of the overlay
	 */
	private Overlay makeOverlay(ImageProcessor ip, Color color, Double opacity) {
		ip.setLut(LUT.createLutFromColor(color));
		ImageRoi roi = new ImageRoi(0, 0, ip);
		roi.setZeroTransparent(false);
		roi.setOpacity(opacity);
		return new Overlay(roi);
	}

	/**
	 * Transform an ImageProcessor into an Overlay with a given color and opacity
	 *
	 * @param ip The ImageProcessor that should serve as overlay
	 * @param color The color to create the Overlay's LUT
	 * @param opacity The opacity of the overlay
	 */
	private Overlay makeDoubleOverlay(ImageProcessor ip1, ImageProcessor ip2, Color color1, Color color2, Double opacity) {
		ip1.setLut(LUT.createLutFromColor(color1));
		ip2.setLut(LUT.createLutFromColor(color2));
		ImageRoi roi1 = new ImageRoi(0, 0, ip1);
		ImageRoi roi2 = new ImageRoi(0, 0, ip2);
		roi1.setZeroTransparent(true);
		roi2.setZeroTransparent(false);
		roi1.setOpacity(opacity);
		roi2.setOpacity(opacity);
		Overlay ovl = new Overlay();
		ovl.add(roi2);
		ovl.add(roi1);
		return ovl;
	}
}
