package sprouting;

/**
 * Sprout segmentation plugin for Fiji and ImageJ
 * 
 * @version 1.7
 * 
 * (C) 2012-2013 Jan Eglinger
 * Heinrich-Heine University Düsseldorf
 */

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.TextField;

import java.awt.image.BufferedImage;

import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;

import ij.measure.Calibration;
import ij.measure.ResultsTable;

import ij.plugin.PlugIn;
import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ParticleAnalyzer;

import ij.process.ImageProcessor;
import ij.process.AutoThresholder;
import ij.process.FloodFiller;

import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.ImageRoi;

import skeleton_analysis.AnalyzeSkeleton_;
import skeleton_analysis.SkeletonResult;

/**
 * Sprout segmentation plugin
 * 
 * @author Jan Eglinger (jan.eglinger at gmail.com)
 */
public class Sprout_Analyzer implements PlugIn, DialogListener {

	/*
	 * Private variables
	 */

	/*  Configuration options   */
	private double bead_radius;
	private double blur_bead, blur_sprout;
	private double min_plexus_area, min_sprout_area, min_bead_radius, min_nuc_area;
	private double bead_radius_multiplier, max_hole_area;
	private String thr_sprout, thr_bead;

	/*  Channel numbers   */
	private int ch_bead, ch_sprout, ch_nuc, ch_endo, ch_peri;

	/* Analysis options */
	private boolean use_bead_mask, use_sprout_mask, do_recover, do_exclude_borders, classify_ec, pericyte_area;
	private boolean do_n_beads, do_n_sprouts, do_sprout_area, do_total_length, do_n_cells;
	private boolean do_avg_length, do_avg_width, do_cell_density, do_pericyte;

	/*  Image-dependent variables */
	private double pixel_size;
	private Calibration cal;

	/*  Results   */
	private int num_beads, num_sprouts, num_nuc, num_peri;
	private double sprout_area, avg_sprout_length, avg_sprout_width, cell_density, totalLength, peri_area;

	/**
	 * Run method.
	 * 
	 * This method creates a GenericDialog and runs the plugin.
	 */
	public void run(String arg) {
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		cal = imp.getCalibration();
		pixel_size = cal.getX(1.0);

		GenericDialog gd = initializeDialog(imp.getNChannels());
		gd.showDialog();

		if (gd.wasOKed()) {
			/* Channel choices  */
			ch_bead = gd.getNextChoiceIndex() + 1;
			thr_bead = gd.getNextChoice(); // Threshold method for use in bead detection
			ch_sprout = gd.getNextChoiceIndex() + 1;
			thr_sprout = gd.getNextChoice(); // Threshold method for use in sprout/cell detection
			ch_nuc = gd.getNextChoiceIndex() + 1;
			ch_endo = gd.getNextChoiceIndex() + 1;
			ch_peri = gd.getNextChoiceIndex() + 1;
			/* Analysis parameters  */
			bead_radius = gd.getNextNumber();
			blur_bead = gd.getNextNumber();
			bead_radius_multiplier = gd.getNextNumber();
			blur_sprout = gd.getNextNumber();
				//min_plexus_area
				//min_sprout_area
				//min_nuc_area
				//max_hole_area
			/* Predefined image masks */
			use_bead_mask = gd.getNextBoolean();
			use_sprout_mask = gd.getNextBoolean();
			do_recover = gd.getNextBoolean();
			do_exclude_borders = gd.getNextBoolean();
			classify_ec = gd.getNextBoolean();
			pericyte_area = gd.getNextBoolean();
			do_n_beads = gd.getNextBoolean();
			do_avg_length = gd.getNextBoolean();
			do_n_sprouts = gd.getNextBoolean();
			do_avg_width = gd.getNextBoolean();
			do_sprout_area = gd.getNextBoolean();
			do_cell_density = gd.getNextBoolean();
			do_total_length = gd.getNextBoolean();
			do_pericyte = gd.getNextBoolean();
			do_n_cells = gd.getNextBoolean();

			processAndShow(imp); // this starts the actual processing

			/* Save settings to Prefs */
			Prefs.set("sprout_analyzer.beads", ch_bead);
			Prefs.set("sprout_analyzer.predefined_bead_mask", use_bead_mask);
			Prefs.set("sprout_analyzer.minimum_bead_radius", bead_radius);
			Prefs.set("sprout_analyzer.blur_radius_for_bead", blur_bead);
			Prefs.set("sprout_analyzer.sprouts", ch_sprout);
			Prefs.set("sprout_analyzer.predefined_sprout_mask", use_sprout_mask);
			Prefs.set("sprout_analyzer.recover_interrupted_structures", do_recover);
			Prefs.set("sprout_analyzer.bead_threshold", thr_bead);
			Prefs.set("sprout_analyzer.threshold", thr_sprout);
			Prefs.set("sprout_analyzer.blur_radius_for_sprout", blur_sprout);
			Prefs.set("sprout_analyzer.exclude_well-border_artefacts", do_exclude_borders);
			Prefs.set("sprout_analyzer.dilate_beads", bead_radius_multiplier);
			Prefs.set("sprout_analyzer.nucleus_marker", ch_nuc);
			Prefs.set("sprout_analyzer.classify_endothelial", classify_ec);
			Prefs.set("sprout_analyzer.endothelial_cell_nuclei", ch_endo);
			Prefs.set("sprout_analyzer.quantify_pericyte_area", pericyte_area);
			Prefs.set("sprout_analyzer.pericytes", ch_peri);
			Prefs.set("sprout_analyzer.number_of_beads", do_n_beads);
			Prefs.set("sprout_analyzer.number_of_sprouts", do_n_sprouts);
			Prefs.set("sprout_analyzer.total_sprout_area", do_sprout_area);
			Prefs.set("sprout_analyzer.total_network_length", do_total_length);
			Prefs.set("sprout_analyzer.number_of_cells", do_n_cells);
			Prefs.set("sprout_analyzer.average_sprout_length", do_avg_length);
			Prefs.set("sprout_analyzer.average_sprout_width", do_avg_width);
			Prefs.set("sprout_analyzer.cell_density", do_cell_density);
			Prefs.set("sprout_analyzer.pericyte_coverage", do_pericyte);
		}
	}

	/**
	 * Does the actual processing.
	 * 
	 * @param imp
	 */
	private void processAndShow(ImagePlus imp) {
		/* Private intermediate images */
		ImagePlus bead_imp, sprout_imp, ssp_imp, skel_imp, nuc_imp, endo_imp, peri_imp;
		ImageStack result_stack;
		/* Segmentation */
		if (use_bead_mask) {
			bead_imp = new Duplicator().run(imp, ch_bead, ch_bead, 1, 1, 1, 1);
			IJ.run(bead_imp, "Convert to Mask", ""); 
		} else {
			bead_imp = findBeads(imp, ch_bead);
		}

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
		if (do_pericyte & classify_ec) result.addValue("Pericytes per total cells", (double) num_peri / num_nuc);
		if (do_pericyte & pericyte_area) result.addValue("Pericyte area fraction", peri_area / sprout_area);
		result.show("Results");

		/* Subtract beads from skeleton   */
		ic.run("Subtract", skel_imp, bead_imp);
		/* Show results stack */
		result_stack = bead_imp.getStack();
		result_stack.addSlice(ssp_imp.getProcessor());
		result_stack.addSlice(sprout_imp.getProcessor());
		result_stack.addSlice(skel_imp.getProcessor());
		// result_stack.addSlice(nuc_imp.getProcessor()); // disabled for screencast
		if (classify_ec) {
			endo_imp.setSlice(1);
			result_stack.addSlice(endo_imp.getProcessor().duplicate());
			endo_imp.setSlice(2);
			result_stack.addSlice(endo_imp.getProcessor().duplicate());
		}
		if (pericyte_area) {
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
	 * Finds beads in a given channel of an ImagePlus.
	 * 
	 * @param imp
	 * @param channel
	 */
	public ImagePlus findBeads(ImagePlus imp, int channel) {
	 	ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
		IJ.run(output, "Gaussian Blur...", "sigma=" + blur_bead + " scaled");
		IJ.setAutoThreshold(output, thr_bead + " dark");
		IJ.run(output, "Convert to Mask", "");
		ResultsTable rt = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.IN_SITU_SHOW + ParticleAnalyzer.INCLUDE_HOLES, 0, rt, 10, 50000);
		pa.analyze(output);
		IJ.showStatus("Finding beads...");
		IJ.run(output, "Minimum...", "radius=" + IJ.d2s(bead_radius / pixel_size));
		IJ.run(output, "Maximum...", "radius=" + IJ.d2s(bead_radius_multiplier * bead_radius / pixel_size));
		return output;
	 }

	/**
	 * Finds sprouts in a given channel of an ImagePlus with a bead mask.
	 * 
	 * @param imp
	 * @param channel
	 * @param beads
	 */
	public ImagePlus findSprouts(ImagePlus imp, int channel, ImagePlus beads) {
		ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
		IJ.run(output, "Gaussian Blur...", "sigma=" + blur_sprout + " scaled");
		IJ.setAutoThreshold(output, thr_sprout + " dark"); // Use combined threshold here??
		// Additional step of analyze particles, fclose(10)
		ResultsTable rt = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 200, Double.POSITIVE_INFINITY);
		pa.analyze(output);
		IJ.showStatus("Finding sprouts...");

		/* Single close operation
			float morphRad = Math.round(2 * 10 / pixel_size) / 2; // Morphology accepts values .0 and .5
			StructureElement se1 = new StructureElement(StructureElement.CIRCLE, 1, morphRad, StructureElement.OFFSET0);
			MorphoProcessor mp1 = new MorphoProcessor(se1);
			mp1.fclose(output.getProcessor());
		*/

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
	 * Counts the number of objects in a segmented binary image.
	 * 
	 * @param imp Segmented binary image
	 */
	public int count(ImagePlus imp) {
		ResultsTable rt = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE, 0, rt, 0, Double.POSITIVE_INFINITY);
		pa.analyze(imp);
		return rt.getCounter();
	}

	/**
	 * Skeletonizes a given sprout image and removes unimportant branches.
	 * 
	 * @param sprouts binary image containing sprout segmentation, the image to be skeletonized
	 * @param beads binary image containing bead segmentation
	 */
	public ImagePlus getCleanSkeleton(ImagePlus sprouts, ImagePlus beads) {
		ImagePlus output = new Duplicator().run(sprouts);
		IJ.run(output, "Skeletonize (2D/3D)", "");
		// TODO: remove "short" branches (pruning algorithm?)
		return output;
	}

	/**
	 * Populates result variables with values from skeleton analysis.
	 * 
	 * @param skeleton
	 * @param beads
	 * 
	 *  Determines the number of sprouts (num_sprouts) from the intersections
	 *  of the skeleton with the bead frames.
	 *  Determines the average sprout length (avg_sprout_length) by summing up
	 *  the lengths of all branches, and dividing by the determined number of
	 *  sprouts.
	 *  
	 *  Fill the following variables:
	 *   int	num_sprouts
	 *   double	avg_sprout_length
	 */
	public boolean analyzeSproutSkeleton(ImagePlus skeleton, ImagePlus beads) {
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
	public ImagePlus getNucleusMask(ImagePlus imp, ImagePlus sprouts, int channel) {
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

	/**
	 * Classify endothelial cells based on EC-specific nuclear staining
	 * 
	 * @param imp
	 * @param nuclei
	 * @param channel
	 */
	public ImagePlus classifyEC(ImagePlus imp, ImagePlus nuclei, int channel) {
	 	/* Create EC-positive mask */
	 	ImagePlus temp = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
	 	IJ.run(temp, "Subtract Background...", "rolling=50");
		IJ.run(temp, "Gaussian Blur...", "sigma=2");
		ImageCalculator ic = new ImageCalculator();
		ImagePlus output = ic.run("Multiply create 32-bit", temp, nuclei);
		IJ.setAutoThreshold(output, "Default dark");
		ResultsTable rt = new ResultsTable();		
		ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 20, Double.POSITIVE_INFINITY);
		pa.analyze(output);
		IJ.run(output, "Watershed", "");
		/* Create EC-negative mask */
		ImagePlus big_ec = new Duplicator().run(output);
		IJ.run(big_ec, "Dilate", "");
		IJ.run(big_ec, "Dilate", "");
		ImagePlus ec_neg = new Duplicator().run(nuclei);
		ic.run("Subtract", ec_neg, big_ec);
		// Watershed ??
		ParticleAnalyzer pa2 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 0, rt, 30, Double.POSITIVE_INFINITY, 0.6, 1.0);
		pa2.analyze(ec_neg);
		output.getStack().addSlice(ec_neg.getProcessor());
		return output;
	}

	/**
	 * Quantifies the pericyte coverage (area fraction)
	 * 
	 * @param imp
	 * @param sprouts
	 * @param channel
	 */
	public ImagePlus quantPericytes(ImagePlus imp, ImagePlus sprouts, int channel) {
		ImagePlus output = new Duplicator().run(imp, channel, channel, 1, 1, 1, 1);
		IJ.setAutoThreshold(output, "Default dark"); // TODO: factor out the threshold method
		IJ.run(output, "Convert to Mask", "");
		// mask with sprouts
		ImageCalculator ic = new ImageCalculator();
		ic.run("AND", output, sprouts);
		ResultsTable rt = new ResultsTable();
		IJ.setThreshold(output, 1, 255);
		Analyzer an = new Analyzer(output, Analyzer.AREA + Analyzer.LIMIT, rt);  
		an.measure();
		peri_area = rt.getValueAsDouble(rt.getLastColumn(), rt.getCounter()-1);
		return output;
	}

	/**
	 * Creates a transparent Roi from a binary image for use as overlay
	 * 
	 * @param imp
	 * @param fg - Foreground color with alpha coding
	 */
	public Roi makeTransparentRoi(ImagePlus imp, Color fg) {
		ImageProcessor ip = imp.getProcessor();
		int width = ip.getWidth();
		int height = ip.getHeight();
		//ip = ip.convertToByte(true);
		BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics g = bi.getGraphics();
		for (int x=0; x<width; x++) {
			for (int y=0; y<height; y++) {
				int v = ip.get(x, y);
				if (v>1) {
					g.setColor(fg); // color here!
					g.fillRect(x, y, 1, 1);
				}
			}
		}
        	Roi roi = new ImageRoi(0, 0, bi);
		return roi;
	}

	/**
	 * Initializes and builds a GenericDialog with values from ij.Prefs
	 * 
	 * @param nChannels - the number of channels in the image to be processed
	 */
	public GenericDialog initializeDialog(int nChannels) {
		GenericDialog gd = new GenericDialog("Sprout Analyzer - multiparametric sprout morphometry");
		String[] channels = new String[nChannels];
		for (int i=0; i<nChannels; i++) {
			channels[i] = "Channel " + (i+1);
		}
		/* Initialize settings from Prefs or default values */
		ch_bead = 		(int)Prefs.get("sprout_analyzer.beads", 1);
		if (ch_bead > nChannels) ch_bead = 1;
		use_bead_mask =		Prefs.get("sprout_analyzer.predefined_bead_mask", false);
		bead_radius =		Prefs.get("sprout_analyzer.minimum_bead_radius", 60);
		blur_bead =		Prefs.get("sprout_analyzer.blur_radius_for_bead", 2);
		bead_radius_multiplier = Prefs.get("sprout_analyzer.dilate_beads", 1);
		ch_sprout =		(int)Prefs.get("sprout_analyzer.sprouts", 2);
		if (ch_sprout > nChannels) ch_sprout = 1;
		use_sprout_mask =	Prefs.get("sprout_analyzer.predefined_sprout_mask", false);
		do_recover =		Prefs.get("sprout_analyzer.recover_interrupted_structures", true);
		thr_bead =		Prefs.get("sprout_analyzer.bead_threshold", "Triangle");
		thr_sprout =		Prefs.get("sprout_analyzer.threshold", "Li");
		blur_sprout =		Prefs.get("sprout_analyzer.blur_radius_for_sprout", 4);
		do_exclude_borders =	Prefs.get("sprout_analyzer.exclude_well-border_artefacts", true);
		ch_nuc =		(int)Prefs.get("sprout_analyzer.nucleus_marker", 1);
		if (ch_nuc > nChannels) ch_nuc = 1;
		classify_ec =		Prefs.get("sprout_analyzer.classify_endothelial", true);
		ch_endo =		(int)Prefs.get("sprout_analyzer.endothelial_cell_nuclei", 1);
		if (ch_endo > nChannels) ch_endo = 1;
		pericyte_area =		Prefs.get("sprout_analyzer.quantify_pericyte_area", true);
		ch_peri =		(int)Prefs.get("sprout_analyzer.pericytes", 1);
		if (ch_peri > nChannels) ch_peri = 1;

		do_n_beads =		Prefs.get("sprout_analyzer.number_of_beads", true);
		do_n_sprouts =		Prefs.get("sprout_analyzer.number_of_sprouts", true);
		do_sprout_area =	Prefs.get("sprout_analyzer.total_sprout_area", true);
		do_total_length =	Prefs.get("sprout_analyzer.total_network_length", true);
		do_n_cells =		Prefs.get("sprout_analyzer.number_of_cells", true);
		do_avg_length =		Prefs.get("sprout_analyzer.average_sprout_length", true);
		do_avg_width =		Prefs.get("sprout_analyzer.average_sprout_width", true);
		do_cell_density =	Prefs.get("sprout_analyzer.cell_density", true);
		do_pericyte =		Prefs.get("sprout_analyzer.pericyte_coverage", true);
		
		/* Header */
		//gd.setInsets(0, 0, 0);
		//gd.addMessage("Sprout analyzer - multiparametric sprout morphometry", new Font("my", Font.BOLD, 16));
		/* Bead detection */
		gd.setInsets(0, 0, 0);
		gd.addMessage("Bead detection______________________________________________________________", new Font("my", Font.BOLD, 12));
		gd.setInsets(5, 20, 0);
		gd.addChoice("Beads are in", channels, channels[ch_bead - 1]); /* Choice choice0 */
		gd.setInsets(10, 170, 10);
		gd.addCheckbox("Predefined_bead_mask", use_bead_mask); /* Checkbox box0 */
		gd.setInsets(0, 20, 0);
		gd.addChoice("Bead_finding_threshold method", AutoThresholder.getMethods(), thr_bead); /* Choice choice1 */
		gd.setInsets(0, 20, 0);
		gd.addNumericField("Minimum_bead_radius:", bead_radius, 1, 5, cal.getUnits() ); /* TextField num0 */
		gd.setInsets(0, 20, 0);
		gd.addNumericField("Blur_radius_for_bead detection:", blur_bead, 1, 5, cal.getUnits() ); /* TextField num1 */
		gd.setInsets(0, 20, 0);
		gd.addNumericField("Dilate_beads by factor:", bead_radius_multiplier, 1, 5, "" ); /* TextField num2 */
		/* Sprout detection */
		gd.setInsets(15, 0, 0);
		gd.addMessage("Sprout detection____________________________________________________________", new Font("my", Font.BOLD, 12));
		gd.setInsets(5, 20, 0);
		gd.addChoice("Sprouts are in", channels, channels[ch_sprout - 1]); /* Choice choice2 */
		gd.setInsets(10, 170, 10);
		gd.addCheckbox("Predefined_sprout_mask", use_sprout_mask); /* Checkbox box1 */
		gd.setInsets(0, 250, 0);
		gd.addCheckbox("Recover_interrupted_structures", do_recover); /* Checkbox box2 */
		gd.setInsets(0, 20, 0);
		gd.addChoice("Threshold method", AutoThresholder.getMethods(), thr_sprout); /* Choice choice3 */
		gd.setInsets(0, 20, 0);
		gd.addNumericField( "Blur_radius_for_sprout detection:", blur_sprout, 1, 5, cal.getUnits() ); /* TextField num3 */
		gd.setInsets(0, 250, 0);
		gd.addCheckbox("Exclude_well-border_artefacts", do_exclude_borders); /* Checkbox box3 */
		/* Cell analysis */
		gd.setInsets(15, 0, 0);
		gd.addMessage("Cell analysis_______________________________________________________________", new Font("my", Font.BOLD, 12));
		gd.setInsets(5, 20, 0);
		gd.addChoice("Nucleus_marker", channels, channels[ch_nuc -1]); /* Choice choice4 */
		gd.setInsets(10, 170, 10);
		gd.addCheckbox("Classify_endothelial cell nuclei", classify_ec); /* Checkbox box4 */
		gd.setInsets(0, 20, 0);
		gd.addChoice("Endothelial_cell_nuclei", channels, channels[ch_endo - 1]); /* Choice choice5 */
		gd.setInsets(10, 170, 10);
		gd.addCheckbox("Quantify_pericyte_area", pericyte_area); /* Checkbox box5 */
		gd.setInsets(0, 20, 0);
		gd.addChoice("Pericytes", channels, channels[ch_peri - 1]); /* Choice choice6 */
		/* Output */
		gd.setInsets(15, 0, 0);
		gd.addMessage("Output_____________________________________________________________________", new Font("my", Font.BOLD, 12));
		String[] labels = {
			"Number_of_beads",		"Average_sprout_length",
			"Number_of_sprouts",		"Average_sprout_width",
			"Total_sprout_area",		"Cell_density",
			"Total_network_length",		"Pericyte_coverage",
			"Number_of_cells",		"",
		};
		boolean[] defaults = {
			do_n_beads,			do_avg_length,
			do_n_sprouts,			do_avg_width,
			do_sprout_area,			do_cell_density,
			do_total_length,		do_pericyte,
			do_n_cells,			false
		};
		gd.setInsets(10, 50, 0);
		gd.addCheckboxGroup(5, 2, labels, defaults);
		gd.setInsets(20, 50, 0);
		gd.addCheckbox("Display_result_image_stack", true);
		// TODO: add option to display result stack, result image with overlay, or nothing
		gd.addDialogListener(this);
		dialogItemChanged(gd, null); // disable the appropriate fields
		return gd;		
	}

	/**
	 * DialogListener method to disable irrelevant fields.
	 * 
	 * @param gd
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		Vector<?> boxes = gd.getCheckboxes();
		Checkbox box0 = (Checkbox) boxes.get(0);
		Checkbox box1 = (Checkbox) boxes.get(1);
		Checkbox box2 = (Checkbox) boxes.get(2);
		Checkbox box3 = (Checkbox) boxes.get(3);
		Checkbox box4 = (Checkbox) boxes.get(4);
		Checkbox box5 = (Checkbox) boxes.get(5);
		Checkbox box13 = (Checkbox) boxes.get(13);
		Vector<?> numbers = gd.getNumericFields();
		TextField num0 = (TextField) numbers.get(0);
		TextField num1 = (TextField) numbers.get(1);
		TextField num2 = (TextField) numbers.get(2);
		TextField num3 = (TextField) numbers.get(3);
		Vector<?> choices = gd.getChoices();
		Choice choice1 = (Choice) choices.get(1);
		Choice choice3 = (Choice) choices.get(3);
		Choice choice5 = (Choice) choices.get(5);
		Choice choice6 = (Choice) choices.get(6);
		if(box0.getState()) {
			num0.setEnabled(false);
			choice1.setEnabled(false);
			num1.setEnabled(false);
			num2.setEnabled(false);
		} else {
			num0.setEnabled(true);
			choice1.setEnabled(true);
			num1.setEnabled(true);
			num2.setEnabled(true);
		}
		if(box1.getState()) {
			box2.setEnabled(false);
			choice3.setEnabled(false);
			num3.setEnabled(false);
			box3.setEnabled(false);			
		} else {
			box2.setEnabled(true);
			choice3.setEnabled(true);
			num3.setEnabled(true);
			box3.setEnabled(true);
		}
		if(box4.getState()) {
			choice5.setEnabled(true);
		} else {
			choice5.setEnabled(false);
		}
		if(box5.getState()) {
			choice6.setEnabled(true);
		} else {
			choice6.setEnabled(false);
		}
		if(box4.getState() | box5.getState()) {
			box13.setEnabled(true);
		} else {
			box13.setEnabled(false);
		}

		return true;
	}

}
