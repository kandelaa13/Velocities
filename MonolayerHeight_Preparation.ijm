dir = "C:/Users/kballerini/Desktop/Lab/11112022/20221111_noY27/20221111_Fluo_Infrared_Stack_10x/";
file = "20221111_Fluo_Infrared_Stack_10x_f";
position = 7;

slices = 40; frames = 25;
pixelsize = 0.65; pixeldepth = 2;

run("Image Sequence...", "open="+dir+file+IJ.pad(position,4)+"_t0000.ome.tif file=f"+IJ.pad(position,4)+"_ sort");

run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+d2s(slices,0)+" frames="+d2s(frames,0)+" display=Color");

Dialog.create("How many ROIs?");
Dialog.addNumber("Number", 2);

num_rectangles = Dialog.getNumber();

run("Properties...", "channels=1 slices="+d2s(slices,0)+" frames="+d2s(frames,0)+" unit=µm pixel_width="+d2s(pixelsize,2)+" pixel_height="+d2s(pixelsize,2)+" voxel_depth="+d2s(pixeldepth,0));

for (i = 0; i < num_rectangles; i++) {

	makeRectangle(396, 468, 1275, 168);
	waitForUser("Place ROI, then press OK");
	
	run("Reslice [/]...", "output=2.000 start=Top flip");
	run("Z Project...", "projection=[Average Intensity] all");
	run("Smooth", "stack");
	run("Threshold...");
	waitForUser("Ajust the threshold, normally 106, 107");
	run("Convert to Mask", "method=Otsu background=Dark calculate");
	run("Analyze Particles...", "size=300-Infinity add stack");
	roiManager("Deselect");
	
	roiManager("Save", "C:/Users/kballerini/Desktop/Lab/11112022/20221111_noY27/RoiSet_f"+IJ.pad(position,4)+"_"+d2s(i,0)+".zip");
	selectWindow("20221111_Fluo_Infrared_Stack_10x");
}
run("Close All");
