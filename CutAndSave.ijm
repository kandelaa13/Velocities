setBatchMode(true);

for (i = 0; i < 3; i++) {

run("Image Sequence...", "open=D:/Bendings/11112022/20221111_noY27/20221111_BF_Stack_10x/20221111_BF_Stack_10x_f0000_t0000.ome.tif file=f"+IJ.pad(i, 4)+" sort");

// Convert into a hyperstack of 15 z-slices, 26 timepoints
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=15 frames=26 display=Color");

// Keep only slice number 10, since it's the one in good focus
run("Make Substack...", "slices=10 frames=1-26");

saveAs("Tiff", "D:/Bendings/11112022/20221111_noY27/20221111_BF_10x/20221111_BF_10x_f"+ IJ.pad(i, 4) +".tif");

run("Close All");

}

for (i = 3; i < 8; i++) {

run("Image Sequence...", "open=D:/Bendings/11112022/20221111_noY27/20221111_BF_Stack_10x/20221111_BF_Stack_10x_f0000_t0000.ome.tif file=f"+IJ.pad(i, 4)+" sort");

// Convert into a hyperstack of 15 z-slices, 26 timepoints
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=15 frames=25 display=Color");

// Keep only slice number 10, since it's the one in good focus
run("Make Substack...", "slices=10 frames=1-25");

saveAs("Tiff", "D:/Bendings/11112022/20221111_noY27/20221111_BF_10x/20221111_BF_10x_f"+ IJ.pad(i, 4) +".tif");

run("Close All");

}