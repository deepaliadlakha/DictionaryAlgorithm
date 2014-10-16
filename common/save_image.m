function save_image(img, outp_file_name, colmap)
myNumOfColors = 200;
myColorScale = [(0:1/(myNumOfColors-1):1)' , (0:1/(myNumOfColors-1):1)' , (0:1/(myNumOfColors-1):1)'];
imshow(img); 
colormap (myColorScale);
if colmap==1
    colormap jet;
else
    colormap gray;
end 
daspect ([1 1 1]); 
axis tight;
colorbar;
export_fig(outp_file_name, '-png')
