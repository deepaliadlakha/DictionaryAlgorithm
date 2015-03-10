function compdriver()
    [D_new1,coeff1]=nnsc_face(1,2,1,100,0.01,0.1);
    [D_new2,coeff2]=nnsc_face(1,2,2,100,0.01,0.1);
    
    
    height = 92;
    width = 112;
        
    Y_new2=D_new2*coeff2; 
    Y_new1=D_new1*coeff1; 
    
   
    for i=1:size(Y_new2,2),
       img2=reshape(Y_new2(:,i),width,height);
       img1=reshape(Y_new1(:,i),width,height);
       imagesc([ img1 img2 ]); axis equal tight off;colorbar;colormap gray;
       pause
    end

end