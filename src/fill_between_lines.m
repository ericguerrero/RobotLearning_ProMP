function fill_between_lines(X,Y1,Y2,C,alpha )
% This function fills the area between two lines Y1 & Y2 
%  where X is the x values of the lines. The area will have
% the color specified by C (RGB coordinates) and 
% a transparency specified by alpha.

% check if the formast is correct: patch function needs vector as 1xn
 if (size(X,1)==1 || size(X,2)==1 )
     if(size(X,2)==1)
         X=X';
     end
 else
     error('fill_between_lines is supposed to work only with one dimensional vectors. X is not one dimensional.');
 end
 
 if(size(Y1,1)==1 || size(Y1,2)==1 )
     if(size(Y1,2)==1)
         Y1=Y1';
     end
 else
     error('fill_between_lines is supposed to work only with one dimensional vectors. Y1 is not one dimensional.');
 end
 
 if(size(Y2,1)==1 || size(Y2,2)==1 )
     if(size(Y2,2)==1)
         Y2=Y2';
     end
 else
     error('fill_between_lines is supposed to work only with one dimensional vectors. Y2 is not one dimensional.');
 end
 
 if(length(X) == length(Y1) & length(X)==length(Y2))
     patch( [X fliplr(X)],  [Y1 fliplr(Y2)], C ,'FaceAlpha',alpha);
 else
     error('The input vectors have not the same length');
 end
end