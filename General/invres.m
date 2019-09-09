function top=invres(bottom,rankbottom,level)
%
%INPUT: Column "bottom", array "rankbottom" and int "level"
%
%OUTPUT: Column "top"
%
%DESCRIPTION: The inverse of restriction (for free Mackey functors). 
%
%Given "bottom" at "level" return "top" one level higher with Res(top)=bottom (assuming that "bottom" is in the image of the Res)
%
%Sloppy code but it's used so so rarely that it really has no performance impact.

top=[];
trackhor=0;
transferlevel=2*level;
for i=1:size(rankbottom,2)
    next=[];
    if rankbottom(i)>4/transferlevel
        next=bottom(trackhor+1:trackhor+rankbottom(i)/2);
    else
        next=bottom(trackhor+1:trackhor+rankbottom(i));
    end
    top=[top;next];
    trackhor=trackhor+rankbottom(i);
end
end