function element=Homologyelement(element,SmithVariables)

%Q0=SmithVariables{1};
Q0i=SmithVariables{1};

zerovectors=SmithVariables{2};
P1=SmithVariables{3};
modoutcompletely=SmithVariables{4};
Homology=SmithVariables{5};

for i=1:size(element,2)
    if isempty(element) || isequal(element{i},0)  
        continue
    elseif isempty(P1)
        element{i}=0; %If some P1 is provided empty that means Generator=Homology=0
        continue
    end   
    
 %   element{i}=round(Q0\element{i});%Our element in terms of x,y,z,.. and the zero vectors  %Rounding is unfortunately necessary in some ranges
    element{i}=Q0i*element{i};%Our element in terms of x,y,z,.. and the zero vectors  %Rounding is unfortunately necessary in some ranges

    element{i}(zerovectors)=[]; %Just like the kernel
    element{i}=P1*element{i}; %Expressing the element in terms of u,v,w,...
    element{i}(modoutcompletely)=[]; %Mod out the image.
    if isempty(element{i})
        element{i}=0;
    end
    for j=1:size(Homology,2)
        if Homology(j)~=1
            element{i}(j)=mod(element{i}(j), Homology(j));
        end
    end
    element{i}=element{i}'; %Have the element as an array, reflecting the generator
end
end