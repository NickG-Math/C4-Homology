function [howmanytimesthegenerator,varargout]=C4genmult(level,Data,varargin) 
%Given a symbolic product of asigmas ulambdas etc., get the degree of that product, and multiply the generator of that degree with the inverse of the negative exponents
%The result is howmanytimesthegenerator(1) times the generator.
%Next multiply the elements of the symbolic product that exist, the result
%is howmanytimesthegenerator(2) times the generator.


choice=[]; %The choice of generator. Only matters if there are multiple:
if size(varargin,2)==8 %This means there are two generators at some spot, so we need to choose
    choice=varargin{8};
    varargin(8)=[]; %Delete it from varargin so it doesn't bother us
end

if level==4
    basicgenerators{1}=[0,1,0]; %asigma
    basicgenerators{2}=[2,2,0]; %u2sigma
    basicgenerators{3}=[0,0,1]; %alambda
    basicgenerators{4}=[2,0,1]; %ulambda
    basicgenerators{5}=[-3,-3,0]; %w3
    basicgenerators{6}=[-3,-1,-1]; %x11
    basicgenerators{7}=[-3,0,-2]; %s3
elseif level==2
    basicgenerators{1}=[1,1,0]; %usigma
    basicgenerators{2}=[0,0,1]; %bar alambda
    basicgenerators{3}=[2,0,1]; %bar ulambda
    basicgenerators{4}=[-3,0,-2]; %bar s3
else
    basicgenerators{1}=[1,1,0]; %bar usigma
    basicgenerators{2}=[2,0,1]; %bar bar ulambda
end





generator=[0,0,0];
A=varargin;
element=cell(1,size(A,2));
for i=1:size(A,2)
        element{i}=A{i}*basicgenerators{i};
        generator=generator+element{i};
end
varargout{1}=generator;

if level==4
    limitofpositives=min(size(A,2),4);
elseif level==2
    limitofpositives=min(size(A,2),3);
else
    limitofpositives=size(A,2); %There are no negative elements anyway
end

%First we multiply the generator with the negative exponents so as to get
%undo the divisibilities. 
negativeindices=find([A{1:limitofpositives}]<0); %The negative indices
if isempty(negativeindices)
    howmanytimesthegeneratorneg=[];
else
    howmanytimesthegeneratorneg=1;  
    negativeproductwithgen=generator;
    for i=negativeindices
        firstneg=negativeproductwithgen;
        secondneg=-element{i};
        [~,prodneg]=C4mult(level,firstneg,secondneg,Data,choice); %We use the simplified basis. Don't try to compare resolutions with this answer!
        if isempty(prodneg)
            error('A homology group is 0')
        end
        howmanytimesthegeneratorneg=howmanytimesthegeneratorneg*prodneg;        
        negativeproductwithgen=-element{i}+negativeproductwithgen;
    end
end

%Multiply the elements with positive exponents for Euler+Orientation classes and the x_{n,m},s_3
positiveindices=find([A{1:limitofpositives}]>0); %The positive indices
positiveindices=[positiveindices,5:size(A,2)]; %The elements that exist.
if isempty(positiveindices)
    howmanytimesthegeneratorpos=[];
else
    howmanytimesthegeneratorpos=1;    
    positiveproduct=element{positiveindices(1)};
    positiveindices(1)=[];
    for i=positiveindices
        firstpos=positiveproduct;
        secondpos=element{i};
        [~,prodpos]=C4mult(level,firstpos,secondpos,Data);  %We use the simplified basis. Don't try to compare resolutions with this answer!
        if isempty(prodpos)
            error('A homology group is 0')
        end
        howmanytimesthegeneratorpos=howmanytimesthegeneratorpos*prodpos;
        positiveproduct=element{i}+positiveproduct;
    end
end

howmanytimesthegenerator=[howmanytimesthegeneratorneg,howmanytimesthegeneratorpos];

if isempty(howmanytimesthegenerator) %Only possible for a product with only 0 exponents i.e. 1
    howmanytimesthegenerator=1;
end