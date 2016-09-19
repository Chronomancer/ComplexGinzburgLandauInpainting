% Simple implementation of digital image inpainting with Complex Ginzburg-Landau equation in Matlab
% Based on the paper, Digital Inpainting Using the Complex Ginzburg-Landau Equation by Grossauer Harald

function ImgOut = CGL(Img0, mask0, TIME, dt, eps)
% Inputs:
% Img0 - grey scale input image with damages
% mask0 - binary image that marks the areas of the damages
% TIME - maximum number of iterations
% dt - time step for iterative solution
% eps - epislon, less than 1, I used 0.2 ~ 0.5

% Output:
% ImgOut - Output image in grey scale

lambda=0.028; 
maskwidth = 2;

root2 = sqrt(2);

[N M]=size(Img0);

ImgOut=double(Img0);
ImgLaplacian=zeros(N,M);


Imax = max(max(ImgOut));
ImgOut = ImgOut*(2/Imax)-1; %scale to [-1 1]

ImgReal = ImgOut;
ImgImaginary = sqrt(1 - ImgReal.^2);
ImgComplex = complex(ImgReal,ImgImaginary);
ImgOut = ImgComplex;

ImgNew = ImgOut;

maskindexCount = 0;
maskindexX = zeros(1,N*M);
maskindexY = zeros(1,N*M);

mask1 = zeros(N,M);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2+1:N-1+1,2  :M-1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2-1:N-1-1,2  :M-1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2  :N-1  ,2+1:M-1+1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2  :N-1  ,2-1:M-1-1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2+1:N-1+1,2+1:M-1+1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2-1:N-1-1,2-1:M-1-1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2+1:N-1+1,2-1:M-1-1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2-1:N-1-1,2+1:M-1+1);
mask1(2:N-1,2:M-1) = mask1(2:N-1,2:M-1) + mask0(2  :N-1  ,2  :M-1);

for y=2:N-1,
for x=2:M-1,
if mask1(y,x) > 0
	maskindexCount = maskindexCount+1;
	maskindexX(maskindexCount)=x;
	maskindexY(maskindexCount)=y;
end;end;end;


% begin iteration
for t=1:TIME,
	ImgOut = ImgNew;	

	% iterations
	for k = 1:maskindexCount,
		j = maskindexY(k);
		ImgOut = maskindexX(k);
		
		laplacianI = ImgOut(j,ImgOut-1)+ ImgOut(j,ImgOut+1)+ ImgOut(j-1,ImgOut)+ ImgOut(j+1,ImgOut)- (4+4/root2)*ImgOut(j,ImgOut);
		laplacianI = laplacianI+ ImgOut(j+1,ImgOut+1)/root2+ ImgOut(j-1,ImgOut-1)/root2+ ImgOut(j+1,ImgOut-1)/root2+ ImgOut(j-1,ImgOut+1)/root2;
		Isquare = real(ImgOut(j,ImgOut))^2+imag(ImgOut(j,ImgOut))^2;
		
		ImgNew(j,ImgOut) = ImgOut(j,ImgOut) + dt*(laplacianI) + (1-Isquare)*ImgOut(j,ImgOut)* dt/eps^2;
	end;	

	diffsum = max(max(abs(ImgNew-ImgOut)));	
	
	if diffsum < 10^-5
		fprintf('\rt=%d\n',t);
		t=TIME; 
		
		ImgOut = real(ImgNew);
		ImgOut = (ImgOut+1)*(Imax/2);
		return; 
	end;	

end;

ImgOut = real(ImgNew);
ImgOut = (ImgOut+1)*(Imax/2);
