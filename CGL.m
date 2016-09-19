function Img = CGL(Img0, mask0, TIME, dt, eps, h)

lambda=0.028; 
maskwidth = 2;

root2 = sqrt(2);

[N M]=size(Img0);

Img=double(Img0);
ImgLaplacian=zeros(N,M);


Imax = max(max(Img));
Img = Img*(2/Imax)-1; %scale to [-1 1]

ImgReal = Img;
ImgImaginary = sqrt(1 - ImgReal.^2);
ImgComplex = complex(ImgReal,ImgImaginary);
Img = ImgComplex;

ImgNew = Img;

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
	Img = ImgNew;	

	% iterations
	for k = 1:maskindexCount,
		j = maskindexY(k);
		Img = maskindexX(k);
		
		laplacianI = Img(j,Img-1)+ Img(j,Img+1)+ Img(j-1,Img)+ Img(j+1,Img)- (4+4/root2)*Img(j,Img);
		laplacianI = laplacianI+ Img(j+1,Img+1)/root2+ Img(j-1,Img-1)/root2+ Img(j+1,Img-1)/root2+ Img(j-1,Img+1)/root2;
		Isquare = real(Img(j,Img))^2+imag(Img(j,Img))^2;
		
		ImgNew(j,Img) = Img(j,Img) + dt*(laplacianI) + (1-Isquare)*Img(j,Img)* dt/eps^2;
	end;	

	diffsum = max(max(abs(ImgNew-Img)));	
	
	if diffsum < 10^-5
		fprintf('\rt=%d\n',t);
		t=TIME; 
		
		Img = real(ImgNew);
		Img = (Img+1)*(Imax/2);
		return; 
	end;	

end;

Img = real(ImgNew);
Img = (Img+1)*(Imax/2);

