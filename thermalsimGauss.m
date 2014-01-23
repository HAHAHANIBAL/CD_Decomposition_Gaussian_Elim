function [ Temperature ] = thermalsimGauss( p, mediumX, mediumY, leftBound, rightBound, topBound, bottomBound )
%calculate delta x and delta y
del_x=mediumX/length(topBound(:,1));
del_y=mediumY/length(rightBound(:,1));

%calculate size
Lx=length(topBound(:,1));
Ly=length(rightBound(:,1));

%constant conductivity
k=157;

%initialize new boundary array with known parameters
BoundaryArray=zeros(length(rightBound(:,1))+2,length(topBound(:,1))+2);
BoundaryArray(1,2:(end-1))=topBound;
BoundaryArray(2:(end-1),1)=leftBound;
BoundaryArray(2:(end-1),end)=rightBound;
BoundaryArray(end,2:(end-1))=bottomBound;

%initialize coefficient matrix A
A=zeros(Lx*Ly);
B=zeros(Lx*Ly,1);




%update coefficients matrix A (exclude boundary)
%m*n matrix, m is rows, n is columns
for m = 2:(Ly-1)
	for n = 2:(Lx-1)
		%Locate the elements 
		i = (n-1)*Ly + m;
		%Coefficient i+1,j
		A(i,i+Ly) = (1/(del_x^2))*k;
		%Coefficient i-1,j
		A(i,i-Ly) = (1/(del_x^2))*k;
		%Coefficient i,j+1
		A(i,i+1) = (1/(del_y^2))*k;
		%Coefficient i,j-1
		A(i,i-1) = (1/(del_y^2))*k;
		%Coefficient i,j
		A(i,i) = (-2/(del_y^2)-2/(del_x^2))*k;
		%Update Matrix B
		B(i)=-p(n,m);
	end
end

%update coeffcients matrix A (boundary)

%left and right bound
for m=2:(Ly-1)
	i=m;
	A(i,i+Ly) = (1/(del_x^2))*k;
	A(i,i+1) = (1/(del_y^2))*k;
	A(i,i-1) = (1/(del_y^2))*k;
	A(i,i) = (-2/(del_y^2)-2/(del_x^2))*k;
	B(i)=-p(1,m)-(k*BoundaryArray(end-m+1,1))/(del_x^2);
	j=m+(Lx-1)*Ly;
	A(j,j-Ly) = (1/(del_x^2))*k;
	A(j,j+1) = (1/(del_y^2))*k;
	A(j,j-1) = (1/(del_y^2))*k;
	A(j,j) = (-2/(del_y^2)-2/(del_x^2))*k;
	B(j)=-p(end,m)-(k*BoundaryArray(end-m+1,end))/(del_x^2);
end

%bottom and top bound
for n=2:(Lx-1)
	i=(n-1)*Ly+1;
	A(i,i+Ly) = (1/(del_x^2))*k;
	A(i,i-Ly) = (1/(del_x^2))*k;
	A(i,i+1) = (1/(del_y^2))*k;
	A(i,i) = (-2/(del_y^2)-2/(del_x^2))*k;
	B(i)=-p(n,end)-(k*BoundaryArray(end,end-n+1))/(del_y^2);
	j=n*Ly;
	A(j,j+Ly) = (1/(del_x^2))*k;
	A(j,j-Ly) = (1/(del_x^2))*k;
	A(j,j-1) = (1/(del_y^2))*k;
	A(j,j) = (-2/(del_y^2)-2/(del_x^2))*k;
	B(j)=-p(n,1)-(k*BoundaryArray(1,end-n+1))/(del_y^2);
end

%update coeffcients matrix A (corners)
%bottom left corner
A(1,1)=(-2/(del_y^2)-2/(del_x^2))*k;
A(1,2)=(1/(del_y^2))*k;
A(1,1+Ly)=(1/(del_x^2))*k;
B(1)=-p(1,1)-(k*BoundaryArray(end,2))/(del_y^2)-(k*BoundaryArray(end-1,1))/(del_x^2);

%Bottom right corner
A((Lx-1)*Ly+1,(Lx-1)*Ly+2)=(1/(del_y^2))*k;
A((Lx-1)*Ly+1,(Lx-1)*Ly+1-Ly)=(1/(del_x^2))*k;
A((Lx-1)*Ly+1,(Lx-1)*Ly+1)=(-2/(del_y^2)-2/(del_x^2))*k;
B((Lx-1)*Ly+1)=-p(end,1)-(k*BoundaryArray(end,end-1))/(del_y^2)-(k*BoundaryArray(end-1,end))/(del_x^2);

%top left corner
A(Ly,Ly+Ly)=(1/(del_x^2))*k;
A(Ly,Ly-1)=(1/(del_y^2))*k;
A(Ly,Ly)=(-2/(del_y^2)-2/(del_x^2))*k;
B(Ly)=-p(1,end)-(k*BoundaryArray(1,2))/(del_y^2)-(k*BoundaryArray(2,1))/(del_x^2);

%top right corner
A(Lx*Ly,Lx*Ly-1)=(1/(del_y^2))*k;
A(Lx*Ly,Lx*Ly)=(-2/(del_y^2)-2/(del_x^2))*k;
A(Lx*Ly,Lx*Ly-Ly)=(1/(del_x^2))*k;
B(Lx*Ly)=-p(end,end)-(k*BoundaryArray(1,end-1))/(del_y^2)-(k*BoundaryArray(2,end))/(del_x^2);



%answer=A\B;


%forward GE elimination
GE=cat(2,A,B);
rows=size(GE,1);
cols=size(GE,2);
T=zeros(rows,1);


for i=1:rows-1
	for j=i+1:rows		
		%initiate multiplier
		m=GE(j,i)/GE(i,i);
		GE(j,i:cols)=GE(j,i:cols)-m*GE(i,i:cols);
	end
end
	
%update B,A
B=GE(1:rows,rows+1);
A=GE(1:rows,1:rows);

%backward subsitution

%caculate the first T value
T(rows)=B(rows)/A(rows,rows);

%backward subsitute
for k=rows-1:-1:1
	T(k)=(B(k)-A(k,k+1:rows)*T(k+1:rows))/A(k,k);
end

%parse the data into p(i,j) format, plot the graph 
Temperature=zeros(size(p));
k=0;
for i=1:size(p,1)
for j=1:size(p,2)
Temperature(i,j)=T(k*size(p,2)+j);
end
k=k+1;
end
thermalplot(Temperature);



end
