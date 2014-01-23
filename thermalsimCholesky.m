function [ Temperature ] = thermalsimCholesky( p, mediumX, mediumY, leftBound, rightBound, topBound, bottomBound )
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

%Cholesky

rows=size(A,1);
T=zeros(rows,1);
V=zeros(rows,1);
L=zeros(rows);

%-A is positive definite
A=-A;

%a=eig(A);
%L=chol(A);
%LT=L';

%conversion for A matrix into L lower(old method)

%first component
%A(1,1)=sqrt(A(1,1));

%for i=2:rows
%	A(i,1)=A(i,1)/A(1,1);
%end

%for i=2:rows-1
%	A(i,i)=sqrt(A(i,i)-sum(A(i,1:i-1).^2));
%	for j=i+1:rows
%		A(j,i)=(A(j,i)-sum(A(j,1:i-1).*A(i,1:i-1)))/A(i,i);
%	end
%end

%A(rows,rows)=sqrt(A(rows,rows)-sum(A(rows,1:rows-1).^2));

%for i=1:rows
%	for j=i+1:rows
%		A(i,j)=0;
%	end
%end
%L=A;

%faster L lower generation

for i=1:rows
	%Calculate diagonal values
	L(i,i)=sqrt(A(i,i)-sum(L(i,:).^2));
	for j=i+1:rows
		tmp=(A(j,i)-sum(L(j,:).*L(i,:)));
		L(j,i)=tmp/L(i,i);
	end
end

%calculate L transpose
LT=L';


%calculate first V value
V(1)=B(1)/L(1,1);

%forward subsitute
for k=2:1:rows
	V(k)=(B(k)-L(k,1:k-1)*V(1:k-1))/L(k,k);
end

%caculate the first T value
T(rows)=V(rows)/LT(rows,rows);

%backward subsitute
for k=rows-1:-1:1
	T(k)=(V(k)-LT(k,k+1:rows)*T(k+1:rows))/LT(k,k);
end



T=-T;

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
