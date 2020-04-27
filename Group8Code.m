%change to xlsread
%For example, A = sym('a',[1,3]) creates a row vector A = [a1,a2,a3].
% check for local stifness matrix
% take Iy and Iz alao from table
%% first we read the necessary data

% InM is member information matrix that has member no.(1) , frst_joint(2),
% sec_joint(3), Length(4), Area(5), depth(6), width(7), Ix(8), E(9), G(10)
M = input('Input Member file name:');
InfoMem = xlsread(M);      %change to xlsread
% Rearranging the 2nd and 3rd collumn to have smaller joint number in
% to have smaller joint no. in 2nd collumn
for l=1:size(InfoMem,1)
    if InfoMem(l,3) < InfoMem(l,2)
        temp = InfoMem(l,3);
        InfoMem(l,3) = InfoMem(l,2);
        InfoMem(l,2) = temp;
    end
end
k = input('Number of Joint restrained:( )'); 
% InJ is the Joint information Matrix that contain Co-ordinates, forces and momenrts of the joints.
J = input('Input Joint file name:');
InJ = xlsread(J);

%% Making the global stiffness matrix for each member

for l=1:size(InfoMem,1)
    %% Making the local stiffness matrix

    E = InfoMem(l,8);   % Elasticity
    % Ix and Iy have interchanged formula because of taking y in horizontal
    % direction
    G=InfoMem(l,9); J=InfoMem(l,6); Iy=InfoMem(l,7); Iz=InfoMem(l,8);%% take this also from table
    L = InfoMem(l,4);      % Length of member 
    A = InfoMem(l,5);      % Area of cross section
    % Allocating some variables to save them in memory and prevent calculation everytime
    U = A*E/L; V = 12*Iz*E/L^3; W = 12*Iy*E/L^3;  %%     dobara dekhna h Iy aur Iz bcoz Y taken horizontal while in notebook it is vertical
    X = 6*Iz*E/L^2; Y = 6*Iy*E/L^2; Z = G*J/L;
    % Taking all elements zero at first as max. elements are zeros
    Km(:,:,l) = zeros(12,12);
    % Putting the values in the upper triangular part
    Km(1,1,l) = U; Km(1,7,l) = -U; 
    Km(2,2,l) = V; Km(2,6,l) = X; Km(2,8,l) = -V; Km(2,12,l) = X;
    Km(3,3,l) = W; Km(3,5,l) = -Y; Km(3,9,l) = -W; Km(3,11,l) = -Y;
    Km(4,4,l) = Z; Km(4,10,l) = -Z;
    Km(5,5,l) = 4*E*Iy/L; Km(5,9,l) = Y; Km(5,12,l) = 2*E*Iy/L;
    Km(6,6,l) = 4*E*Iz/L; Km(6,8,l) = -X; Km(5,12,l) = 2*E*Iz/L;
    Km(7,7,l) = U;
    Km(8,8,l) = V; Km(8,12,l) = -X;
    Km(9,9,l) = W; Km(9,11,l) = Y;
    Km(10,10,l) = Z;
    Km(11,11,l) = 4*E*Iy/L;
    Km(12,12,l) = 4*E*Iz/L;
    % Making the symetric part by transforming processes
    P = Km(:,:,l) - Km(:,:,l).*eye(12);
    Km(:,:,l) = Km(:,:,l) + P';
    
    %% Making the transformation matrix
    % Extracting near joint as nj and far joint fj
    nj = InfoMem(l,2); fj = InfoMem(l,3);
    % Finding the axial x-vector of local member
    VectorX = [InJ(fj,2)-InJ(nj,2) InJ(fj,3)-InJ(nj,3) InJ(fj,4)-InJ(nj,4)];
    % Converting to unit matrix
    VectorX = VectorX./norm(VectorX);
    % Saving cos of thetas as xlambdas in X, Y, Z directions
    xlambda = VectorX';
    % if x in direction of X then dt is identity 
    if xlambda(1) == 1 
        dt = eye(12);
        % if x opp. of X then dt as follows and similarly others
    elseif xlambda(1) == -1
        D = [-1 0 0;0 -1 0;0 0 1];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == 1
        D = [0 1 0;-1 0 0;0 0 1];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == -1
        D = [0 -1 0;1 0 0;0 0 1];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(3) == 1
        D = [0 0 1;0 1 0;-1 0 0];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(3) == -1
        D = [0 0 -1;0 1 0;1 0 0];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(1) == 0 && xlambda(2) > 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; -1 0 0; 0 -xlambda(3) xlambda(2)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(1) == 0 && xlambda(2) < 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; 1 0 0; 0 xlambda(3) -xlambda(2)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == 0 && xlambda(1) > 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; 0 1 0;-xlambda(3) 0 xlambda(1)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == 0 && xlambda(1) < 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; 0 -1 0;xlambda(3) 0 -xlambda(1)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    end
    Dt(:,:,l) = dt;
    % The global matrix of each member is ready by mulltiplying with dt in
    % following manner
    K(:,:,l) = Dt(:,:,l)'*Km(:,:,l)*Dt(:,:,l);
    
end
%% Making the global stifness matrix for arranged known and unknown joint number
%Again pre-allocating as in local stiffness matrix
Kg = zeros(6*size(InJ,1),6*size(InJ,1));
%% 
%adding all the small global matrices to find the bigger global matrices
for j=1:size(InJ,1)         %add condition of being equal joint make loop for all joint combination
    for b=1:size(InJ,1)
        if j==b
            i = find(InfoMem(:,2) == j);
            for a = 1:size(i)
                Kg(j*6-5:j*6,j*6-5:j*6) = Kg(j*6-5:j*6,j*6-5:j*6) + K(1:6,1:6,i(a));
            end
            i = find(InfoMem(:,3) == j);
            for a = 1:size(i)
                Kg(j*6-5:j*6,j*6-5:j*6) = Kg(j*6-5:j*6,j*6-5:j*6) + K(7:12,7:12,i(a));
            end
        elseif j > b
            i = find(InfoMem(:,2) == b);
            f = find(InfoMem(:,3) == j);
            for q = 1:size(i)
                for m = 1:size(f)
                    if i(q) == f(m)                        
                        Kg(j*6-5:j*6,b*6-5:b*6) = K(7:12,1:6,i(q));                        
                        Kg(b*6-5:b*6,j*6-5:j*6) = K(1:6,7:12,i(q));
                    end
                end
            end
        end
    end
end
%%
% Getting all the forces 
z = InJ(k+1:end,5:10)';
%converting the forces into vectors
Pc = z(:);
%getting unknown displacements by using inverse
dis = (Kg(6*k+1:end,6*k+1:end))\Pc;
%net displacement vector
d = [zeros(6*k,1); dis(:)];
% x,y,z disp and rotations
Xd = reshape(d, 6, size(InJ,1))';
%net forces and reaction vector
P = Kg*d;   
% the values of all forces and reaction in order
fr = reshape(P,6,size(InJ,1))';

%finding local end moments and forces
for l=1:size(InfoMem,1)
     % Extracting near joint as nj and far joint fj
    nj = InfoMem(l,2); fj = InfoMem(l,3);
    as = [Xd(nj,:) Xd(fj,:)]';
    fx(:,l) = K(:,:,l)*Dt(:,:,l)*as(:);
