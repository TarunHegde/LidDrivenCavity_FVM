%Lid Driven Cavity using Staggered Grid FVM, Explicit Method
%Input Parameters
Re = 100;                %Reynolds No.
dt = 0.001;              %Timestep
m = 129;                 %Grid Size
n = 129;
w = 1.9;                 %SOR Parameter
Length = 1;
Height = 1;
U_Lid = 1;
%Convergence Criteria
ErrPC = 0.00000001;
ErrVC = 0.000000001;
ErrUC = 0.000000001;

dx = 1/(m-1);
dy = 1/(n-1);
r = dx/dy;

U = zeros(m+1,n+1);
Uo = zeros(m+1,n+1);
V = zeros(m+1,n+1);
Vo = zeros(m+1,n+1);
P = zeros(m+1,n+1);
Po = zeros(m+1,n+1);
F = zeros(m+1,n+1);
G = zeros(m+1,n+1);
H = zeros(m+1,n+1);

count = 0;
t = 0;

while(1)
    count = count+1;
    t = t+dt;
    %Left and Right Walls
    for j=2:n
        Uo(1,j) = 0;
        Uo(m,j) = 0;
    end
    for j=1:n
        Vo(1,j) = -Vo(2,j);
        Vo(m+1,j) = -Vo(m,j);
    end
    %Top and Bottom Walls
    for i=2:m-1
        Uo(i,n+1) = 2 - Uo(i,n);
    end
    for i=1:m
        Uo(i,1) =  -Uo(i,2);
    end
    for i=2:m
        Vo(i,1) = 0;
        Vo(i,n) = 0;
    end
    %BC for Pressure
    for i=2:m
        Po(i,1) = Po(i,2);
        Po(i,n+1) = Po(i,n);
    end
    for j=2:n
        Po(1,j) = Po(2,j);
        Po(m+1,j) = Po(m,j);
    end
    %F
    for i=2:m-1
        for j=2:n
            Ue = (Uo(i,j) + Uo(i+1,j))/2;
            Uw = (Uo(i,j) + Uo(i-1,j))/2;
            Un = (Uo(i,j) + Uo(i,j+1))/2;
            Us = (Uo(i,j) + Uo(i,j-1))/2;
            Vn = (Vo(i,j) + Vo(i+1,j))/2;
            Vs = (Vo(i,j-1) + Vo(i+1,j-1))/2;
            dudx = (Uo(i+1,j) - 2*Uo(i,j) + Uo(i-1,j))/(dx*dx);
            dudy = (Uo(i,j+1) - 2*Uo(i,j) + Uo(i,j-1))/(dy*dy);
            F(i,j) = -((Ue*Ue - Uw*Uw)/dx)-((Un*Vn-Us*Vs)/dy)+(1/Re)*(dudx+dudy);
        end
    end
    %G
    for i=2:m
        for j=2:n-1
            Ue = (Uo(i,j) + Uo(i,j+1))/2;
            Uw = (Uo(i-1,j) + Uo(i-1,j+1))/2;
            Ve = (Vo(i,j) + Vo(i+1,j))/2;
            Vw = (Vo(i,j) + Vo(i-1,j))/2;
            Vn = (Vo(i,j) + Vo(i,j+1))/2;
            Vs = (Vo(i,j) + Vo(i,j-1))/2;
            dvdx = (Vo(i+1,j) - 2*Vo(i,j) + Vo(i-1,j))/(dx*dx);
            dvdy = (Vo(i,j+1) - 2*Vo(i,j) + Vo(i,j-1))/(dy*dy);
            G(i,j) = -((Ue*Ve - Uw*Vw)/dx) - ((Vn*Vn - Vs*Vs)/dy) + (1/Re)*(dvdx+dvdy);
        end
    end
    itrp = 0;
    while(1)
        itrp = itrp + 1;
        ErrorP = 0;
        Po = P;
        for i=2:m
            for j=2:n
                Aw = 1;
                Ae = 1;
                An = r*r;
                As = r*r;
                if i==2
                    Aw = 0;
                    F(i,j) = 0;
                end
                if i==m
                    Ae = 0;
                    F(i,j) = 0;
                end
                if j==2
                    As = 0;
                    G(i,j) = 0;
                end
                if j==n
                    An = 0;
                    G(i,j) = 0;
                end
                Ap = -(Aw+Ae+As+An);
                H(i,j) = dx*dx*((F(i,j)-F(i-1,j))/dx + (G(i,j) - G(i,j-1))/dy);
                P(i,j) = w*(H(i,j) - (Aw*P(i-1,j)+Ae*P(i+1,j)+As*P(i,j-1)+An*P(i,j+1)))/Ap + (1-w)*Po(i,j);
                Error = abs(P(i,j) - Po(i,j));
                if Error>ErrorP
                    ErrorP = Error;
                end
                Po(i,j) = P(i,j);
            end
        end
        if ErrorP<ErrPC
            break;
        end
    end
    ErrorV = 0;
    ErrorU = 0;
    %U
    for i=2:m-1
        for j=2:n
            U(i,j) = Uo(i,j) + dt*F(i,j) - (dt/dx)*(Po(i+1,j)-Po(i,j));
            Error = abs(U(i,j) - Uo(i,j));
            if Error>ErrorU
                ErrorU = Error;
            end
            Uo(i,j) = U(i,j);
        end
    end
    for i=2:m
        for j=2:n-1
            V(i,j) = Vo(i,j) + dt*G(i,j) - (dt/dy)*(Po(i,j+1)-Po(i,j));
            Error = abs(V(i,j) - Vo(i,j));
            if Error>ErrorV
                ErrorV = Error;
            end
            Vo(i,j) = V(i,j);
        end
    end
    %Left and Right Walls
    for j=2:n
        Uo(1,j) = 0;
        Uo(m,j) = 0;
    end
    for j=1:n
        Vo(1,j) = -Vo(2,j);
        Vo(m+1,j) = -Vo(m,j);
    end
    %Top and Bottom Walls
    for i=2:m-1
        Uo(i,n+1) = 2 - Uo(i,n);
    end
    for i=1:m
        Uo(i,1) = -Uo(i,2);
    end
    for i=2:m
        Vo(i,1) = 0;
        Vo(i,n) = 0;
    end
    %BC for Pressure
    for i=2:m
        Po(i,1) = Po(i,2);
        Po(i,n+1) = Po(i,n);
    end
    for j=2:n
        Po(1,j) = Po(2,j);
        Po(m+1,j) = Po(m,j);
    end
    if ErrorV<ErrVC
        if ErrorU<ErrUC
            break;
        end
    end
end

P_Grid = zeros(m,n);
U_Grid = zeros(m,n);
V_Grid = zeros(m,n);

for i=1:m
    for j=1:n
        U_Grid(i,j) = (Uo(i,j)+Uo(i,j+1))/2;
        V_Grid(i,j) = (Vo(i,j)+Vo(i+1,j))/2;
        P_Grid(i,j) = (Po(i,j)+Po(i+1,j)+Po(i,j+1)+Po(i+1,j+1))/4;
    end
end

P = P_Grid;
P = transpose(P);

figure(1);
clabel(contour([0:dx:1],[0:dy:1],P));
title('Isobars');

figure(2);
contourf([0:dx:1],[0:dy:1],P,'edgecolor','none');
title('Pressure Contours');
U_Grid = transpose(U_Grid);
V_Grid = transpose(V_Grid);

figure(3);
xlim([0 1]);
ylim([0 1]);
VEL = streamslice([0:dx:1],[0:dy:1],U_Grid,V_Grid,5);
set(VEL,'Color','k');